#!/usr/bin/env python
"""Stitch every well of a phenix-icechunk plate with multiview-stitcher.

For each well:

1. Read: every field comes from the repo as an ``xr.DataArray`` (``open_datatree``). The
   requested z projection is used if the repo stores it (``--zproject``), otherwise it is computed
   here; ``--projection none`` keeps the 3D stacks. Flatfield follows the repo default unless set.
2. Cache: the preprocessed fields are written once to a chunked zarr (``OUTPUT_DIR/_tiles``), so
   registration and fusion read small chunks instead of whole source planes, each only once.
3. Register: fields start at their stage positions; overlapping neighbours are aligned by phase
   correlation and resolved globally. Pairs below ``--quality-threshold`` are dropped, every group
   of registered fields is centred on the stage frame, and a field whose shift disagrees with its
   neighbours' by more than ``--outlier-um`` follows them instead. Wells whose fields don't overlap
   are placed by stage position. 3D stacks are registered in y/x on their max projection.
4. Fuse and write:

   * ``<plate>_<well>[_z<method>].ome.zarr``: one OME-Zarr per well with a resolution pyramid,
     channel names/colours and the per-field shifts in ``attrs["phenix_stitch"]``
   * ``<plate>[_z<method>]_stitched.icechunk``: the same images as an OME-Zarr 0.5 HCS plate
     (``/<row>/<col>/0``), one icechunk commit per well
   * ``<plate>[_z<method>]_field_shifts.csv``: per-field shifts of every stitched well
   * ``<plate>_<well>[_z<method>]_qc.png``: shift map and channel previews

A well is written to ``<name>.partial`` and renamed when complete. Without ``--overwrite``,
finished wells are skipped, so an interrupted run resumes where it stopped.

Needs multiview-stitcher and matplotlib next to phenix-icechunk::

    uv pip install --python .venv/bin/python multiview-stitcher matplotlib

Example::

    .venv/bin/python scripts/stitch_wells.py data/plate.icechunk /path/to/stitched             # every well
    .venv/bin/python scripts/stitch_wells.py data/plate.icechunk /path/to/stitched --well A4   # one well
"""

from __future__ import annotations

import argparse
import copy
import logging
import re
import shutil
import time
import warnings
from dataclasses import dataclass, fields
from importlib.metadata import version
from pathlib import Path
from typing import Any

import dask
import dask.array as da
import icechunk as ic
import matplotlib
import networkx as nx
import numpy as np
import pandas as pd
import xarray as xr
import zarr
from dask.diagnostics import ProgressBar
from icechunk.xarray import to_icechunk
from multiview_stitcher import fusion, msi_utils, ngff_utils, param_utils, registration
from multiview_stitcher import spatial_image_utils as si_utils

from phenix_icechunk import open_datatree, open_plate
from phenix_icechunk.convert import open_repo

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

log = logging.getLogger("stitch_wells")

PROJECTIONS = ("max", "min", "mean", "median", "sum", "std")
YX = {"x_in": ["y", "x"], "x_out": "1"}  # translation part of a 2D affine


@dataclass
class Options:
    repo: str
    output_dir: Path
    wells: list[str] | None = None  # None = every well
    projection: str | None = "max"  # stored one if the repo has it, else computed; None = 3D
    flatfield: bool | None = None  # None = repo default
    output_dtype: str = "uint16"
    register: bool = True
    reg_channel: str | None = None  # None = first channel
    registration_binning: int = 2
    quality_threshold: float = 0.2
    min_overlap_um: float = 5.0
    outlier_um: float = 10.0
    fix_outliers: bool = True
    chunk: int = 2048
    ngff_version: str = "0.5"
    workers: int = 16
    keep_tile_cache: bool = False
    overwrite: bool = False
    icechunk: Path | None = None  # None = OUTPUT_DIR/<plate>[_z<method>]_stitched.icechunk


DEFAULTS = {f.name: f.default for f in fields(Options)}


def well_path(well: str) -> str:
    # "A4" -> "A/4"
    if "/" in well:
        return well
    row = well.rstrip("0123456789")
    return f"{row}/{well[len(row):]}"


def spacing(img: xr.DataArray, dim: str) -> float:
    return float(img[dim][1] - img[dim][0]) if img.sizes[dim] > 1 else 1.0


def stage_overlaps(origin: np.ndarray, size: np.ndarray) -> np.ndarray:
    # Overlap (µm, along the narrower axis) of every pair of fields whose stage boxes intersect in y and x
    hi = origin + size
    overlap = (np.minimum(hi[:, None], hi[None]) - np.maximum(origin[:, None], origin[None])).min(-1)
    i, j = np.triu_indices(len(origin), k=1)
    overlap = overlap[i, j]
    return overlap[overlap > 0]


def neighbours(origin: np.ndarray, size: np.ndarray) -> list[np.ndarray]:
    # Each field's stage-grid neighbours: origins less than 1.5 field sizes apart in y and x
    near = (np.abs(origin[:, None] - origin[None]) < 1.5 * size).all(-1)
    np.fill_diagonal(near, False)
    return [np.flatnonzero(row) for row in near]


def ome_metadata(group: zarr.Group) -> tuple[dict[str, Any], dict[str, Any] | None]:
    # (multiscale, omero) of an OME-Zarr image group, NGFF 0.4 or 0.5
    attrs = group.attrs.asdict()
    ome = attrs.get("ome", attrs)
    multiscale = {k: v for k, v in ome["multiscales"][0].items() if k != "version"}
    return multiscale, ome.get("omero")


class PlateStitcher:
    def __init__(self, opts: Options):
        self.opts = opts
        self.plate = open_plate(opts.repo)
        self.meta = self.plate.attrs["phenix"]
        projections = self.meta.get("projections", {})
        self.stored = bool(opts.projection) and f"z{opts.projection}" in projections
        if self.stored and opts.flatfield is not None:
            applied = projections[f"z{opts.projection}"]["flatfield_applied"]
            if opts.flatfield != applied:
                raise ValueError(
                    f"stored z{opts.projection} was computed with flatfield={applied}; leave --flatfield unset"
                )
        self.computed = bool(opts.projection) and not self.stored
        self.tree = open_datatree(
            self.plate,
            zproject=opts.projection if self.stored else None,
            flatfield=None if self.stored else opts.flatfield,
            include_table=False,
        )
        if self.stored:
            self.source = f"stored z{opts.projection} projection"
            self.flatfield = projections[f"z{opts.projection}"]["flatfield_applied"]
        else:
            self.source = f"z{opts.projection} computed here" if opts.projection else "3D stacks"
            self.flatfield = None  # known per field (attrs["flatfield"])
        self.source_snapshot = open_repo(opts.repo).lookup_branch("main")

        self.plate_name = re.sub(
            r"[^A-Za-z0-9._-]+", "_", self.meta.get("plate_name") or self.meta.get("plate_id") or "plate"
        )
        self.suffix = f"_z{opts.projection}" if opts.projection else ""
        all_wells = [w["path"] for w in self.plate.attrs["ome"]["plate"]["wells"]]
        self.wells = all_wells if opts.wells is None else [well_path(w) for w in opts.wells]
        missing = sorted(set(self.wells) - set(all_wells))
        if missing:
            raise KeyError(f"wells {missing} not in the plate; available: {all_wells}")
        self.icechunk_path = opts.icechunk or opts.output_dir / f"{self.plate_name}{self.suffix}_stitched.icechunk"
        self.shifts_csv = opts.output_dir / f"{self.plate_name}{self.suffix}_field_shifts.csv"

    # ---- reading -------------------------------------------------------------------------------

    def well_fields(self, well: str) -> list[xr.DataArray]:
        # Raw fields of one well, dask-backed: one chunk per plane, or per z stack when projecting here.
        # The explicit token skips hashing each lazily indexed array (~80 ms per field otherwise).
        node = self.tree[well]
        z_chunk = -1 if self.computed else 1
        return [
            node[key]["image"].chunk(
                {"t": 1, "c": 1, "z": z_chunk},
                token=f"{self.opts.repo}:{well}/{key}:{self.source}:{self.opts.flatfield}",
            )
            for key in sorted(node.children, key=int)
        ]

    def is_3d(self, img: xr.DataArray) -> bool:
        return not self.opts.projection and img.sizes["z"] > 1

    def preprocess(self, img: xr.DataArray) -> xr.DataArray:
        # Project (if computed here) or drop a single z plane, then cast. Lazy.
        if self.computed:
            img = getattr(img, self.opts.projection)("z")
        elif not self.is_3d(img):
            img = img.isel(z=0, drop=True)
        if self.opts.output_dtype == "uint16" and img.dtype != np.uint16:
            return img.round().clip(0, 65535).astype(np.uint16)
        return img.astype(self.opts.output_dtype)

    def cache_tiles(self, raw: list[xr.DataArray], path: Path) -> xr.DataArray:
        # All fields of a well as one (field, t, c, [z,] y, x) array on disk. Stage origins become
        # per-field coords; y/x are split into equal tiles of at most 600 px (2160 -> 540).
        images = []
        for img in raw:
            img = self.preprocess(img).reset_coords(drop=True)
            images.append(
                img.drop_vars(["y", "x"]).assign_coords(origin_y=float(img.y[0]), origin_x=float(img.x[0]))
            )
        tiles = xr.concat(
            images, dim="field", coords=["origin_y", "origin_x"], compat="override", join="override",
            combine_attrs="drop",
        )
        tile = {d: -(-n // -(-n // 600)) for d, n in (("y", tiles.sizes["y"]), ("x", tiles.sizes["x"]))}
        tiles = tiles.chunk({"field": 1, "t": 1, "c": 1, "z": 1, **tile} if "z" in tiles.dims else
                            {"field": 1, "t": 1, "c": 1, **tile})
        tiles.attrs = {"spacing_y": spacing(raw[0], "y"), "spacing_x": spacing(raw[0], "x")}
        with ProgressBar(dt=10):
            tiles.to_dataset(name="tiles").to_zarr(path, mode="w", zarr_format=3, consolidated=False)
        return xr.open_zarr(path, consolidated=False)["tiles"]

    @staticmethod
    def to_msims(tiles: xr.DataArray) -> list:
        sdims = [d for d in ("z", "y", "x") if d in tiles.dims]
        msims = []
        for i in range(tiles.sizes["field"]):
            tile = tiles.isel(field=i)
            scale = {"y": tiles.attrs["spacing_y"], "x": tiles.attrs["spacing_x"]}
            origin = {"y": float(tile.origin_y), "x": float(tile.origin_x)}
            if "z" in sdims:
                scale["z"], origin["z"] = spacing(tiles, "z"), float(tiles.z[0])
            sim = si_utils.get_sim_from_array(
                tile.transpose("t", "c", *sdims).data,
                dims=["t", "c", *sdims],
                scale=scale,
                translation=origin,
                transform_key="stage",
                c_coords=[str(c) for c in tiles.c.values],
                t_coords=tiles.t.values,
            )
            msims.append(msi_utils.get_msim_from_sim(sim, scale_factors=[]))
        return msims

    # ---- registration --------------------------------------------------------------------------

    def register(self, msims: list, tiles: xr.DataArray, reg_channel: str) -> tuple[xr.Dataset, dict]:
        # y/x registration of 2D fields; sets transform "registered" and returns per-field shifts and QC
        o = self.opts
        result = registration.register(
            msims,
            transform_key="stage",
            new_transform_key="registered",
            reg_channel=reg_channel,
            registration_binning={"y": o.registration_binning, "x": o.registration_binning},
            pre_registration_pruning_method="keep_axis_aligned",
            post_registration_do_quality_filter=True,
            post_registration_quality_threshold=o.quality_threshold,
            n_parallel_pairwise_regs=o.workers,
            return_dict=True,
        )
        qualities = {
            pair: float(np.min(q)) for pair, q in result["pairwise_registration"]["metrics"]["qualities"].items()
        }
        good = nx.Graph()
        good.add_nodes_from(range(len(msims)))
        good.add_edges_from(pair for pair, q in qualities.items() if q >= o.quality_threshold)

        affines = [
            msi_utils.get_transform_from_msim(m, "registered").transpose("t", "x_in", "x_out").copy() for m in msims
        ]
        shift = xr.concat([a.sel(YX, drop=True) for a in affines], dim="field")  # (field, t, x_in)

        # Each connected group of fields is registered relative to an arbitrary field of the group.
        # Centre every group on zero shift so groups and unregistered fields share the stage frame.
        groups = list(nx.connected_components(good))
        for group in groups:
            idx = sorted(group)
            shift[idx] = shift[idx] - shift[idx].median("field")

        # A field whose shift differs from its neighbours' median by more than outlier_um is most
        # likely misregistered (little structure in its overlaps); fix_outliers moves it with them.
        origin, size = self.stage_boxes(tiles)
        near = neighbours(origin, size)
        expected = xr.concat(
            [shift[n].median("field") if len(n) else shift[i] for i, n in enumerate(near)], dim="field"
        )
        deviation = np.sqrt(((shift - expected) ** 2).sum("x_in")).max("t")
        outlier = deviation > o.outlier_um
        if o.fix_outliers:
            shift = xr.where(outlier, expected, shift)

        for msim, affine, s in zip(msims, affines, shift.transpose("field", "t", "x_in")):
            affine.loc[YX] = s.values
            msi_utils.set_affine_transform(msim, affine, transform_key="registered")

        qc = {
            "pairs": len(qualities),
            "good_pairs": good.number_of_edges(),
            "median_quality": float(np.median(list(qualities.values()))) if qualities else np.nan,
            "fields_at_stage_position": sum(1 for g in groups if len(g) == 1),
            "outliers": int(outlier.sum()),
        }
        shifts = xr.Dataset(
            {
                "shift_y_um": shift.sel(x_in="y", drop=True),
                "shift_x_um": shift.sel(x_in="x", drop=True),
                "deviation_um": deviation,
                "outlier": outlier,
            }
        )
        return shifts, qc

    @staticmethod
    def stage_boxes(tiles: xr.DataArray) -> tuple[np.ndarray, np.ndarray]:
        origin = np.stack([tiles.origin_y.values, tiles.origin_x.values], axis=-1)
        size = np.array([tiles.sizes["y"] * tiles.attrs["spacing_y"], tiles.sizes["x"] * tiles.attrs["spacing_x"]])
        return origin, size

    @staticmethod
    def lift_to_3d(msims_2d: list, msims_3d: list) -> None:
        # Copy the y/x registration of the projections onto the 3D stacks (z untouched: fields share their planes)
        yx1 = {"x_in": ["y", "x", "1"], "x_out": ["y", "x", "1"]}
        for m2, m3 in zip(msims_2d, msims_3d):
            affine2 = msi_utils.get_transform_from_msim(m2, "registered").transpose("t", "x_in", "x_out")
            affine3 = param_utils.identity_transform(3, t_coords=affine2.t.values).transpose("t", "x_in", "x_out")
            affine3.loc[yx1] = affine2.sel(yx1).values
            msi_utils.set_affine_transform(m3, affine3, transform_key="registered")

    # ---- one well ------------------------------------------------------------------------------

    def output_path(self, well: str) -> Path:
        return self.opts.output_dir / f"{self.plate_name}_{well.replace('/', '')}{self.suffix}.ome.zarr"

    def stitch_well(self, well: str, reg_channel: str) -> tuple[dict, xr.Dataset]:
        o = self.opts
        out_path = self.output_path(well)
        name = out_path.name.removesuffix(".ome.zarr")
        raw = self.well_fields(well)
        three_d = self.is_3d(raw[0])
        timings = {}

        t0 = time.perf_counter()
        cache = o.output_dir / "_tiles" / f"{name}.zarr"
        cache.parent.mkdir(parents=True, exist_ok=True)
        log.info("%s: caching %d fields (%s) -> %s", well, len(raw), self.source, cache)
        tiles = self.cache_tiles(raw, cache)
        msims = self.to_msims(tiles)
        timings["cache_s"] = time.perf_counter() - t0

        t0 = time.perf_counter()
        overlaps = stage_overlaps(*self.stage_boxes(tiles))
        do_register = o.register and len(raw) > 1 and bool((overlaps >= o.min_overlap_um).any())
        if do_register:
            log.info("%s: registering on %s%s", well, reg_channel, " (max projection of the stacks)" if three_d else "")
            if three_d:
                reg_msims = self.to_msims(tiles.max("z", keep_attrs=True))
                shifts, qc = self.register(reg_msims, tiles, reg_channel)
                self.lift_to_3d(reg_msims, msims)
            else:
                shifts, qc = self.register(msims, tiles, reg_channel)
            transform_key = "registered"
        else:
            reason = "--no-register" if not o.register else f"no field overlap >= {o.min_overlap_um} µm"
            log.info("%s: not registering (%s), fields stay at their stage positions", well, reason)
            zeros = xr.zeros_like(tiles.origin_y.expand_dims(t=tiles.t.values), dtype=float).reset_coords(drop=True)
            shifts = xr.Dataset({
                "shift_y_um": zeros, "shift_x_um": zeros,
                "deviation_um": zeros.isel(t=0, drop=True), "outlier": zeros.isel(t=0, drop=True).astype(bool),
            })
            qc = {"pairs": 0, "good_pairs": 0, "median_quality": np.nan, "fields_at_stage_position": len(raw),
                  "outliers": 0}
            transform_key = "stage"
        shifts = shifts.assign_coords(
            field=np.arange(len(raw)), well=well,
            stage_y_um=("field", tiles.origin_y.values), stage_x_um=("field", tiles.origin_x.values),
        )
        timings["register_s"] = time.perf_counter() - t0

        t0 = time.perf_counter()
        partial = out_path.with_name(out_path.name + ".partial")
        chunks = {"y": o.chunk, "x": o.chunk} | ({"z": tiles.sizes["z"]} if three_d else {})
        fused = fusion.fuse(
            images=[msi_utils.get_sim_from_msim(m) for m in msims], transform_key=transform_key, output_chunksize=chunks
        )
        log.info("%s: fusing %s -> %s", well, dict(fused.sizes), out_path)
        ngff_utils.write_sim_to_ome_zarr(
            fused, str(partial), overwrite=True, ngff_version=o.ngff_version, batch_options={"n_batch": None}
        )
        first = next(iter(self.tree[well].children.values()))["image"]
        record = {
            "source_repo": str(Path(o.repo).resolve()),
            "source_snapshot": self.source_snapshot,
            "well": well,
            "fields": len(raw),
            "pixels": self.source,
            "flatfield": self.flatfield if self.stored else "flatfield" in first.attrs,
            "transform_key": transform_key,
            "registration": {
                "channel": reg_channel, "binning": o.registration_binning, "quality_threshold": o.quality_threshold,
                "outlier_um": o.outlier_um, "fix_outliers": o.fix_outliers,
            } | qc if do_register else None,
            "field_shifts_um": shifts.to_dataframe().reset_index().drop(columns="well").to_dict(orient="records"),
            "multiview_stitcher": version("multiview-stitcher"),
            "phenix_icechunk": version("phenix-icechunk"),
        }
        group = zarr.open_group(str(partial), mode="r+")
        attrs = group.attrs.asdict()
        omero = first.attrs.get("ome", {}).get("omero")
        if omero:
            (attrs["ome"] if o.ngff_version == "0.5" else attrs)["omero"] = omero
        attrs["phenix_stitch"] = record
        group.attrs.update(attrs)
        if out_path.exists():
            shutil.rmtree(out_path)
        partial.rename(out_path)
        timings["fuse_s"] = time.perf_counter() - t0

        if not o.keep_tile_cache:
            shutil.rmtree(cache)
        return {"well": well, "registered": do_register, **qc, **timings}, shifts

    # ---- icechunk ------------------------------------------------------------------------------

    def open_icechunk(self) -> ic.Repository:
        storage = ic.local_filesystem_storage(str(self.icechunk_path))
        if ic.Repository.exists(storage):
            return ic.Repository.open(storage)
        return ic.Repository.create(storage)

    def icechunk_wells(self) -> set[str]:
        if not ic.Repository.exists(ic.local_filesystem_storage(str(self.icechunk_path))):
            return set()
        session = self.open_icechunk().readonly_session(branch="main")
        try:
            root = zarr.open_group(session.store, mode="r", zarr_format=3)
        except FileNotFoundError:
            return set()
        return {w["path"] for w in root.attrs.get("ome", {}).get("plate", {}).get("wells", [])}

    def write_icechunk(self, well: str) -> str:
        # Copy one stitched well into the plate repo (/<row>/<col>/0, OME-Zarr 0.5) as one commit
        t0 = time.perf_counter()
        src = zarr.open_group(str(self.output_path(well)), mode="r")
        multiscale, omero = ome_metadata(src)
        record = src.attrs["phenix_stitch"]
        axes = [a["name"] for a in multiscale["axes"]]

        repo = self.open_icechunk()
        session = repo.writable_session("main")
        root = zarr.open_group(session.store, mode="a", zarr_format=3)
        if well in root:
            del root[well]
        image = f"{well}/0"
        root.require_group(well)
        root.require_group(image)
        for dataset in multiscale["datasets"]:
            level = xr.DataArray(da.from_zarr(src[dataset["path"]]), dims=axes, name=dataset["path"])
            to_icechunk(level.to_dataset(), session, group=image, mode="a")

        # attrs last: xarray resets the attrs of the group it writes to
        ome = {"version": "0.5", "multiscales": [multiscale]} | ({"omero": omero} if omero else {})
        root[image].attrs.put({"ome": ome, "phenix_stitch": record})
        root[well].attrs.put({
            "ome": {"version": "0.5", "well": {"images": [{"path": "0", "acquisition": 0}]}},
            "phenix": dict(self.plate[well].attrs.get("phenix", {})),
        })
        done = {w["path"] for w in root.attrs.get("ome", {}).get("plate", {}).get("wells", [])} | {well}
        plate = copy.deepcopy(self.plate.attrs["ome"]["plate"])
        plate["name"] = f"{plate.get('name', self.plate_name)} stitched{self.suffix}"
        plate["wells"] = [w for w in plate["wells"] if w["path"] in done]
        plate["field_count"] = 1
        root.attrs.put({
            "ome": {"version": "0.5", "plate": plate},
            "phenix_stitch": {
                "source_repo": record["source_repo"],
                "pixels": record["pixels"],
                "wells": sorted(done),
                "multiview_stitcher": record["multiview_stitcher"],
            },
        })
        message = f"Stitched well {well} ({record['fields']} fields, {record['pixels']}, {record['transform_key']})"
        metadata = {"well": well} | ({"source_snapshot": record["source_snapshot"]} if "source_snapshot" in record else {})
        snapshot = session.commit(message, metadata=metadata)
        log.info("%s: committed %s to %s in %.0f s", well, snapshot, self.icechunk_path, time.perf_counter() - t0)
        return snapshot

    # ---- QC ------------------------------------------------------------------------------------

    def plot_qc(self, well: str) -> Path:
        # Shift map (arrows ×10, outliers circled) next to every channel of a ~1000 px pyramid level
        out_path = self.output_path(well)
        group = zarr.open_group(str(out_path), mode="r")
        multiscale, omero = ome_metadata(group)
        record = group.attrs["phenix_stitch"]
        levels = [d["path"] for d in multiscale["datasets"]]
        level = next((p for p in reversed(levels) if group[p].shape[-1] >= 1000), levels[0])
        axes = [a["name"] for a in multiscale["axes"]]
        img = xr.DataArray(group[level][:], dims=axes).isel(t=0)
        if "z" in img.dims:
            img = img.isel(z=img.sizes["z"] // 2)
        names = [ch.get("label", str(i)) for i, ch in enumerate((omero or {}).get("channels", []))]

        n = img.sizes["c"] + 1
        aspect = img.sizes["y"] / img.sizes["x"]
        fig, axs = plt.subplots(1, n, figsize=(5 * n, 5 * max(aspect, 0.6)), squeeze=False)
        ax = axs[0, 0]
        s = pd.DataFrame(record["field_shifts_um"])
        s = s[s["t"] == s["t"].iloc[0]]
        length = np.hypot(s["shift_y_um"], s["shift_x_um"])
        if record["transform_key"] == "registered":
            q = ax.quiver(s["stage_x_um"], s["stage_y_um"], s["shift_x_um"], s["shift_y_um"], length,
                          angles="xy", scale_units="xy", scale=0.1, cmap="viridis")
            fig.colorbar(q, ax=ax, shrink=0.7, label="shift (µm)")
            out = s[s["outlier"]]
            ax.scatter(out["stage_x_um"], out["stage_y_um"], s=120, facecolors="none", edgecolors="r")
            ax.set_title(f"{well}: field shifts ×10")
        else:
            ax.set_title(f"{well}: stage positions (not registered)")
        ax.scatter(s["stage_x_um"], s["stage_y_um"], s=4, c="k")
        ax.set(aspect="equal", xlabel="x (µm)", ylabel="y (µm)")
        ax.invert_yaxis()
        for ax, c in zip(axs[0, 1:], range(img.sizes["c"])):
            plane = img.isel(c=c).values
            nonzero = plane[plane > 0]
            lo, hi = np.percentile(nonzero, [1, 99.5]) if nonzero.size else (0, 1)
            ax.imshow(plane, cmap="gray", vmin=lo, vmax=hi)
            ax.set_title(names[c] if c < len(names) else f"channel {c}")
            ax.axis("off")
        fig.suptitle(f"{out_path.name} (level {level}, {img.sizes['y']}×{img.sizes['x']})")
        png = out_path.with_name(out_path.name.removesuffix(".ome.zarr") + "_qc.png")
        fig.savefig(png, dpi=100, bbox_inches="tight")
        plt.close(fig)
        return png

    # ---- plate ---------------------------------------------------------------------------------

    def run(self) -> pd.DataFrame:
        o = self.opts
        o.output_dir.mkdir(parents=True, exist_ok=True)
        log.info("%s: %d of %d wells, %s -> %s", self.plate_name, len(self.wells),
                 len(self.plate.attrs["ome"]["plate"]["wells"]), self.source, o.output_dir)

        channels = None
        for well in self.wells:
            raw = self.well_fields(well)
            origin = np.array([[float(img.y[0]), float(img.x[0])] for img in raw])
            size = np.array([raw[0].sizes["y"] * spacing(raw[0], "y"), raw[0].sizes["x"] * spacing(raw[0], "x")])
            overlaps = stage_overlaps(origin, size)
            dims = {d: raw[0].sizes[d] for d in ("t", "c", "z", "y", "x") if d != "z" or self.is_3d(raw[0])}
            log.info(
                "%s: %d fields %s, %.4f µm/px, field %.1f µm, %d overlapping pairs (median %.1f µm)",
                well, len(raw), dims, spacing(raw[0], "x"), size[1],
                int((overlaps >= o.min_overlap_um).sum()), float(np.median(overlaps)) if len(overlaps) else 0.0,
            )
            channels = channels or [str(c) for c in raw[0].c.values]
        reg_channel = o.reg_channel or channels[0]
        if reg_channel not in channels:
            raise KeyError(f"--reg-channel {reg_channel!r} not in {channels}")
        log.info("channels %s, registering on %s", channels, reg_channel if o.register else "nothing (--no-register)")

        in_icechunk = self.icechunk_wells()
        results, all_shifts = [], []
        for well in self.wells:
            out_path = self.output_path(well)
            stitched = out_path.exists() and not o.overwrite
            if stitched and well in in_icechunk:
                log.info("%s: %s and the icechunk copy exist, skipped (--overwrite to redo)", well, out_path.name)
                continue
            t0 = time.perf_counter()
            result = {"well": well}
            if stitched:
                log.info("%s: %s exists, only adding it to the icechunk repo", well, out_path.name)
            else:
                result, shifts = self.stitch_well(well, reg_channel)
                all_shifts.append(shifts)
            result["snapshot"] = self.write_icechunk(well)
            result["qc_png"] = str(self.plot_qc(well))
            result["total_s"] = time.perf_counter() - t0
            log.info("%s: done in %.0f s", well, result["total_s"])
            results.append(result)

        if all_shifts:
            table = pd.concat([s.to_dataframe().reset_index() for s in all_shifts], ignore_index=True)
            if self.shifts_csv.exists():
                old = pd.read_csv(self.shifts_csv)
                table = pd.concat([old[~old["well"].isin(table["well"])], table], ignore_index=True)
            table.to_csv(self.shifts_csv, index=False)
            log.info("per-field shifts -> %s", self.shifts_csv)
        tiles_dir = o.output_dir / "_tiles"
        if tiles_dir.exists() and not any(tiles_dir.iterdir()):
            tiles_dir.rmdir()
        return pd.DataFrame(results)


def parse_args(argv: list[str] | None = None) -> Options:
    d = DEFAULTS
    p = argparse.ArgumentParser(
        description="Stitch every well of a phenix-icechunk plate with multiview-stitcher into OME-Zarr "
        "and an icechunk HCS plate.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("repo", help="icechunk repo written by phenix-icechunk")
    p.add_argument("output_dir", type=Path, help="gets one <plate>_<well>[_z<method>].ome.zarr per well, "
                   "the stitched icechunk repo, the shifts CSV and QC plots")
    p.add_argument("--wells", "--well", nargs="+", default=d["wells"], metavar="WELL",
                   help="only stitch these wells, e.g. --well A4 for one well or --wells A4 B/3; "
                   "every well if not given")
    p.add_argument("--projection", choices=[*PROJECTIONS, "none"], default=d["projection"],
                   help="z projection to stitch (the stored one if the repo has it, else computed); none = 3D stacks")
    p.add_argument("--flatfield", action=argparse.BooleanOptionalAction, default=d["flatfield"],
                   help="flatfield correction; the repo default if not given (a stored projection keeps its own)")
    p.add_argument("--output-dtype", choices=("uint16", "float32"), default=d["output_dtype"],
                   help="uint16 is rounded and clipped")
    p.add_argument("--register", action=argparse.BooleanOptionalAction, default=d["register"],
                   help="register overlapping fields; --no-register places them at their stage positions")
    p.add_argument("--reg-channel", default=d["reg_channel"], help="channel to register on; the first channel if not given")
    p.add_argument("--registration-binning", type=int, default=d["registration_binning"],
                   help="y/x binning for registration")
    p.add_argument("--quality-threshold", type=float, default=d["quality_threshold"],
                   help="drop pairwise registrations with a lower correlation")
    p.add_argument("--min-overlap-um", type=float, default=d["min_overlap_um"],
                   help="register only wells whose neighbouring fields overlap at least this much")
    p.add_argument("--outlier-um", type=float, default=d["outlier_um"],
                   help="a field whose shift differs this much from its neighbours' median shift is an outlier")
    p.add_argument("--fix-outliers", action=argparse.BooleanOptionalAction, default=d["fix_outliers"],
                   help="move outliers along with their neighbours (--no-fix-outliers only flags them)")
    p.add_argument("--chunk", type=int, default=d["chunk"], help="output chunk size in y/x")
    p.add_argument("--ngff-version", choices=("0.4", "0.5"), default=d["ngff_version"],
                   help="NGFF version of the .ome.zarr files (the icechunk repo is always 0.5)")
    p.add_argument("--workers", type=int, default=d["workers"], help="dask threads")
    p.add_argument("--keep-tile-cache", action="store_true", default=d["keep_tile_cache"],
                   help="keep OUTPUT_DIR/_tiles/<name>.zarr after a well is stitched")
    p.add_argument("--overwrite", action="store_true", default=d["overwrite"],
                   help="restitch wells that already have an output (otherwise they are skipped)")
    p.add_argument("--icechunk", type=Path, default=d["icechunk"],
                   help="stitched icechunk repo; OUTPUT_DIR/<plate>[_z<method>]_stitched.icechunk if not given")
    args = p.parse_args(argv)
    args.projection = None if args.projection == "none" else args.projection
    return Options(**vars(args))


def main(argv: list[str] | None = None) -> None:
    opts = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    logging.getLogger("multiview_stitcher").setLevel(logging.WARNING)
    # multiview-stitcher's pyramid writer rechunks to the zarr chunks; dask warns about it on every level
    warnings.filterwarnings("ignore", category=da.PerformanceWarning)
    zarr.config.set({"async.concurrency": 16})
    dask.config.set(scheduler="threads", num_workers=opts.workers)
    log.info("multiview-stitcher %s, phenix-icechunk %s", version("multiview-stitcher"), version("phenix-icechunk"))
    results = PlateStitcher(opts).run()
    if len(results):
        print(results.set_index("well").round(2).to_string())


if __name__ == "__main__":
    main()
