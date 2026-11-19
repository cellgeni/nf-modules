#!/usr/bin/env python
"""Collect stitched wells into one OME-Zarr 0.5 HCS plate in a new icechunk repo.

The images are per-well stitching results of phenix-icechunk plates, either

* OME-TIFFs (e.g. ASHLAR, one image registered across several rounds), or
* OME-Zarrs (e.g. multiview-stitcher, one image per round, 2D or 3D),

listed in a JSON manifest::

    {
      "sources": [{"repo": "r1.icechunk", "round_index": "1"}, ...],
      "images":  [{"path": "A1.ome.tif", "well": "A1", "rounds": ["1", "2"]}, ...]
    }

``rounds`` names the source repos (by ``round_index``) whose channels the image holds, in channel
order. Every distinct ``rounds`` becomes one plate acquisition, so a well can hold one image per
round (``/<row>/<col>/<i>``) as well as a registered image of all rounds. Pixels are written as
real chunks, normalised to t/c/z/y/x, with a resolution pyramid (the input's, or built here).

The plate keeps the source layout and metadata: rows/columns and the plate/well ``phenix``
attributes come from the first source repo. The raw-data provenance of every round is copied from
its source repo (metadata only, no pixels), so the new repo stays traceable on its own::

    /sources/r<round>                    plate ``phenix`` attrs + source snapshot and commit metadata
    /sources/r<round>/phenix             Index.xml contents: images table (one row per plane with its
                                         TIFF file and tags) and the export's non-TIFF files
    /sources/r<round>/OME                generated OME-XML of the source plate
    /sources/r<round>/<row>/<col>/<f>    ``phenix`` attrs of the source well and its fields

Each image records in ``attrs["phenix_stitch"]`` its input file and, per round, the source
snapshot, Index.xml, measurement directory and the ``/sources/r<round>`` group, so every well can
be traced back to the raw TIFFs it came from. One commit is made per round's provenance and per well.
"""

from __future__ import annotations

import argparse
import json
import logging
import time
from importlib.metadata import version
from pathlib import Path
from typing import Any

import dask
import dask.array as da
import icechunk as ic
import tifffile
import xarray as xr
import zarr
from icechunk.xarray import to_icechunk

log = logging.getLogger("stitched_plate")

AXES = ("t", "c", "z", "y", "x")
UNITS = {"t": ("time", None), "c": ("channel", None), "z": ("space", "micrometer"),
         "y": ("space", "micrometer"), "x": ("space", "micrometer")}


def well_path(well: str) -> str:
    # "A4" -> "A/4"
    if "/" in well:
        return well
    row = well.rstrip("0123456789")
    return f"{row}/{int(well[len(row):])}"


def open_source(repo_path: str) -> dict[str, Any]:
    # Metadata only, so the virtual chunks (and the source export) are never read.
    # Commit metadata is merged oldest first: later commits (e.g. z projections) lack index_xml.
    repo = ic.Repository.open(ic.local_filesystem_storage(repo_path))
    ancestry = list(repo.ancestry(branch="main"))
    commit = {}
    for snapshot in reversed(ancestry):
        commit |= {k: v for k, v in (snapshot.metadata or {}).items() if not k.startswith("__")}
    plate = zarr.open_group(repo.readonly_session(branch="main").store, mode="r", zarr_format=3)
    files = plate["phenix/files"].attrs.get("root") if "phenix/files" in plate else None
    return {"plate": plate, "snapshot": ancestry[0].id,
            "commit": commit | ({"measurement_dir": files} if files else {})}


def copy_tree(src: zarr.Group, dst: zarr.Group) -> None:
    # Groups, attributes and arrays with their data, e.g. the /phenix metadata tables
    dst.attrs.put(src.attrs.asdict())
    for name, member in src.members():
        if isinstance(member, zarr.Group):
            copy_tree(member, dst.require_group(name))
            continue
        arr = dst.create_array(
            name, shape=member.shape, dtype=member.metadata.data_type, chunks=member.chunks,
            fill_value=member.fill_value, attributes=member.attrs.asdict(),
            dimension_names=member.metadata.dimension_names,
        )
        if member.size:
            arr[...] = member[...]


def copy_provenance(repo: ic.Repository, source: dict[str, Any], round_index: str, wells: list[str]) -> str:
    # Raw-data provenance of one round under /sources/r<round>, see the module docstring
    session = repo.writable_session("main")
    root = zarr.open_group(session.store, mode="a", zarr_format=3)
    plate = source["plate"]
    group = root.require_group(f"sources/r{round_index}")
    group.attrs.put({
        "phenix": dict(plate.attrs.get("phenix", {})),
        "source_snapshot": source["snapshot"],
        "commit": source["commit"],
    })
    for name in ("phenix", "OME"):
        if name in plate:
            copy_tree(plate[name], group.require_group(name))
    for well in wells:
        if well not in plate:
            continue
        group.require_group(well).attrs.put({"phenix": dict(plate[well].attrs.get("phenix", {}))})
        for field, image in plate[well].groups():
            group.require_group(f"{well}/{field}").attrs.put({"phenix": dict(image.attrs.get("phenix", {}))})
    return session.commit(f"Source metadata of round {round_index}", metadata={"round_index": round_index,
                                                                                "source_snapshot": source["snapshot"]})


def channels(plate: zarr.Group, well: str) -> list[dict[str, Any]]:
    if well_path(well) not in plate:
        return []
    image = next(iter(plate[well_path(well)].group_keys()))
    return list(plate[f"{well_path(well)}/{image}"].attrs.get("ome", {}).get("omero", {}).get("channels", []))


def to_tczyx(arr: da.Array, axes: str) -> xr.DataArray:
    img = xr.DataArray(arr, dims=list(axes.lower()))
    unknown = set(img.dims) - set(AXES)
    if unknown:
        raise ValueError(f"unsupported axes {sorted(unknown)} in {axes}")
    return img.expand_dims([a for a in AXES if a not in img.dims]).transpose(*AXES)


def read_ome_tiff(path: Path) -> tuple[list[xr.DataArray], list[float], list[str]]:
    # Levels of the first series, its z/y/x pixel size and channel names
    with tifffile.TiffFile(path) as tif:
        series = tif.series[0]
        levels = [
            to_tczyx(da.from_zarr(zarr.open(series.aszarr(level=i), mode="r")), series.axes)
            for i in range(len(series.levels))
        ]
        pixels = {}
        if tif.ome_metadata:
            image = tifffile.xml2dict(tif.ome_metadata)["OME"]["Image"]
            pixels = (image[0] if isinstance(image, list) else image)["Pixels"]
    names = pixels.get("Channel", [])
    names = [ch.get("Name", "") for ch in (names if isinstance(names, list) else [names])]
    size = [float(pixels.get(f"PhysicalSize{a}", 1.0)) for a in "ZYX"]
    return levels, size, names


def read_ome_zarr(path: Path) -> tuple[list[xr.DataArray], list[list[float]], dict | None, dict | None]:
    # Levels, their t/c/z/y/x scales, the omero block and the stitching record
    group = zarr.open_group(str(path), mode="r")
    attrs = group.attrs.asdict()
    ome = attrs.get("ome", attrs)
    multiscale = ome["multiscales"][0]
    axes = "".join(a["name"] for a in multiscale["axes"])
    levels, scales = [], []
    for dataset in multiscale["datasets"]:
        levels.append(to_tczyx(da.from_zarr(group[dataset["path"]]), axes))
        scale = next(t["scale"] for t in dataset["coordinateTransformations"] if t["type"] == "scale")
        by_axis = dict(zip(axes, scale))
        scales.append([float(by_axis.get(a, 1.0)) for a in AXES])
    return levels, scales, ome.get("omero"), attrs.get("phenix_stitch")


def channel_block(sources: dict[str, dict], well: str, rounds: list[str], n: int, fallback: list[str]) -> dict:
    # Channel labels/colours of the source wells, suffixed by round when an image spans several
    per_round = [(r, ch) for r in rounds for ch in channels(sources[r]["plate"], well)]
    if len(per_round) == n:
        chans = [dict(ch, label=f"{ch.get('label', i)}_r{r}" if len(rounds) > 1 else ch.get("label", str(i)))
                 for i, (r, ch) in enumerate(per_round)]
    else:
        log.warning("%s: %d channels in the image, %d in rounds %s; using the image's channel names",
                    well, n, len(per_round), rounds)
        chans = [{"label": fallback[i] if i < len(fallback) and fallback[i] else f"c{i}", "active": True}
                 for i in range(n)]
    return {"channels": chans}


def multiscale(name: str, scales: list[list[float]]) -> dict:
    return {
        "name": name,
        "axes": [{"name": a, "type": UNITS[a][0]} | ({"unit": UNITS[a][1]} if UNITS[a][1] else {}) for a in AXES],
        "datasets": [
            {"path": str(i), "coordinateTransformations": [{"type": "scale", "scale": s}]}
            for i, s in enumerate(scales)
        ],
    }


def write_levels(session: ic.Session, image: str, levels: list[xr.DataArray], chunk: int) -> int:
    # Write the input levels; a single level gets a pyramid, halving y/x until it fits one chunk
    root = zarr.open_group(session.store, mode="a", zarr_format=3)
    for n, level in enumerate(levels):
        write_level(session, image, n, level, chunk)
    n = len(levels)
    while len(levels) == 1 and max(root[f"{image}/{n - 1}"].shape[-2:]) > chunk:
        prev = root[f"{image}/{n - 1}"]
        level = xr.DataArray(da.from_zarr(prev), dims=AXES).coarsen(y=2, x=2, boundary="trim").mean()
        write_level(session, image, n, level.astype(prev.dtype), chunk)
        n += 1
    return n


def write_level(session: ic.Session, image: str, n: int, level: xr.DataArray, chunk: int) -> None:
    level = level.chunk({"t": 1, "c": 1, "z": 1, "y": chunk, "x": chunk}).rename(str(n))
    to_icechunk(level.to_dataset(), session, group=image, mode="a")


def extend(scales: list[list[float]], n: int) -> list[list[float]]:
    # Scales of the levels built by write_levels: y/x doubled per level
    last = scales[-1]
    return scales + [last[:3] + [last[3] * 2 ** j, last[4] * 2 ** j] for j in range(1, n - len(scales) + 1)]


def run(opts: argparse.Namespace) -> dict:
    manifest = json.loads(Path(opts.manifest).read_text())
    sources = {}
    for src in manifest["sources"]:
        sources[str(src["round_index"])] = open_source(src["repo"]) | {"repo": src["repo"]}
    first = next(iter(sources.values()))["plate"]
    source_plate = first.attrs["ome"]["plate"]
    phenix = dict(first.attrs.get("phenix", {}))

    # One acquisition per distinct set of rounds, ordered by first round
    order = list(sources)
    keys = sorted({tuple(str(r) for r in img["rounds"]) for img in manifest["images"]},
                  key=lambda k: [order.index(r) for r in k])
    acquisitions = [
        {"id": i, "name": f"round {k[0]}" if len(k) == 1 else f"rounds {', '.join(k)} registered",
         "description": "; ".join(sources[r]["commit"].get("index_xml", r) for r in k)}
        for i, k in enumerate(keys)
    ]

    by_well: dict[str, list[dict]] = {}
    for img in manifest["images"]:
        by_well.setdefault(well_path(img["well"]), []).append(img)
    rows, cols = [r["name"] for r in source_plate["rows"]], [c["name"] for c in source_plate["columns"]]
    wells = sorted(by_well, key=lambda w: (rows.index(w.split("/")[0]), cols.index(w.split("/")[1])))

    repo = ic.Repository.create(ic.local_filesystem_storage(opts.output))
    for round_index, source in sources.items():
        snapshot = copy_provenance(repo, source, round_index, wells)
        log.info("round %s: copied source metadata in %s", round_index, snapshot)
    plate_name = f"{source_plate.get('name', phenix.get('plate_name', 'plate'))} stitched"
    done = []
    for well in wells:
        t0 = time.perf_counter()
        session = repo.writable_session("main")
        root = zarr.open_group(session.store, mode="a", zarr_format=3)
        images = sorted(by_well[well], key=lambda img: keys.index(tuple(str(r) for r in img["rounds"])))
        entries = []
        for i, img in enumerate(images):
            rounds = [str(r) for r in img["rounds"]]
            path = Path(img["path"])
            image = f"{well}/{i}"
            root.require_group(well)
            root.require_group(image)
            if path.suffix == ".zarr" or path.is_dir():
                levels, scales, omero, record = read_ome_zarr(path)
                fallback = [ch.get("label", "") for ch in (omero or {}).get("channels", [])]
            else:
                levels, size, fallback = read_ome_tiff(path)
                factors = [levels[0].sizes["y"] / lvl.sizes["y"] for lvl in levels]
                scales = [[1.0, 1.0, size[0], size[1] * f, size[2] * f] for f in factors]
                record = None
            scales = extend(scales, write_levels(session, image, levels, opts.chunk))
            omero = channel_block(sources, img["well"], rounds, levels[0].sizes["c"], fallback)
            # attrs last: xarray resets the attrs of the group it writes to
            root[image].attrs.put({
                "ome": {"version": "0.5", "multiscales": [multiscale(f"{well.replace('/', '')}_{i}", scales)],
                        "omero": omero},
                "phenix_stitch": {
                    "input": path.name,
                    "rounds": [
                        {"round_index": r, "source_snapshot": sources[r]["snapshot"], "metadata": f"sources/r{r}"}
                        | sources[r]["commit"]
                        for r in rounds
                    ],
                } | ({"stitching": record} if record else {}),
            })
            entries.append({"path": str(i), "acquisition": keys.index(tuple(rounds))})
        root[well].attrs.put({
            "ome": {"version": "0.5", "well": {"images": entries}},
            "phenix": dict(first[well].attrs.get("phenix", {})) if well in first else {},
        })
        done.append(well)
        root.attrs.put({
            "ome": {"version": "0.5", "plate": {
                "name": plate_name,
                "rows": source_plate["rows"],
                "columns": source_plate["columns"],
                "acquisitions": acquisitions,
                "field_count": max(len(by_well[w]) for w in done),
                "wells": [
                    {"path": w, "rowIndex": rows.index(w.split("/")[0]), "columnIndex": cols.index(w.split("/")[1])}
                    for w in done
                ],
            }},
            "phenix": phenix,
            "phenix_stitch": {
                "sources": [
                    {"round_index": r, "source_snapshot": s["snapshot"], "metadata": f"sources/r{r}"} | s["commit"]
                    for r, s in sources.items()
                ],
                "wells": done,
                "phenix_icechunk": version("phenix-icechunk"),
            },
        })
        snapshot = session.commit(f"Stitched well {well} ({len(images)} image(s))", metadata={"well": well})
        log.info("%s: committed %s in %.0f s", well, snapshot, time.perf_counter() - t0)
    return {"repo": opts.output, "plate": plate_name, "wells": len(done), "acquisitions": acquisitions}


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("output", help="icechunk repo to create")
    p.add_argument("--manifest", required=True, help="JSON manifest of the source repos and stitched images")
    p.add_argument("--chunk", type=int, default=1024, help="output chunk size in y/x; built pyramids stop once a level fits one chunk")
    p.add_argument("--workers", type=int, default=1, help="dask threads")
    opts = p.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(levelname)s %(message)s")
    with dask.config.set(scheduler="threads", num_workers=opts.workers):
        print(json.dumps(run(opts), indent=2))


if __name__ == "__main__":
    main()
