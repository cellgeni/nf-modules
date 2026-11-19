#!/usr/bin/env python
"""Align an ordered series of cell sections to one reference with STalign (LDDMM) into one SpatialData object.

Follows the STalign Xenium-to-Xenium tutorial
(https://jef.works/STalign/notebooks/xenium-xenium-alignment.html) for every section: the cell
centroids of the section and of the reference are rasterized into density images, LDDMM fits an
affine map A plus a diffeomorphism, and the section's cells are mapped into the reference frame.
Each section is aligned to the reference independently.

The source data are not warped. SpatialData transformations are affine only, so each fitted map
is split into:

* an ``Affine`` from the section's coordinate system (named after the section) to the shared
  reference one (``--target-cs``), registered on the section's elements with
  ``spatialdata.transformations``. By default it is LDDMM's affine component A
  (``--affine stalign``); ``--affine similarity`` registers ``get_transformation_between_landmarks``
  (rotation, uniform scale, translation, and a reflection if needed) fitted between the section's
  cells and their LDDMM positions instead;
* the exact LDDMM positions of the cells (``<section>_cells_aligned`` and
  ``table.obsm['spatial_aligned']``), cheap because only the centroids are mapped;
* ``<section>_displacement``: the residual of the full map after the affine (``dx``, ``dy`` in
  reference units) on a grid over the section, so that any point p of the section (e.g. a
  transcript) maps to ``Affine(p) + displacement(p)`` without STalign.

Inputs (``--reference``, ``--sources``; each can be any of):

* a cell table (.csv, .csv.gz, .tsv[.gz], .parquet), e.g. Xenium's ``cells.csv.gz``;
* a Xenium output folder (its ``cells.parquet`` or ``cells.csv.gz`` is read);
* a SpatialData ``.zarr`` (e.g. from ``spatialdata_io.xenium``). Its cells are the centroids of a
  shapes element (``--sdata-shapes``; default: the one annotated by a table, else
  ``cell_circles``, else the only one), in the element's own coordinates unless ``--sdata-cs``
  is given, with the obs of the table annotating it. Its elements (images, labels, points, shapes,
  tables; ``--sdata-elements``) are copied, unwarped, into the output as ``<section>_<element>``
  and registered to the reference through their transformation to the cells' frame followed by
  the section's affine, so that e.g. the morphology image and the transcripts of every section are
  aligned in ``--target-cs`` too. The reference's elements become ``reference_<element>``.

Coordinates are assumed to be in micrometers (STalign's defaults depend on it). Landmarks are
given in each section's cell coordinates.

Sections need not overlap the reference to begin with (e.g. sections cut out of one slide keep
their slide positions): without landmarks, each fit starts from the translation that puts the
section's mean cell position on the reference's (``--init centroid``), and STalign's pre-fit
searches rotations from there. Sections with landmarks start from the landmarks' affine.

Sections are named with ``--names``, else after their file names (without extensions and a
trailing ``cells``, e.g. ``Rep1_cells.csv.gz`` -> ``Rep1``) when these are unique, else ``section_01``, ``section_02``, ... in input order.

Outputs (``PREFIX`` = ``--prefix``):

* ``PREFIX.zarr``: SpatialData with, for every section, shapes ``<section>_cells`` and
  ``<section>_cells_aligned`` and images ``<section>_density`` and ``<section>_displacement``; the
  reference's ``reference_cells`` and ``reference_density``; and one ``table`` of all the cells
  (``obs['section']``, ``obs['section_index']``: 1..N in input order and 0 for the reference;
  ``obsm['spatial']``, ``obsm['spatial_aligned']``, ``obsm['spatial_affine']``); plus the copied
  elements of SpatialData inputs (``<section>_<element>``, ``reference_<element>``).
  ``attrs['stalign']`` holds the summary.
* ``PREFIX_aligned_cells.csv.gz``: the input tables of all sections and the reference with
  ``section``, ``section_index``, ``aligned_x/y`` (LDDMM), ``affine_x/y`` and
  ``stalign_match_weight`` (weight of the matching class at the aligned position, 0-1).
* ``PREFIX_<section>_transform.npz``: STalign's ``A`` (row/column order), ``v``, ``xv0``, ``xv1``
  (to reuse the map with ``STalign.transform_points_source_to_target(xv, v, A, points_yx)``) and
  the registered affine as a 3x3 x/y matrix (``affine_xy``).
* ``PREFIX_stalign.json``: parameters, engine, timings and residuals of every section.
* ``PREFIX_<section>_qc.png``: overlays before alignment, after the affine and after LDDMM.

Example::

    stalign_align.py --reference rep2_cells.csv.gz --sources rep1_cells.csv.gz rep3_cells.csv.gz \\
        --prefix breast --threads 8
    stalign_align.py --reference rep2.zarr --sources rep1.zarr rep3.zarr --names rep1 rep3 --prefix breast
"""

from __future__ import annotations

import argparse
import json
import logging
import re
import time
from importlib.metadata import version
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import spatialdata as sd
import torch
from spatialdata.models import Image2DModel, PointsModel, ShapesModel, TableModel, get_table_keys
from spatialdata.transformations import (
    Affine,
    BaseTransformation,
    Identity,
    Scale,
    Sequence,
    Translation,
    get_transformation,
    get_transformation_between_coordinate_systems,
    get_transformation_between_landmarks,
    set_transformation,
)

import STalign

log = logging.getLogger("stalign_align")

REFERENCE = "reference"


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--reference", type=Path, required=True,
                   help="Reference (target) section: cell table (.csv, .csv.gz, .tsv[.gz], .parquet), "
                        "Xenium output folder or SpatialData .zarr")
    p.add_argument("--sources", type=Path, nargs="+", required=True,
                   help="Ordered source sections, any of the same formats (they can be mixed); each is "
                        "aligned to the reference")
    p.add_argument("--names", nargs="+", default=None, help="Section names, one per source, in order")
    p.add_argument("--landmarks", type=Path, default=None,
                   help="CSV of corresponding points with columns source_x, source_y, target_x, target_y and, "
                        "with several sources, section (a section name or its 1-based index). Landmarks of a "
                        "section initialize its affine (L, T) and add a landmark term to LDDMM, as in the "
                        "tutorial; sections without landmarks are aligned from the images alone")
    p.add_argument("--prefix", default="stalign", help="Output file prefix")

    io = p.add_argument_group("input columns")
    io.add_argument("--x-col", default="x_centroid", help="x coordinate column (micrometers)")
    io.add_argument("--y-col", default="y_centroid", help="y coordinate column (micrometers)")
    io.add_argument("--id-col", default="cell_id", help="Cell id column (used if present)")
    io.add_argument("--area-col", default="cell_area",
                    help="Cell area column for the circle radius in the SpatialData shapes (used if present)")
    io.add_argument("--radius", type=float, default=5.0, help="Circle radius when there is no area column")

    sdi = p.add_argument_group("SpatialData input")
    sdi.add_argument("--sdata-shapes", default=None,
                     help="Shapes element holding the cells (circles: centres, polygons: centroids). Default: "
                          "the shapes element annotated by a table, else cell_circles, else the only one")
    sdi.add_argument("--sdata-cs", default=None,
                     help="Coordinate system to read the cell positions in (must be in micrometers). Default: "
                          "the shapes element's own coordinates, which are micrometers for spatialdata_io.xenium "
                          "(whose 'global' system is in morphology pixels)")
    sdi.add_argument("--sdata-elements", default="all",
                     help="Elements of SpatialData inputs to copy into the output, registered to the reference: "
                          "all, none, or a comma-separated list of element names (e.g. "
                          "morphology_focus,transcripts,table)")

    m = p.add_argument_group("alignment")
    m.add_argument("--method", choices=("align_cells", "lddmm"), default="align_cells",
                   help="align_cells: STalign.align_cells, coarse-to-fine LDDMM on blurred density images "
                        "(tutorial's recommended route). lddmm: rasterize + a single STalign.LDDMM call "
                        "(tutorial's first route)")
    m.add_argument("--dx", type=float, default=None,
                   help="Raster pixel size in micrometers (default: 15 for align_cells, 30 for lddmm)")
    m.add_argument("--blur-um", type=float, default=None,
                   help="Rasterization kernel standard deviation in micrometers (default: 30 for align_cells, "
                        "STalign.rasterize's default for lddmm)")
    m.add_argument("--levels", default="4:500,2:500",
                   help="align_cells pyramid as factor:niter pairs, coarsest first")
    m.add_argument("--niter", type=int, default=None, help="lddmm iterations (STalign default if unset)")
    m.add_argument("--mixture", choices=("em", "none"), default=None,
                   help="Mixture model of the matching weights (default: em for align_cells, none for lddmm). "
                        "Use none when the reference section is complete")
    m.add_argument("--init", choices=("centroid", "identity"), default="centroid",
                   help="Initial affine of sections without landmarks: centroid, the translation that puts the "
                        "section's mean cell position on the reference's (needed when sections do not overlap "
                        "the reference, e.g. sections of one slide in slide coordinates); identity, STalign's default")
    m.add_argument("--prefit", choices=("auto", "true", "false"), default="auto",
                   help="Coarse rotation search before LDDMM (auto: unless landmarks are given)")
    for name, help_ in (("a", "Velocity field smoothness scale"), ("sigmaM", "Matching term std"),
                        ("sigmaB", "Background term std"), ("sigmaA", "Artifact term std"),
                        ("sigmaR", "Regularization std"), ("sigmaP", "Landmark term std"),
                        ("epV", "Velocity step size"), ("epL", "Linear step size"), ("epT", "Translation step size"),
                        ("sigma_disp", "Prior displacement scale in micrometers")):
        m.add_argument(f"--{name.replace('_', '-')}", dest=name, type=float, default=None,
                       help=f"{help_} (STalign default if unset)")

    r = p.add_argument_group("runtime")
    r.add_argument("--backend", choices=("cpp", "auto", "torch"), default="cpp",
                   help="LDDMM engine. cpp: compiled C++ engine (CPU, float64), fails if it cannot run; "
                        "auto: cpp when possible, else torch")
    r.add_argument("--threads", type=int, default=None, help="CPU threads (default: STALIGN_NUM_THREADS or all)")

    o = p.add_argument_group("output")
    o.add_argument("--affine", choices=("stalign", "similarity"), default="stalign",
                   help="Transformation registered from each section to the reference coordinate system: "
                        "STalign's affine component A, or spatialdata's get_transformation_between_landmarks "
                        "(similarity) between the section's cells and their LDDMM positions")
    o.add_argument("--target-cs", default="aligned", help="Coordinate system of the reference (aligned) frame")
    o.add_argument("--field-dx", type=float, default=50.0,
                   help="Pixel size of the <section>_displacement fields in micrometers")
    o.add_argument("--no-qc", action="store_true", help="Skip the QC plots")
    return p.parse_args(argv)


# ----------------------------------------------------------------------------- input


def read_table(path: Path) -> pd.DataFrame:
    name = path.name.lower()
    if name.endswith((".parquet", ".pq")):
        return pd.read_parquet(path)
    sep = "\t" if ".tsv" in name or ".txt" in name else ","
    return pd.read_csv(path, sep=sep)


XENIUM_CELLS = ("cells.parquet", "cells.csv.gz")
# columns of the output tables; input columns of these names are dropped
RESERVED = ("section", "section_index", "region", "instance_id", "aligned_x", "aligned_y", "affine_x", "affine_y",
            "stalign_match_weight")


def is_spatialdata(path: Path) -> bool:
    return path.is_dir() and any((path / f).exists() for f in ("zarr.json", ".zgroup", ".zattrs"))


def read_cells(path: Path, args) -> dict:
    """Cells of a section: the table df (all rows) with coordinates x, y and the finite ones ok."""
    if is_spatialdata(path):
        cells = read_sdata_cells(path, args)
    elif path.is_dir():
        found = [path / f for f in XENIUM_CELLS if (path / f).exists()]
        if not found:
            raise SystemExit(f"{path}: a folder that is neither a SpatialData .zarr nor a Xenium output folder "
                             f"(with {' or '.join(XENIUM_CELLS)})")
        log.info("%s: Xenium output folder, reading %s", path, found[0].name)
        cells = dict(kind="xenium", df=read_table(found[0]))
    else:
        cells = dict(kind="table", df=read_table(path))
    df = cells["df"]
    for c in df.columns:  # Xenium cells.parquet of older versions stores ids as bytes
        if df[c].dtype == object and len(df) and isinstance(df[c].iloc[0], bytes):
            df[c] = df[c].str.decode("utf-8")
    df.drop(columns=[c for c in RESERVED if c in df.columns], inplace=True)
    missing = [c for c in (args.x_col, args.y_col) if c not in df.columns]
    if missing:
        raise SystemExit(f"{path}: missing column(s) {missing}; columns are {list(df.columns)}. "
                         "Set --x-col/--y-col.")
    x = df[args.x_col].to_numpy(dtype=np.float64)
    y = df[args.y_col].to_numpy(dtype=np.float64)
    ok = np.isfinite(x) & np.isfinite(y)
    if not ok.all():
        log.warning("%s: %d of %d cells have non-finite coordinates and are not aligned", path, (~ok).sum(), len(df))
    if ok.sum() < 3:
        raise SystemExit(f"{path}: fewer than 3 cells with finite coordinates")
    log.info("%s: %d cells, x %.1f-%.1f, y %.1f-%.1f", path, ok.sum(),
             x[ok].min(), x[ok].max(), y[ok].min(), y[ok].max())
    cells.update(path=path, x=x, y=y, ok=ok)
    return cells


def table_regions(table) -> list[str]:
    if "spatialdata_attrs" not in table.uns:
        return []
    region = get_table_keys(table)[0]
    return [region] if isinstance(region, str) else list(region)


def pick_cells_element(sdata, path: Path, wanted: str | None) -> str:
    shapes = list(sdata.shapes)
    if wanted is not None:
        if wanted not in sdata.shapes:
            raise SystemExit(f"{path}: no shapes element {wanted!r}; shapes are {shapes}")
        return wanted
    annotated = sorted({r for t in sdata.tables.values() for r in table_regions(t) if r in sdata.shapes})
    for candidates in (annotated, [s for s in shapes if s == "cell_circles"], shapes):
        if len(candidates) == 1:
            return candidates[0]
    raise SystemExit(f"{path}: cannot tell which shapes element holds the cells (shapes: {shapes}, annotated by a "
                     f"table: {annotated}); set --sdata-shapes")


def read_sdata_cells(path: Path, args) -> dict:
    """Cells of a SpatialData object: centroids of a shapes element and the obs of its table."""
    sdata = sd.read_zarr(path)
    el = pick_cells_element(sdata, path, args.sdata_shapes)
    gdf = sdata.shapes[el]
    to_frame = Identity()
    if args.sdata_cs is not None:
        systems = get_transformation(gdf, get_all=True)
        if args.sdata_cs not in systems:
            raise SystemExit(f"{path}: {el} is not in coordinate system {args.sdata_cs!r} ({list(systems)})")
        to_frame = systems[args.sdata_cs]
    centroids = gdf.geometry.centroid
    x, y = transform_xy(to_frame, centroids.x.to_numpy(float), centroids.y.to_numpy(float))
    if "radius" in gdf.columns:
        radius = gdf["radius"].to_numpy(float)
    else:
        radius = np.sqrt(gdf.geometry.area.to_numpy(float) / np.pi)
    lin = xy_matrix(to_frame)[:2, :2]
    radius = radius * np.sqrt(abs(np.linalg.det(lin)))

    df, table_name = pd.DataFrame(index=pd.RangeIndex(len(gdf))), None
    for tname, table in sdata.tables.items():
        if el not in table_regions(table):
            continue
        _, rkey, ikey = get_table_keys(table)
        obs = table.obs[table.obs[rkey].astype(str) == el]
        pos = pd.Index(gdf.index.astype(str)).get_indexer(obs[ikey].astype(str))
        obs = obs.loc[pos >= 0].drop(columns=[rkey]).set_axis(pos[pos >= 0])
        df = obs.reindex(pd.RangeIndex(len(gdf)))
        table_name = tname
        log.info("%s: obs of table %r matched for %d of %d cells", path, tname, (pos >= 0).sum(), len(gdf))
        break
    df = df.reset_index(drop=True)
    if args.id_col not in df.columns:
        df.insert(0, args.id_col, gdf.index.astype(str).to_numpy())
    df[args.x_col], df[args.y_col] = x, y
    log.info("%s: SpatialData, cells from shapes %r in %s, %s", path, el,
             f"coordinate system {args.sdata_cs!r}" if args.sdata_cs else "its own coordinates",
             f"obs from table {table_name!r}" if table_name else "no table")
    return dict(kind="spatialdata", df=df, radius=radius, sdata=sdata, element=el, table=table_name)


def stem(path: Path) -> str:
    """File name without extensions or a trailing 'cells' (Rep1_cells.csv.gz -> Rep1)."""
    name = path.name
    for ext in (".zarr", ".gz", ".parquet", ".pq", ".csv", ".tsv", ".txt"):
        name = name.removesuffix(ext)
    return re.sub(r"[_.-]?cells$", "", name) or name


def section_names(args) -> list[str]:
    n = len(args.sources)
    if args.names is not None:
        if len(args.names) != n:
            raise SystemExit(f"--names has {len(args.names)} names for {n} sources")
        names = list(args.names)
    else:
        names = [stem(p) for p in args.sources]
        if len(set(names)) != n:
            names = [f"section_{i:0{max(2, len(str(n)))}d}" for i in range(1, n + 1)]
    # element and coordinate system names: letters, digits, '_', '-', '.'
    names = [re.sub(r"[^\w.-]", "_", s) for s in names]
    if len(set(names)) != n:
        raise SystemExit(f"section names are not unique: {names}")
    if REFERENCE in names or args.target_cs in names:
        raise SystemExit(f"section names must differ from {REFERENCE!r} and --target-cs {args.target_cs!r}")
    return names


def read_landmarks(path: Path, names: list[str]) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Landmarks per section name, as (pointsI, pointsJ) in STalign's (y, x) order."""
    lm = pd.read_csv(path)
    cols = ["source_x", "source_y", "target_x", "target_y"]
    if any(c not in lm.columns for c in cols):
        raise SystemExit(f"{path}: landmarks need columns {cols}, got {list(lm.columns)}")
    if "section" not in lm.columns:
        if len(names) > 1:
            raise SystemExit(f"{path}: with several sources the landmarks need a 'section' column")
        lm["section"] = names[0]
    by_index = {str(i): s for i, s in enumerate(names, 1)}
    lm["section"] = [s if s in names else by_index.get(s, s) for s in lm["section"].astype(str)]
    unknown = sorted(set(lm["section"]) - set(names))
    if unknown:
        raise SystemExit(f"{path}: unknown section(s) {unknown}; sections are {names} (or 1..{len(names)})")
    out = {}
    for s, g in lm.groupby("section", sort=False):
        g = g[cols].astype(float).dropna()
        if len(g):
            out[s] = (g[["source_y", "source_x"]].to_numpy(), g[["target_y", "target_x"]].to_numpy())
            log.info("%s: %d landmarks (%s)", s, len(g), "affine" if len(g) >= 3 else "translation only")
    return out


# ----------------------------------------------------------------------------- alignment


def numpy(x) -> np.ndarray:
    if isinstance(x, torch.Tensor):
        return x.detach().cpu().numpy()
    return np.asarray(x)


def lddmm_kwargs(args, landmarks, init=None) -> dict:
    kw = dict(device="cpu", dtype=torch.float64, backend=args.backend, draw=False)
    if args.threads:
        kw["num_threads"] = args.threads
    for name in ("a", "sigmaM", "sigmaB", "sigmaA", "sigmaR", "sigmaP", "epV", "epL", "epT", "sigma_disp"):
        if getattr(args, name) is not None:
            kw[name] = getattr(args, name)
    if args.mixture is not None:
        kw["mixture"] = None if args.mixture == "none" else "em"
    if args.prefit != "auto":
        kw["prefit"] = args.prefit == "true"
    if landmarks is not None:
        pointsI, pointsJ = landmarks
        L, T = STalign.L_T_from_points(pointsI, pointsJ)
        kw.update(L=torch.tensor(L), T=torch.tensor(T), pointsI=torch.tensor(pointsI), pointsJ=torch.tensor(pointsJ))
    elif init is not None:
        L, T = STalign.L_T_from_points(*init)
        kw.update(L=torch.tensor(L), T=torch.tensor(T))
    return kw


def centroid_init(xI, yI, xJ, yJ) -> tuple[np.ndarray, np.ndarray]:
    """Point pairs (pointsI, pointsJ; STalign's (y, x) order) of the translation from the source's median
    cell position to the target's. Only L and T are taken from them: they are not landmarks, so the
    pre-fit still searches rotations, about the target's centroid, from this start."""
    cI = np.array([np.mean(yI), np.mean(xI)])
    cJ = np.array([np.mean(yJ), np.mean(xJ)])
    unit = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    return cI + unit, cJ + unit


def parse_levels(text: str) -> list[tuple[int, int]]:
    try:
        return [(int(f), int(n)) for f, n in (lvl.split(":") for lvl in text.split(","))]
    except ValueError as e:
        raise SystemExit(f"--levels must look like 4:500,2:500, got {text!r}") from e


def run_alignment(xI, yI, xJ, yJ, args, landmarks) -> dict:
    init = centroid_init(xI, yI, xJ, yJ) if landmarks is None and args.init == "centroid" else None
    kw = lddmm_kwargs(args, landmarks, init)
    if args.method == "align_cells":
        dx = args.dx or 15.0
        blur_um = args.blur_um if args.blur_um is not None else 30.0
        out = STalign.align_cells(xI, yI, xJ, yJ, dx=dx, blur_um=blur_um, levels=parse_levels(args.levels), **kw)
        out["xJ_weights"] = out["levels"][-1]["xJ"]
    else:
        dx = args.dx or 30.0
        rkw = {} if args.blur_um is None else {"blur_um": args.blur_um}
        XI, YI, I = STalign.rasterize(xI, yI, dx=dx, draw=0, **rkw)
        XJ, YJ, J = STalign.rasterize(xJ, yJ, dx=dx, draw=0, **rkw)
        if args.niter is not None:
            kw["niter"] = args.niter
        out = STalign.LDDMM([YI, XI], I, [YJ, XJ], J, **kw)
        out["rasters"] = dict(XI=XI, YI=YI, I=I, XJ=XJ, YJ=YJ, J=J)
        out["xJ_weights"] = [YJ, XJ]
    out["dx"] = dx
    return out


def source_to_target(out, x, y) -> tuple[np.ndarray, np.ndarray]:
    """Full LDDMM map of source points (x, y) into the target frame; NaN stays NaN."""
    pts = np.stack([y, x], 1)
    t = numpy(STalign.transform_points_source_to_target(out["xv"], out["v"], out["A"], pts))
    return t[:, 1], t[:, 0]


def match_weight(out, ax, ay) -> np.ndarray:
    """Matching-class weight WM at the aligned positions (as in the tutorial); 0 outside the target."""
    WM = torch.as_tensor(out["WM"], dtype=torch.float64)
    WM = WM.reshape(-1, *WM.shape[-2:])
    xJ = [torch.as_tensor(numpy(v), dtype=torch.float64) for v in out["xJ_weights"]]
    phii = torch.as_tensor(np.stack([ay, ax]), dtype=torch.float64)[:, None, :]
    w = numpy(STalign.interp(xJ, WM, phii))[0, 0]
    return np.where(np.isfinite(ax), w, np.nan)


# ----------------------------------------------------------------------------- spatialdata transformations


def points(x, y, transformations=None):
    return PointsModel.parse(np.c_[x, y], transformations=transformations)


def transform_xy(t: BaseTransformation, x, y) -> tuple[np.ndarray, np.ndarray]:
    """Coordinates (x, y) mapped by a spatialdata transformation, with spatialdata.transform."""
    moved = sd.transform(points(x, y, {"_to": t}), to_coordinate_system="_to")[["x", "y"]].compute()
    return moved["x"].to_numpy(), moved["y"].to_numpy()


def xy_matrix(t: BaseTransformation) -> np.ndarray:
    return t.to_affine_matrix(input_axes=("x", "y"), output_axes=("x", "y"))


def registered_transformations(out, sx, sy, ax, ay) -> dict[str, BaseTransformation]:
    """The source -> target affine candidates as spatialdata transformations.

    stalign: LDDMM's affine A, which STalign defines on (row, column) = (y, x) coordinates.
    similarity: get_transformation_between_landmarks with every cell as a landmark (moving: its
    source position, reference: its LDDMM position).
    """
    ok = np.isfinite(ax) & np.isfinite(ay)
    A_yx = Affine(numpy(out["A"]).astype(np.float64), input_axes=("y", "x"), output_axes=("y", "x"))
    return {
        # same map, stored with x/y axes like the elements it is registered on
        "stalign": Affine(xy_matrix(A_yx), input_axes=("x", "y"), output_axes=("x", "y")),
        "similarity": get_transformation_between_landmarks(references_coords=points(ax[ok], ay[ok]),
                                                           moving_coords=points(sx[ok], sy[ok])),
    }


def grid_transform(x0, y0, step) -> Sequence:
    """Pixel (row, col) -> coordinates, with pixel centres at x0 + col*step, y0 + row*step."""
    return Sequence([
        Scale([step, step], axes=("y", "x")),
        Translation([y0 - step / 2, x0 - step / 2], axes=("y", "x")),
    ])


def rms(dx, dy) -> float:
    return float(np.sqrt(np.nanmean(dx**2 + dy**2)))


def displacement_field(out, t, x, y, step):
    """Residual of the full map after the registered transformation, on a grid over the source cells.

    Returns the (2, H, W) [dx, dy] field and the centre (x0, y0) of pixel (0, 0).
    """
    pad = 0.05 * max(np.ptp(x), np.ptp(y))
    gx = np.arange(np.min(x) - pad, np.max(x) + pad + step, step)
    gy = np.arange(np.min(y) - pad, np.max(y) + pad + step, step)
    GY, GX = np.meshgrid(gy, gx, indexing="ij")
    lx, ly = source_to_target(out, GX.ravel(), GY.ravel())
    fx, fy = transform_xy(t, GX.ravel(), GY.ravel())
    field = np.stack([(lx - fx).reshape(GX.shape), (ly - fy).reshape(GX.shape)]).astype(np.float32)
    return field, gx[0], gy[0]


def align_section(name, src, ref, args, landmarks) -> dict:
    """Align one section to the reference; adds the results to src and returns its summary."""
    sx, sy, sok = src["x"], src["y"], src["ok"]
    t0 = time.perf_counter()
    out = run_alignment(sx[sok], sy[sok], ref["x"][ref["ok"]], ref["y"][ref["ok"]], args, landmarks)
    t_fit = time.perf_counter() - t0
    log.info("%s: %s finished in %.1f s", name, args.method, t_fit)

    ax, ay, w = np.full_like(sx, np.nan), np.full_like(sy, np.nan), np.full_like(sx, np.nan)
    ax[sok], ay[sok] = source_to_target(out, sx[sok], sy[sok])
    w[sok] = match_weight(out, ax[sok], ay[sok])

    candidates = registered_transformations(out, sx[sok], sy[sok], ax[sok], ay[sok])
    resid = {}
    for k, cand in candidates.items():
        cx, cy = transform_xy(cand, sx[sok], sy[sok])
        resid[k] = rms(ax[sok] - cx, ay[sok] - cy)
    log.info("%s: RMS distance between the LDDMM map and the affine part: %s", name,
             ", ".join(f"{k} {v:.2f}" for k, v in resid.items()))
    t = candidates[args.affine]
    fx, fy = np.full_like(sx, np.nan), np.full_like(sy, np.nan)
    fx[sok], fy[sok] = transform_xy(t, sx[sok], sy[sok])

    src.update(out=out, t=t, ax=ax, ay=ay, fx=fx, fy=fy, w=w,
               field=displacement_field(out, t, sx[sok], sy[sok], args.field_dx))
    summary = {
        "source": str(src["path"]),
        "input": src["kind"],
        **({"cells_element": src["element"], "table": src["table"]} if src["kind"] == "spatialdata" else {}),
        "n_cells": int(sok.sum()),
        "init": "landmarks" if landmarks is not None else args.init,
        "seconds_fit": t_fit,
        "dx": out["dx"],
        "affine_xy": {k: xy_matrix(v).tolist() for k, v in candidates.items()},
        "rms_lddmm_vs_affine": resid,
        "match_weight_median": float(np.nanmedian(w)),
        "settings": jsonable(out.get("settings", {})),
        "prefit": jsonable({k: v for k, v in out.get("prefit", {}).items() if k in ("ran", "reason", "best", "seconds")}),
    }
    if landmarks is not None:
        pointsI, pointsJ = landmarks
        tl = numpy(STalign.transform_points_source_to_target(out["xv"], out["v"], out["A"], pointsI))
        err = np.linalg.norm(tl - pointsJ, axis=1)
        summary["landmark_error"] = {"n": len(err), "mean": float(err.mean()), "max": float(err.max())}
        log.info("%s: landmark error after alignment: mean %.2f, max %.2f", name, err.mean(), err.max())
    return summary


# ----------------------------------------------------------------------------- outputs


def cell_radius(cells, args) -> np.ndarray:
    """Circle radius of the cells with finite coordinates: the input shapes', else from the area column."""
    df, ok = cells["df"], cells["ok"]
    if "radius" in cells:
        r = cells["radius"][ok]
    elif args.area_col in df.columns:
        r = np.sqrt(np.clip(df.loc[ok, args.area_col].to_numpy(dtype=np.float64), 0, None) / np.pi)
    else:
        return np.full(int(ok.sum()), args.radius)
    return np.where(np.isfinite(r) & (r > 0), r, args.radius)


def section_frame(name, index, cells, args) -> pd.DataFrame:
    """The input table with the alignment columns (all rows, NaN where not aligned)."""
    df = cells["df"].copy()
    df.insert(0, "section_index", index)
    df.insert(0, "section", name)
    df["aligned_x"], df["aligned_y"] = cells["ax"], cells["ay"]
    df["affine_x"], df["affine_y"] = cells["fx"], cells["fy"]
    df["stalign_match_weight"] = cells["w"]
    return df


def cells_table(frames: dict[str, pd.DataFrame], cells: dict[str, dict], order: list[str], args) -> ad.AnnData:
    """One AnnData of the cells with finite coordinates of every section, annotating <section>_cells."""
    parts = []
    for name in order:
        df = frames[name].loc[cells[name]["ok"]].copy()
        df["region"] = f"{name}_cells"
        df["instance_id"] = np.arange(len(df))
        parts.append(df)
    obs = pd.concat(parts, ignore_index=True)
    if args.id_col in obs.columns and (obs["section"].astype(str) + ":" + obs[args.id_col].astype(str)).is_unique:
        obs.index = (obs["section"].astype(str) + ":" + obs[args.id_col].astype(str)).to_numpy()
    else:
        obs.index = obs.index.astype(str)
    obs.index.name = None
    # plain object strings: anndata writes them reliably whatever pandas' default string dtype
    for c in obs.columns:
        if pd.api.types.is_string_dtype(obs[c]) or obs[c].dtype == object:
            obs[c] = obs[c].astype(object)
    obs["section"] = pd.Categorical(obs["section"], categories=order, ordered=True)
    obs["region"] = pd.Categorical(obs["region"], categories=[f"{s}_cells" for s in order])
    adata = ad.AnnData(obs=obs.drop(columns=["aligned_x", "aligned_y", "affine_x", "affine_y"]))
    adata.obsm["spatial"] = obs[[args.x_col, args.y_col]].to_numpy(float)
    adata.obsm["spatial_aligned"] = obs[["aligned_x", "aligned_y"]].to_numpy(float)
    adata.obsm["spatial_affine"] = obs[["affine_x", "affine_y"]].to_numpy(float)
    return adata


def density_image(rasters, keys, dx, transformations_of):
    """STalign's density raster (keys: x axis, y axis, image) with transformations_of(pixel -> coordinates)."""
    X, Y, I = (numpy(rasters[k]) for k in keys)
    grid = grid_transform(float(X[0]), float(Y[0]), dx)
    return Image2DModel.parse(I.astype(np.float32), dims=("c", "y", "x"), transformations=transformations_of(grid))


ELEMENT_TYPES = ("images", "labels", "points", "shapes", "tables")


def selected_elements(args) -> set[str] | None:
    """Names of the SpatialData input elements to copy (None: all)."""
    text = args.sdata_elements.strip()
    if text == "all":
        return None
    return set() if text == "none" else {e.strip() for e in text.split(",") if e.strip()}


def unique_name(base: str, taken: set[str]) -> str:
    name, i = base, 2
    while name in taken:
        name, i = f"{base}_{i}", i + 1
    taken.add(name)
    return name


def carry_elements(prefix, cells, args, transformations_of, taken) -> tuple[dict, dict]:
    """Elements of a SpatialData input, renamed <prefix>_<name> and registered with transformations_of(T),
    T being the element's transformation to the frame the cells were read in.

    Returns {element type: {new name: element}} and {old name: new name}.
    """
    sdata, path, wanted = cells["sdata"], cells["path"], selected_elements(args)
    frame = args.sdata_cs if args.sdata_cs is not None else sdata.shapes[cells["element"]]
    # all paths to the frame first: registering an element changes the transformation graph
    to_frame = {}
    for etype, name, el in sdata.gen_spatial_elements():
        if wanted is not None and name not in wanted:
            continue
        if args.sdata_cs is None and name == cells["element"]:
            to_frame[name] = (etype, el, Identity())
            continue
        try:
            to_frame[name] = (etype, el, get_transformation_between_coordinate_systems(sdata, el, frame))
        except Exception as e:  # noqa: BLE001 - spatialdata raises several types for a missing path
            log.warning("%s: %s %r has no transformation to the cells' frame (%s); not copied", path, etype, name, e)
    out, renamed = {k: {} for k in ELEMENT_TYPES}, {}
    for name, (etype, el, t) in to_frame.items():
        set_transformation(el, transformations_of(t), set_all=True)
        renamed[name] = unique_name(f"{prefix}_{name}", taken)
        out[etype][renamed[name]] = el
    for name, table in sdata.tables.items():
        if wanted is not None and name not in wanted:
            continue
        regions = table_regions(table)
        if regions:
            keep = [r for r in regions if r in renamed]
            if not keep:
                log.warning("%s: table %r annotates no copied element (%s); not copied", path, name, regions)
                continue
            _, rkey, ikey = get_table_keys(table)
            table = table[table.obs[rkey].astype(str).isin(keep)].copy()
            table.obs[rkey] = table.obs[rkey].astype(str).map(renamed)
            del table.uns["spatialdata_attrs"]
            table = TableModel.parse(table, region=[renamed[r] for r in keep], region_key=rkey, instance_key=ikey)
        renamed[name] = unique_name(f"{prefix}_{name}", taken)
        out["tables"][renamed[name]] = table
    if renamed:
        log.info("%s: copied %s", path, ", ".join(f"{k} -> {v}" for k, v in renamed.items()))
    return out, renamed


def write_spatialdata(path, args, names, sections, ref, frames, summary):
    tcs = args.target_cs
    shapes, images = {}, {}
    for name in names:
        s = sections[name]
        sok, t = s["ok"], s["t"]
        r = cell_radius(s, args)
        aok = np.isfinite(s["ax"][sok]) & np.isfinite(s["ay"][sok])
        shapes[f"{name}_cells"] = ShapesModel.parse(
            np.c_[s["x"][sok], s["y"][sok]], geometry=0, radius=r, index=np.arange(sok.sum()),
            transformations={name: Identity(), tcs: t})
        shapes[f"{name}_cells_aligned"] = ShapesModel.parse(
            np.c_[s["ax"][sok][aok], s["ay"][sok][aok]], geometry=0, radius=r[aok], index=np.flatnonzero(aok),
            transformations={tcs: Identity()})
        out = s["out"]
        images[f"{name}_density"] = density_image(
            out["rasters"], ("XI", "YI", "I"), out["dx"], lambda g, n=name, t=t: {n: g, tcs: Sequence([g, t])})
        f, fx0, fy0 = s["field"]
        images[f"{name}_displacement"] = Image2DModel.parse(
            f, dims=("c", "y", "x"), c_coords=["dx", "dy"],
            transformations={name: grid_transform(fx0, fy0, args.field_dx)})

    rok = ref["ok"]
    shapes[f"{REFERENCE}_cells"] = ShapesModel.parse(
        np.c_[ref["x"][rok], ref["y"][rok]], geometry=0, radius=cell_radius(ref, args),
        index=np.arange(rok.sum()), transformations={tcs: Identity()})
    # every section rasterizes the reference on the same grid; take the first one's
    first = sections[names[0]]["out"]
    images[f"{REFERENCE}_density"] = density_image(
        first["rasters"], ("XJ", "YJ", "J"), first["dx"], lambda g: {tcs: g})

    order = names + [REFERENCE]
    table = cells_table(frames, {**sections, REFERENCE: ref}, order, args)
    table = TableModel.parse(table, region=[f"{s}_cells" for s in order], region_key="region",
                             instance_key="instance_id")
    elements = {k: {} for k in ELEMENT_TYPES}
    elements.update(images=images, shapes=shapes, tables={"table": table})

    # elements of SpatialData inputs, registered through the cells' frame
    taken = {n for group in elements.values() for n in group}
    for name in order:
        cells = ref if name == REFERENCE else sections[name]
        if cells["kind"] != "spatialdata":
            continue
        if name == REFERENCE:
            def transformations_of(t):
                return {tcs: t}
        else:
            def transformations_of(t, n=name, aff=cells["t"]):
                return {n: t, tcs: Sequence([t, aff])}
        carried, renamed = carry_elements(name, cells, args, transformations_of, taken)
        for k in ELEMENT_TYPES:
            elements[k].update(carried[k])
        (summary if name == REFERENCE else summary["sections"][name])["copied_elements"] = renamed

    sdata = sd.SpatialData(**elements, attrs={"stalign": summary})
    sdata.write(path, overwrite=True)
    log.info("wrote %s", path)


def plot_qc(path, name, src, ref, aff_name, max_points=50000):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(0)

    def sub(n):
        return np.sort(rng.choice(n, min(n, max_points), replace=False))

    si, ti = sub(len(src["x"])), sub(len(ref["x"]))
    tx, ty = ref["x"][ti], ref["y"][ti]
    panels = [
        ("before", src["x"][si], src["y"][si], None),
        (f"affine ({aff_name})", src["fx"][si], src["fy"][si], None),
        ("LDDMM", src["ax"][si], src["ay"][si], None),
        ("LDDMM, match weight", src["ax"][si], src["ay"][si], src["w"][si]),
    ]
    # marker area shrinks with the number of points drawn: ~4 for a few thousand cells, 0.3 for 50k
    ms = float(np.clip(20000 / max(len(si), len(ti)), 0.3, 4.0))
    fig, axes = plt.subplots(1, 4, figsize=(24, 6.5))
    for ax, (title, x, y, c) in zip(axes, panels):
        ax.scatter(tx, ty, s=ms, c="0.6" if c is not None else "tab:blue", alpha=0.5, lw=0, label=REFERENCE)
        if c is None:
            ax.scatter(x, y, s=ms, c="tab:red", alpha=0.5, lw=0, label=name)
        else:
            sc = ax.scatter(x, y, s=ms, c=c, cmap="viridis", vmin=0, vmax=1, lw=0)
            fig.colorbar(sc, ax=ax, fraction=0.04, label="stalign_match_weight")
        ax.set_title(title)
        ax.set_aspect("equal")
        ax.invert_yaxis()
    axes[0].legend(markerscale=4 / np.sqrt(ms), loc="upper right")
    fig.suptitle(f"{name} -> {REFERENCE}")
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)


def jsonable(x):
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, (torch.Tensor, np.ndarray)):
        a = numpy(x)
        return a.tolist() if a.size <= 64 else f"array{list(a.shape)}"
    if isinstance(x, (np.floating, np.integer, np.bool_)):
        return x.item()
    if isinstance(x, (str, int, float, bool)) or x is None:
        return x
    return str(x)


# ----------------------------------------------------------------------------- main


def main(argv=None):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    if args.threads:
        torch.set_num_threads(args.threads)
    info = STalign.backend_info()
    log.info("STalign %s, backends %s", STalign.__version__, info)
    if args.backend == "cpp" and not info.get("cpp"):
        raise SystemExit(f"--backend cpp: the C++ engine is not available ({info.get('cpp_error')}); "
                         "use --backend auto to fall back to torch")
    engine = "cpp" if args.backend == "cpp" or (args.backend == "auto" and info.get("auto") == "cpp") else "torch"

    names = section_names(args)
    landmarks = read_landmarks(args.landmarks, names) if args.landmarks is not None else {}
    ref = read_cells(args.reference, args)
    # the reference is its own aligned frame
    ref.update(ax=ref["x"], ay=ref["y"], fx=ref["x"], fy=ref["y"], w=np.full_like(ref["x"], np.nan))

    summary = {
        "stalign_version": STalign.__version__,
        "spatialdata_version": version("spatialdata"),
        "method": args.method,
        "engine": engine,
        "backend_info": info,
        "reference": str(args.reference),
        "reference_input": ref["kind"],
        "n_reference_cells": int(ref["ok"].sum()),
        "target_cs": args.target_cs,
        "affine_used": args.affine,
        "sections": {},
        "args": jsonable(vars(args)),
    }
    sections, frames = {}, {}
    for i, (name, path) in enumerate(zip(names, args.sources), 1):
        src = read_cells(path, args)
        summary["sections"][name] = {"index": i, **align_section(name, src, ref, args, landmarks.get(name))}
        sections[name] = src
        frames[name] = section_frame(name, i, src, args)

        out = src["out"]
        xv = [numpy(v) for v in out["xv"]]
        np.savez_compressed(f"{args.prefix}_{name}_transform.npz", A=numpy(out["A"]), v=numpy(out["v"]),
                            xv0=xv[0], xv1=xv[1], affine_xy=xy_matrix(src["t"]), WM=numpy(out["WM"]),
                            xJ0=numpy(out["xJ_weights"][0]), xJ1=numpy(out["xJ_weights"][1]))
        if not args.no_qc:
            plot_qc(f"{args.prefix}_{name}_qc.png", name, src, ref, args.affine)
    frames[REFERENCE] = section_frame(REFERENCE, 0, ref, args)

    pd.concat([frames[s] for s in names + [REFERENCE]], ignore_index=True).to_csv(
        f"{args.prefix}_aligned_cells.csv.gz", index=False, compression="gzip")
    write_spatialdata(f"{args.prefix}.zarr", args, names, sections, ref, frames, summary)
    with open(f"{args.prefix}_stalign.json", "w") as fh:
        json.dump(summary, fh, indent=2)


if __name__ == "__main__":
    main()
