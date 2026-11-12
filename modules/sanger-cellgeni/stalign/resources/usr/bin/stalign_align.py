#!/usr/bin/env python
"""Align two sets of cells with STalign (LDDMM) and write the result as a SpatialData object.

Follows the STalign Xenium-to-Xenium tutorial
(https://jef.works/STalign/notebooks/xenium-xenium-alignment.html): cell centroids of a source
and a target section are rasterized into density images, LDDMM fits an affine map A plus a
diffeomorphism, and the source cells are mapped into the target frame.

The source data are not warped. SpatialData transformations are affine only, so the fitted map
is split into:

* an ``Affine`` from the source coordinate system (``--source-cs``) to the target one
  (``--target-cs``), registered on every source element with ``spatialdata.transformations``.
  By default it is LDDMM's affine component A (``--affine stalign``); ``--affine similarity``
  registers ``get_transformation_between_landmarks`` (rotation, uniform scale, translation, and a
  reflection if needed) fitted between the source cells and their LDDMM positions instead;
* the exact LDDMM positions of the source cells (``source_cells_aligned`` and
  ``table.obsm['spatial_aligned']``), cheap because only the centroids are mapped;
* ``stalign_displacement``: the residual of the full map after the affine (``dx``, ``dy`` in
  target units) on a grid over the source section, so that any source point p (e.g. a transcript)
  maps to ``Affine(p) + displacement(p)`` without STalign.

Coordinates are assumed to be in micrometers (STalign's defaults depend on it).

Outputs (``PREFIX`` = ``--prefix``):

* ``PREFIX.zarr``: SpatialData with shapes ``source_cells``, ``source_cells_aligned`` and
  ``target_cells``, images ``source_density``, ``target_density`` and ``stalign_displacement``, and
  tables ``table`` (source cells) and ``target_table``. ``attrs['stalign']`` holds the summary.
* ``PREFIX_aligned_cells.csv.gz``: the source table with ``aligned_x/y`` (LDDMM), ``affine_x/y``
  and ``stalign_match_weight`` (weight of the matching class at the aligned position, 0-1).
* ``PREFIX_transform.npz``: STalign's ``A`` (row/column order), ``v``, ``xv0``, ``xv1`` (to reuse
  the map with ``STalign.transform_points_source_to_target(xv, v, A, points_yx)``) and the
  registered affine as a 3x3 x/y matrix (``affine_xy``).
* ``PREFIX_stalign.json``: parameters, engine, timings and residuals.
* ``PREFIX_qc.png``: overlays before alignment, after the affine and after LDDMM.

Example::

    stalign_align.py rep1_cells.csv.gz rep2_cells.csv.gz --prefix rep1_to_rep2 --threads 8
"""

from __future__ import annotations

import argparse
import json
import logging
import time
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
import spatialdata as sd
import torch
from spatialdata.models import Image2DModel, PointsModel, ShapesModel, TableModel
from spatialdata.transformations import (
    Affine,
    BaseTransformation,
    Identity,
    Scale,
    Sequence,
    Translation,
    get_transformation_between_landmarks,
)

import STalign

log = logging.getLogger("stalign_align")

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("source", type=Path, help="Source cells (.csv, .csv.gz, .tsv[.gz] or .parquet); mapped onto the target")
    p.add_argument("target", type=Path, help="Target cells, same format")
    p.add_argument("--landmarks", type=Path, default=None,
                   help="CSV of corresponding points with columns source_x, source_y, target_x, target_y. "
                        "Initializes the affine (L, T) and adds a landmark term to LDDMM, as in the tutorial")
    p.add_argument("--prefix", default="stalign", help="Output file prefix")

    io = p.add_argument_group("input columns")
    io.add_argument("--x-col", default="x_centroid", help="x coordinate column (micrometers)")
    io.add_argument("--y-col", default="y_centroid", help="y coordinate column (micrometers)")
    io.add_argument("--id-col", default="cell_id", help="Cell id column (used if present)")
    io.add_argument("--area-col", default="cell_area",
                    help="Cell area column for the circle radius in the SpatialData shapes (used if present)")
    io.add_argument("--radius", type=float, default=5.0, help="Circle radius when there is no area column")

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
                        "Use none when the target section is complete")
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
                   help="Transformation registered from the source to the target coordinate system: STalign's "
                        "affine component A, or spatialdata's get_transformation_between_landmarks (similarity) "
                        "between the source cells and their LDDMM positions")
    o.add_argument("--source-cs", default="global", help="Coordinate system of the source data")
    o.add_argument("--target-cs", default="aligned", help="Coordinate system of the target (aligned) frame")
    o.add_argument("--field-dx", type=float, default=50.0,
                   help="Pixel size of the stalign_displacement field in micrometers")
    o.add_argument("--no-qc", action="store_true", help="Skip the QC plot")
    return p.parse_args(argv)


# ----------------------------------------------------------------------------- input


def read_table(path: Path) -> pd.DataFrame:
    name = path.name.lower()
    if name.endswith((".parquet", ".pq")):
        return pd.read_parquet(path)
    sep = "\t" if ".tsv" in name or ".txt" in name else ","
    return pd.read_csv(path, sep=sep)


def read_cells(path: Path, args) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    df = read_table(path)
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
    log.info("%s: %d cells, x %.1f-%.1f, y %.1f-%.1f", path.name, ok.sum(),
             x[ok].min(), x[ok].max(), y[ok].min(), y[ok].max())
    return df, x, y


def read_landmarks(path: Path):
    lm = pd.read_csv(path)
    cols = ["source_x", "source_y", "target_x", "target_y"]
    if any(c not in lm.columns for c in cols):
        raise SystemExit(f"{path}: landmarks need columns {cols}, got {list(lm.columns)}")
    lm = lm[cols].astype(float).dropna()
    if len(lm) == 0:
        raise SystemExit(f"{path}: no landmarks")
    pointsI = lm[["source_y", "source_x"]].to_numpy()
    pointsJ = lm[["target_y", "target_x"]].to_numpy()
    log.info("%d landmarks (%s)", len(lm), "affine" if len(lm) >= 3 else "translation only")
    return pointsI, pointsJ


# ----------------------------------------------------------------------------- alignment


def numpy(x) -> np.ndarray:
    if isinstance(x, torch.Tensor):
        return x.detach().cpu().numpy()
    return np.asarray(x)


def lddmm_kwargs(args, pointsI, pointsJ) -> dict:
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
    if pointsI is not None:
        L, T = STalign.L_T_from_points(pointsI, pointsJ)
        kw.update(L=torch.tensor(L), T=torch.tensor(T), pointsI=torch.tensor(pointsI), pointsJ=torch.tensor(pointsJ))
    return kw


def parse_levels(text: str) -> list[tuple[int, int]]:
    try:
        return [(int(f), int(n)) for f, n in (lvl.split(":") for lvl in text.split(","))]
    except ValueError as e:
        raise SystemExit(f"--levels must look like 4:500,2:500, got {text!r}") from e


def run_alignment(xI, yI, xJ, yJ, args, pointsI, pointsJ) -> dict:
    kw = lddmm_kwargs(args, pointsI, pointsJ)
    if args.method == "align_cells":
        dx = args.dx or 15.0
        blur_um = args.blur_um if args.blur_um is not None else 30.0
        out = STalign.align_cells(xI, yI, xJ, yJ, dx=dx, blur_um=blur_um, levels=parse_levels(args.levels), **kw)
        r = out["rasters"]
        out["xJ_weights"] = out["levels"][-1]["xJ"]
    else:
        dx = args.dx or 30.0
        rkw = {} if args.blur_um is None else {"blur_um": args.blur_um}
        XI, YI, I = STalign.rasterize(xI, yI, dx=dx, draw=0, **rkw)
        XJ, YJ, J = STalign.rasterize(xJ, yJ, dx=dx, draw=0, **rkw)
        if args.niter is not None:
            kw["niter"] = args.niter
        out = STalign.LDDMM([YI, XI], I, [YJ, XJ], J, **kw)
        r = dict(XI=XI, YI=YI, I=I, XJ=XJ, YJ=YJ, J=J)
        out["rasters"] = r
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


def registered_transformations(out, sx, sy, ax, ay) -> dict[str, BaseTransformation]:
    """The source -> target affine candidates as spatialdata transformations.

    similarity: get_transformation_between_landmarks with every cell as a landmark (moving: its
    source position, reference: its LDDMM position). stalign: LDDMM's affine A, which STalign
    defines on (row, column) = (y, x) coordinates.
    """
    ok = np.isfinite(ax) & np.isfinite(ay)
    A_yx = Affine(numpy(out["A"]).astype(np.float64), input_axes=("y", "x"), output_axes=("y", "x"))
    return {
        # same map, stored with x/y axes like the elements it is registered on
        "stalign": Affine(xy_matrix(A_yx), input_axes=("x", "y"), output_axes=("x", "y")),
        "similarity": get_transformation_between_landmarks(references_coords=points(ax[ok], ay[ok]),
                                                           moving_coords=points(sx[ok], sy[ok])),
    }


def xy_matrix(t: BaseTransformation) -> np.ndarray:
    return t.to_affine_matrix(input_axes=("x", "y"), output_axes=("x", "y"))


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


# ----------------------------------------------------------------------------- outputs


def cell_radius(df, ok, args) -> np.ndarray:
    if args.area_col in df.columns:
        area = df.loc[ok, args.area_col].to_numpy(dtype=np.float64)
        r = np.sqrt(np.clip(area, 0, None) / np.pi)
        return np.where(np.isfinite(r) & (r > 0), r, args.radius)
    return np.full(int(ok.sum()), args.radius)


def cells_adata(df, ok, args, region):
    import anndata as ad

    obs = df.loc[ok].copy()
    if args.id_col in obs.columns and obs[args.id_col].is_unique:
        obs.index = obs[args.id_col].astype(str).to_numpy()
    else:
        obs.index = obs.index.astype(str)
    obs.index.name = None
    # plain object strings: anndata writes them reliably whatever pandas' default string dtype
    for c in obs.columns:
        if pd.api.types.is_string_dtype(obs[c]) or obs[c].dtype == object:
            obs[c] = obs[c].astype(object)
    obs["region"] = pd.Categorical([region] * len(obs))
    obs["instance_id"] = np.arange(len(obs))
    adata = ad.AnnData(obs=obs)
    adata.obsm["spatial"] = np.c_[obs[args.x_col].to_numpy(float), obs[args.y_col].to_numpy(float)]
    return adata


def write_spatialdata(path, args, out, src, tgt, t, field, summary):
    scs, tcs = args.source_cs, args.target_cs
    sdf, sx, sy, sok, ax, ay, fx, fy, w = (src[k] for k in ("df", "x", "y", "ok", "ax", "ay", "fx", "fy", "w"))
    tdf, tx, ty, tok = (tgt[k] for k in ("df", "x", "y", "ok"))
    rs, rt = cell_radius(sdf, sok, args), cell_radius(tdf, tok, args)
    # rows that are finite after the LDDMM map (all of them unless points left the velocity grid)
    aok = np.isfinite(ax[sok]) & np.isfinite(ay[sok])

    shapes = {
        "source_cells": ShapesModel.parse(np.c_[sx[sok], sy[sok]], geometry=0, radius=rs,
                                          index=np.arange(sok.sum()),
                                          transformations={scs: Identity(), tcs: t}),
        "source_cells_aligned": ShapesModel.parse(np.c_[ax[sok][aok], ay[sok][aok]], geometry=0, radius=rs[aok],
                                                  index=np.flatnonzero(aok),
                                                  transformations={tcs: Identity()}),
        "target_cells": ShapesModel.parse(np.c_[tx[tok], ty[tok]], geometry=0, radius=rt,
                                          index=np.arange(tok.sum()), transformations={tcs: Identity()}),
    }

    r, dx = out["rasters"], out["dx"]
    s_grid = grid_transform(float(numpy(r["XI"])[0]), float(numpy(r["YI"])[0]), dx)
    t_grid = grid_transform(float(numpy(r["XJ"])[0]), float(numpy(r["YJ"])[0]), dx)
    f, fx0, fy0 = field
    f_grid = grid_transform(fx0, fy0, args.field_dx)
    images = {
        "source_density": Image2DModel.parse(numpy(r["I"]).astype(np.float32), dims=("c", "y", "x"),
                                             transformations={scs: s_grid, tcs: Sequence([s_grid, t])}),
        "target_density": Image2DModel.parse(numpy(r["J"]).astype(np.float32), dims=("c", "y", "x"),
                                             transformations={tcs: t_grid}),
        "stalign_displacement": Image2DModel.parse(f, dims=("c", "y", "x"), c_coords=["dx", "dy"],
                                                   transformations={scs: f_grid}),
    }

    table = cells_adata(sdf, sok, args, "source_cells")
    table.obsm["spatial_aligned"] = np.c_[ax[sok], ay[sok]]
    table.obsm["spatial_affine"] = np.c_[fx[sok], fy[sok]]
    table.obs["stalign_match_weight"] = w[sok]
    tables = {
        "table": TableModel.parse(table, region="source_cells", region_key="region", instance_key="instance_id"),
        "target_table": TableModel.parse(cells_adata(tdf, tok, args, "target_cells"), region="target_cells",
                                         region_key="region", instance_key="instance_id"),
    }

    sdata = sd.SpatialData(images=images, shapes=shapes, tables=tables, attrs={"stalign": summary})
    sdata.write(path, overwrite=True)
    log.info("wrote %s", path)


def plot_qc(path, src, tgt, aff_name, max_points=50000):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(0)

    def sub(n):
        return np.sort(rng.choice(n, min(n, max_points), replace=False))

    si, ti = sub(len(src["x"])), sub(len(tgt["x"]))
    tx, ty = tgt["x"][ti], tgt["y"][ti]
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
        ax.scatter(tx, ty, s=ms, c="0.6" if c is not None else "tab:blue", alpha=0.5, lw=0, label="target")
        if c is None:
            ax.scatter(x, y, s=ms, c="tab:red", alpha=0.5, lw=0, label="source")
        else:
            sc = ax.scatter(x, y, s=ms, c=c, cmap="viridis", vmin=0, vmax=1, lw=0)
            fig.colorbar(sc, ax=ax, fraction=0.04, label="stalign_match_weight")
        ax.set_title(title)
        ax.set_aspect("equal")
        ax.invert_yaxis()
    axes[0].legend(markerscale=4 / np.sqrt(ms), loc="upper right")
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
    if args.source_cs == args.target_cs:
        raise SystemExit("--source-cs and --target-cs must differ")
    if args.threads:
        torch.set_num_threads(args.threads)
    info = STalign.backend_info()
    log.info("STalign %s, backends %s", STalign.__version__, info)
    if args.backend == "cpp" and not info.get("cpp"):
        raise SystemExit(f"--backend cpp: the C++ engine is not available ({info.get('cpp_error')}); "
                         "use --backend auto to fall back to torch")
    engine = "cpp" if args.backend == "cpp" or (args.backend == "auto" and info.get("auto") == "cpp") else "torch"

    sdf, sx, sy = read_cells(args.source, args)
    tdf, tx, ty = read_cells(args.target, args)
    sok, tok = np.isfinite(sx) & np.isfinite(sy), np.isfinite(tx) & np.isfinite(ty)
    pointsI = pointsJ = None
    if args.landmarks is not None:
        pointsI, pointsJ = read_landmarks(args.landmarks)

    t0 = time.perf_counter()
    out = run_alignment(sx[sok], sy[sok], tx[tok], ty[tok], args, pointsI, pointsJ)
    t_fit = time.perf_counter() - t0
    log.info("%s finished in %.1f s (%s engine)", args.method, t_fit, engine)

    ax, ay = np.full_like(sx, np.nan), np.full_like(sy, np.nan)
    ax[sok], ay[sok] = source_to_target(out, sx[sok], sy[sok])
    w = np.full_like(sx, np.nan)
    w[sok] = match_weight(out, ax[sok], ay[sok])

    candidates = registered_transformations(out, sx[sok], sy[sok], ax[sok], ay[sok])
    resid = {}
    for name, cand in candidates.items():
        cx, cy = transform_xy(cand, sx[sok], sy[sok])
        resid[name] = rms(ax[sok] - cx, ay[sok] - cy)
    log.info("RMS distance between the LDDMM map and the affine part: %s",
             ", ".join(f"{k} {v:.2f}" for k, v in resid.items()))
    t = candidates[args.affine]
    fx, fy = np.full_like(sx, np.nan), np.full_like(sy, np.nan)
    fx[sok], fy[sok] = transform_xy(t, sx[sok], sy[sok])

    summary = {
        "stalign_version": STalign.__version__,
        "spatialdata_version": version("spatialdata"),
        "method": args.method,
        "engine": engine,
        "backend_info": info,
        "seconds_fit": t_fit,
        "source": str(args.source),
        "target": str(args.target),
        "n_source_cells": int(sok.sum()),
        "n_target_cells": int(tok.sum()),
        "dx": out["dx"],
        "coordinate_systems": {"source": args.source_cs, "target": args.target_cs},
        "affine_used": args.affine,
        "affine_xy": {k: xy_matrix(v).tolist() for k, v in candidates.items()},
        "rms_lddmm_vs_affine": resid,
        "match_weight_median": float(np.nanmedian(w)),
        "settings": jsonable(out.get("settings", {})),
        "prefit": jsonable({k: v for k, v in out.get("prefit", {}).items() if k in ("ran", "reason", "best", "seconds")}),
        "args": jsonable(vars(args)),
    }
    if pointsI is not None:
        tl = numpy(STalign.transform_points_source_to_target(out["xv"], out["v"], out["A"], pointsI))
        err = np.linalg.norm(tl - pointsJ, axis=1)
        summary["landmark_error"] = {"mean": float(err.mean()), "max": float(err.max())}
        log.info("landmark error after alignment: mean %.2f, max %.2f", err.mean(), err.max())

    res = sdf.copy()
    res["aligned_x"], res["aligned_y"] = ax, ay
    res["affine_x"], res["affine_y"] = fx, fy
    res["stalign_match_weight"] = w
    res.to_csv(f"{args.prefix}_aligned_cells.csv.gz", index=False, compression="gzip")

    xv = [numpy(v) for v in out["xv"]]
    np.savez_compressed(f"{args.prefix}_transform.npz", A=numpy(out["A"]), v=numpy(out["v"]), xv0=xv[0], xv1=xv[1],
                        affine_xy=xy_matrix(t), WM=numpy(out["WM"]),
                        xJ0=numpy(out["xJ_weights"][0]), xJ1=numpy(out["xJ_weights"][1]))
    with open(f"{args.prefix}_stalign.json", "w") as fh:
        json.dump(summary, fh, indent=2)

    src = dict(df=sdf, x=sx, y=sy, ok=sok, ax=ax, ay=ay, fx=fx, fy=fy, w=w)
    tgt = dict(df=tdf, x=tx, y=ty, ok=tok)
    if not args.no_qc:
        plot_qc(f"{args.prefix}_qc.png", src, tgt, args.affine)
    field = displacement_field(out, t, sx[sok], sy[sok], args.field_dx)
    write_spatialdata(f"{args.prefix}.zarr", args, out, src, tgt, t, field, summary)


if __name__ == "__main__":
    main()
