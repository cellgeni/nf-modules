#!/usr/bin/env python3
import argparse
import json
import logging
import warnings

import numpy as np
import scanpy as sc
import segtraq
import spatialdata_io


METRIC_ALIASES = {
    "baseline": "baseline",
    "clustering": "clustering_stability",
    "clusteringstability": "clustering_stability",
    "clustering_stability": "clustering_stability",
    "region": "region_similarity",
    "regionsimilarity": "region_similarity",
    "region_similarity": "region_similarity",
    "volume": "volume",
}
METRIC_ORDER = ["baseline", "region_similarity", "volume", "clustering_stability"]
LOGGER = logging.getLogger("segtraq_qc")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run one or more SegTraQ QC metric groups directly from a Xenium bundle."
    )
    parser.add_argument("--xenium_dir", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument(
        "--metrics",
        default="baseline,clustering_stability,region_similarity,volume",
        help="Comma-separated metric groups: baseline, clustering_stability, region_similarity, volume.",
    )
    parser.add_argument("--images_key", default="morphology_focus")
    parser.add_argument("--table_key", default="table")
    parser.add_argument("--tables_cell_id_key", default="cell_id")
    parser.add_argument(
        "--tables_area_key",
        default="cell_area",
        help="Obs column for cell area; use 'None' to auto-compute.",
    )
    parser.add_argument("--centroid_x_key", default=None)
    parser.add_argument("--centroid_y_key", default=None)
    parser.add_argument("--shapes_key", default="cell_boundaries")
    parser.add_argument("--shapes_cell_id_key", default="cell_id")
    parser.add_argument(
        "--nucleus_shapes_key",
        default="nucleus_boundaries",
        help="Key for nucleus shapes; use 'None' if unavailable.",
    )
    parser.add_argument("--nucleus_shapes_cell_id_key", default="cell_id")
    parser.add_argument(
        "--points_z_key",
        default="z",
        help="Z coordinate column in transcripts; use 'None' for 2D data.",
    )
    parser.add_argument("--points_gene_key", default="feature_name")
    parser.add_argument(
        "--points_background_id",
        default="UNASSIGNED",
        help="Background cell ID in transcripts; use 'None' to disable.",
    )
    parser.add_argument(
        "--min_qv",
        default="20.0",
        help="Minimum quality value for transcript filtering; use 'None' to disable.",
    )
    parser.add_argument("--leiden_resolution", type=float, default=0.5)
    parser.add_argument("--n_pcs", type=int, default=50)
    parser.add_argument("--n_neighbors", type=int, default=15)
    parser.add_argument("--z_quantile", type=float, default=0.30)
    parser.add_argument("--version", action="version", version="segtraq_qc.py 0.1.0")
    return parser.parse_args()


def _none_or_str(val):
    return None if val == "None" else val


def _parse_bg_id(val):
    if val == "None":
        return None
    try:
        return int(val)
    except ValueError:
        return val


def _parse_min_qv(val):
    return None if val == "None" else float(val)


def _parse_metrics(raw_metrics):
    if isinstance(raw_metrics, str):
        parts = raw_metrics.split(",")
    else:
        parts = raw_metrics
    requested = []
    for part in parts:
        key = str(part).strip().lower()
        if not key:
            continue
        if key in ("all", "*"):
            requested.extend(METRIC_ORDER)
            continue
        if key not in METRIC_ALIASES:
            allowed = ", ".join(METRIC_ORDER)
            raise SystemExit(f"Unknown --metrics value '{part}'. Expected one or more of: {allowed}")
        requested.append(METRIC_ALIASES[key])
    return [metric for metric in METRIC_ORDER if metric in set(requested)]


def _fix_shapes_index(sdata, shapes_key, cell_id_key):
    if not shapes_key or shapes_key not in sdata.shapes:
        return
    shape = sdata.shapes[shapes_key]
    if cell_id_key not in shape.columns and shape.index.name != cell_id_key:
        shape.index.name = cell_id_key
        sdata.shapes[shapes_key] = shape


def _build_segtraq(sdata, args):
    nsk = _none_or_str(args.nucleus_shapes_key)
    _fix_shapes_index(sdata, args.shapes_key, args.shapes_cell_id_key)
    if nsk:
        _fix_shapes_index(sdata, nsk, args.nucleus_shapes_cell_id_key)

    return segtraq.SegTraQ(
        sdata,
        images_key=_none_or_str(args.images_key),
        tables_key=args.table_key,
        tables_cell_id_key=args.tables_cell_id_key,
        tables_area_key=_none_or_str(args.tables_area_key),
        tables_centroid_x_key=args.centroid_x_key,
        tables_centroid_y_key=args.centroid_y_key,
        shapes_key=args.shapes_key,
        shapes_cell_id_key=args.shapes_cell_id_key,
        nucleus_shapes_key=nsk,
        nucleus_shapes_cell_id_key=args.nucleus_shapes_cell_id_key,
        points_z_key=_none_or_str(args.points_z_key),
        points_gene_key=args.points_gene_key,
        points_background_id=_parse_bg_id(args.points_background_id),
    )


def _obs_columns(sdata, table_key):
    return set(sdata.tables[table_key].obs.columns)


def _log_added_obs_columns(before_columns, after_columns, metric_name):
    added_columns = sorted(after_columns - before_columns)
    if added_columns:
        LOGGER.info("%s added obs column(s): %s", metric_name, ", ".join(added_columns))
    else:
        LOGGER.info("%s did not add new obs columns", metric_name)


def _write_obs(sdata, table_key, path):
    LOGGER.info("writing obs table to %s", path)
    sdata.tables[table_key].obs.to_csv(path)


def _run_baseline(st, sdata, args):
    LOGGER.info("running baseline metric group")
    before_columns = _obs_columns(sdata, args.table_key)

    LOGGER.info("running baseline metric: num_cells")
    num_cells = st.bl.num_cells()
    LOGGER.info("running baseline metric: perc_unassigned_transcripts")
    perc_unassigned = st.bl.perc_unassigned_transcripts()
    LOGGER.info("running baseline metric: transcripts_per_cell")
    st.bl.transcripts_per_cell()
    LOGGER.info("running baseline metric: genes_per_cell")
    st.bl.genes_per_cell()
    LOGGER.info("running baseline metric: mean_transcripts_per_gene_per_cell")
    st.bl.mean_transcripts_per_gene_per_cell()
    LOGGER.info("running baseline metric: morphological_features")
    st.bl.morphological_features()

    _log_added_obs_columns(before_columns, _obs_columns(sdata, args.table_key), "baseline")
    _write_obs(sdata, args.table_key, f"{args.prefix}_baseline_obs.csv")
    summary = {
        "num_cells": int(num_cells),
        "perc_unassigned_transcripts": float(perc_unassigned),
    }
    summary_path = f"{args.prefix}_baseline_summary.json"
    LOGGER.info("writing baseline summary to %s", summary_path)
    with open(summary_path, "w") as fh:
        json.dump(summary, fh, indent=2)
    LOGGER.info(
        "finished baseline metric group: num_cells=%s, perc_unassigned_transcripts=%s",
        int(num_cells),
        float(perc_unassigned),
    )


def _run_region_similarity(st, sdata, args):
    LOGGER.info("running region_similarity metric group")
    before_columns = _obs_columns(sdata, args.table_key)
    nsk = _none_or_str(args.nucleus_shapes_key)
    if nsk is not None:
        LOGGER.info("running region_similarity metric: similarity_nucleus_cytoplasm")
        try:
            st.rs.similarity_nucleus_cytoplasm()
        except (AssertionError, ValueError) as exc:
            LOGGER.warning("similarity_nucleus_cytoplasm skipped: %s", exc)
            warnings.warn(f"similarity_nucleus_cytoplasm skipped: {exc}")
    else:
        LOGGER.info("skipping similarity_nucleus_cytoplasm because nucleus_shapes_key is None")
    LOGGER.info("running region_similarity metric: similarity_border_neighborhood")
    try:
        st.rs.similarity_border_neighborhood()
    except (KeyError, ValueError) as exc:
        LOGGER.warning("similarity_border_neighborhood skipped: %s", exc)
        warnings.warn(f"similarity_border_neighborhood skipped: {exc}")
    _log_added_obs_columns(before_columns, _obs_columns(sdata, args.table_key), "region_similarity")
    _write_obs(sdata, args.table_key, f"{args.prefix}_region_similarity_obs.csv")
    LOGGER.info("finished region_similarity metric group")


def _run_volume(st, sdata, args):
    LOGGER.info("running volume metric group")
    before_columns = _obs_columns(sdata, args.table_key)
    if _none_or_str(args.points_z_key) is not None:
        LOGGER.info("running volume metric: similarity_top_bottom")
        st.vl.similarity_top_bottom(q=args.z_quantile)
    else:
        LOGGER.info("skipping similarity_top_bottom because points_z_key is None")
    _log_added_obs_columns(before_columns, _obs_columns(sdata, args.table_key), "volume")
    _write_obs(sdata, args.table_key, f"{args.prefix}_volume_similarity_obs.csv")
    LOGGER.info("finished volume metric group")


def _run_clustering_stability(st, sdata, args):
    LOGGER.info("running clustering_stability metric group")
    before_columns = _obs_columns(sdata, args.table_key)
    adata = sdata.tables[args.table_key]
    x = adata.X
    cell_totals = np.asarray(x.sum(axis=1)).flatten() if hasattr(x, "sum") else np.array(x).sum(axis=1)
    expressed_mask = cell_totals > 0
    n_expressed = int(expressed_mask.sum())
    LOGGER.info(
        "clustering_stability found %s/%s cells with non-zero expression",
        n_expressed,
        adata.n_obs,
    )

    connectedness = silhouette = ari = None

    if n_expressed < 10:
        LOGGER.warning("only %s cells have non-zero expression; skipping clustering", n_expressed)
    else:
        if n_expressed < adata.n_obs:
            LOGGER.info("filtering clustering_stability input to non-zero-expression cells")
            sdata.tables[args.table_key] = adata[expressed_mask].copy()
        adata = sdata.tables[args.table_key]

        LOGGER.info("normalizing and log-transforming clustering_stability input")
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        LOGGER.info("computing PCA, neighbors, UMAP, and Leiden clustering")
        sc.pp.pca(adata, n_comps=min(args.n_pcs, adata.n_vars - 1, adata.n_obs - 1))
        sc.pp.neighbors(adata, n_neighbors=min(args.n_neighbors, adata.n_obs - 1))
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=args.leiden_resolution, flavor="igraph", n_iterations=2)

        LOGGER.info("running clustering_stability metric: cluster_connectedness")
        connectedness = st.cs.cluster_connectedness(use_weights=True, leiden_kwargs={"flavor": "igraph"})
        LOGGER.info("running clustering_stability metric: silhouette_score")
        silhouette = st.cs.silhouette_score(leiden_kwargs={"flavor": "igraph"})
        LOGGER.info("running clustering_stability metric: adjusted_rand_index")
        ari = st.cs.adjusted_rand_index(leiden_kwargs={"flavor": "igraph"})
        LOGGER.info("running clustering_stability metric: purity")
        st.cs.purity(leiden_kwargs={"flavor": "igraph"})

    _log_added_obs_columns(before_columns, _obs_columns(sdata, args.table_key), "clustering_stability")
    _write_obs(sdata, args.table_key, f"{args.prefix}_clustering_stability_obs.csv")
    metrics = {
        "cluster_connectedness": float(connectedness) if connectedness is not None else None,
        "silhouette_score": float(silhouette) if silhouette is not None else None,
        "adjusted_rand_index": float(ari) if ari is not None else None,
    }
    metrics_path = f"{args.prefix}_clustering_stability_metrics.json"
    LOGGER.info("writing clustering_stability summary metrics to %s", metrics_path)
    with open(metrics_path, "w") as fh:
        json.dump(metrics, fh, indent=2)
    LOGGER.info("finished clustering_stability metric group")


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    args = parse_args()
    metrics = _parse_metrics(args.metrics)
    if not metrics:
        raise SystemExit("No SegTraQ metric groups selected.")

    LOGGER.info("selected SegTraQ metric group(s): %s", ", ".join(metrics))
    LOGGER.info("reading Xenium bundle from %s", args.xenium_dir)
    sdata = spatialdata_io.xenium(
        args.xenium_dir,
        cells_as_circles=False,
        cells_labels=False,
        nucleus_labels=False,
        morphology_mip=False,
    )
    LOGGER.info("finished reading Xenium bundle")

    LOGGER.info("building SegTraQ object")
    st = _build_segtraq(sdata, args)
    LOGGER.info("filtering control and low-quality transcripts with min_qv=%s", args.min_qv)
    st.filter_control_and_low_quality_transcripts(min_qv=_parse_min_qv(args.min_qv))

    if "baseline" in metrics:
        _run_baseline(st, sdata, args)
    if "region_similarity" in metrics:
        _run_region_similarity(st, sdata, args)
    if "volume" in metrics:
        _run_volume(st, sdata, args)
    if "clustering_stability" in metrics:
        _run_clustering_stability(st, sdata, args)
    LOGGER.info("finished SegTraQ QC")


if __name__ == "__main__":
    main()
