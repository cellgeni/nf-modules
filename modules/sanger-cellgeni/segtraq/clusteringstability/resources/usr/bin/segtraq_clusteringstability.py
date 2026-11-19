#!/usr/bin/env python3
import argparse
import json

import scanpy as sc
import spatialdata
import segtraq


def parse_args():
    parser = argparse.ArgumentParser(description="SegTraQ clustering stability metrics for 10x Xenium")
    parser.add_argument("--zarr_dir", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--images_key", default="morphology_focus")
    parser.add_argument("--table_key", default="table")
    parser.add_argument("--tables_cell_id_key", default="cell_id")
    parser.add_argument("--tables_area_key", default="cell_area",
                        help="Obs column for cell area; use 'None' to auto-compute")
    parser.add_argument("--centroid_x_key", default=None)
    parser.add_argument("--centroid_y_key", default=None)
    parser.add_argument("--shapes_key", default="cell_boundaries")
    parser.add_argument("--shapes_cell_id_key", default="cell_id")
    parser.add_argument("--nucleus_shapes_key", default="nucleus_boundaries",
                        help="Key for nucleus shapes; use 'None' if unavailable")
    parser.add_argument("--nucleus_shapes_cell_id_key", default="cell_id")
    parser.add_argument("--points_z_key", default="z",
                        help="Z coordinate column in transcripts; use 'None' for 2D data")
    parser.add_argument("--points_gene_key", default="feature_name")
    parser.add_argument("--points_background_id", default="UNASSIGNED",
                        help="Background cell ID in transcripts (Xenium default: UNASSIGNED, use 'None' to disable)")
    parser.add_argument("--min_qv", default="20.0",
                        help="Minimum quality value for transcript filtering; use 'None' to disable")
    parser.add_argument("--leiden_resolution", type=float, default=0.5)
    parser.add_argument("--n_pcs", type=int, default=50)
    parser.add_argument("--n_neighbors", type=int, default=15)
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


def _fix_shapes_index(sdata, shapes_key, cell_id_key):
    if shapes_key not in sdata.shapes:
        return
    shape = sdata.shapes[shapes_key]
    if cell_id_key not in shape.columns and shape.index.name != cell_id_key:
        shape.index.name = cell_id_key
        sdata.shapes[shapes_key] = shape


def main():
    args = parse_args()
    bg_id = _parse_bg_id(args.points_background_id)
    min_qv = None if args.min_qv == "None" else float(args.min_qv)

    sdata = spatialdata.read_zarr(args.zarr_dir)

    _fix_shapes_index(sdata, args.shapes_key, args.shapes_cell_id_key)
    nsk = _none_or_str(args.nucleus_shapes_key)
    if nsk:
        _fix_shapes_index(sdata, nsk, args.nucleus_shapes_cell_id_key)

    st = segtraq.SegTraQ(
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
        points_background_id=bg_id,
    )

    st.filter_control_and_low_quality_transcripts(min_qv=min_qv)

    import sys
    import numpy as np

    adata = sdata.tables[args.table_key]
    X = adata.X
    cell_totals = np.asarray(X.sum(axis=1)).flatten() if hasattr(X, "sum") else np.array(X).sum(axis=1)
    expressed_mask = cell_totals > 0
    n_expressed = int(expressed_mask.sum())

    connectedness = silhouette = ari = None

    if n_expressed < 10:
        print(f"WARNING: only {n_expressed} cells have non-zero expression; skipping clustering.", file=sys.stderr)
    else:
        if n_expressed < adata.n_obs:
            sdata.tables[args.table_key] = adata[expressed_mask].copy()
        adata = sdata.tables[args.table_key]

        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)
        sc.pp.pca(adata, n_comps=min(args.n_pcs, adata.n_vars - 1, adata.n_obs - 1))
        sc.pp.neighbors(adata, n_neighbors=min(args.n_neighbors, adata.n_obs - 1))
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=args.leiden_resolution, flavor="igraph", n_iterations=2)

        connectedness = st.cs.cluster_connectedness()
        silhouette = st.cs.silhouette_score()
        ari = st.cs.adjusted_rand_index()
        st.cs.purity()

    adata.obs.to_csv(f"{args.prefix}_clustering_stability_obs.csv")

    metrics = {
        "cluster_connectedness": float(connectedness) if connectedness is not None else None,
        "silhouette_score": float(silhouette) if silhouette is not None else None,
        "adjusted_rand_index": float(ari) if ari is not None else None,
    }
    with open(f"{args.prefix}_clustering_stability_metrics.json", "w") as fh:
        json.dump(metrics, fh, indent=2)


if __name__ == "__main__":
    main()
