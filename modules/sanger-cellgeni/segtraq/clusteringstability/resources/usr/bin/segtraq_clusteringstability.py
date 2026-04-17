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
    parser.add_argument("--images_key", default="image")
    parser.add_argument("--centroid_x_key", default="x_centroid")
    parser.add_argument("--centroid_y_key", default="y_centroid")
    parser.add_argument("--leiden_resolution", type=float, default=0.5)
    parser.add_argument("--n_pcs", type=int, default=50)
    parser.add_argument("--n_neighbors", type=int, default=15)
    return parser.parse_args()


def main():
    args = parse_args()

    sdata = spatialdata.read_zarr(args.zarr_dir)
    st = segtraq.SegTraQ(
        sdata,
        images_key=args.images_key,
        tables_centroid_x_key=args.centroid_x_key,
        tables_centroid_y_key=args.centroid_y_key,
    )

    st.filter_control_and_low_quality_transcripts()

    table_key = next(iter(sdata.tables.keys()))
    adata = sdata.tables[table_key]

    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.pp.pca(adata, n_comps=min(args.n_pcs, adata.n_vars - 1))
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors)
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
