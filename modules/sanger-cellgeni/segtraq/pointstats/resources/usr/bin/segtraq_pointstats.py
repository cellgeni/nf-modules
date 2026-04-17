#!/usr/bin/env python3
import argparse

import anndata as ad
import spatialdata
import segtraq


def parse_args():
    parser = argparse.ArgumentParser(description="SegTraQ point statistics (distance to membrane) for 10x Xenium")
    parser.add_argument("--zarr_dir", required=True)
    parser.add_argument("--labeled_h5ad", required=True, help="H5AD from segtraq/labeltransfer with transferred_cell_type in obs")
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--genes", required=True, help="Comma-separated list of genes to analyse")
    parser.add_argument("--cell_type_key", default="transferred_cell_type")
    parser.add_argument("--cell_type_query", required=True, help="Cell type name to query distance statistics for")
    parser.add_argument("--images_key", default="image")
    parser.add_argument("--centroid_x_key", default="x_centroid")
    parser.add_argument("--centroid_y_key", default="y_centroid")
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

    labeled_adata = ad.read_h5ad(args.labeled_h5ad)
    table_key = next(iter(sdata.tables.keys()))
    sdata.tables[table_key].obs[args.cell_type_key] = labeled_adata.obs[args.cell_type_key].values

    genes = [g.strip() for g in args.genes.split(",")]
    result = st.ps.distance_to_membrane(
        genes=genes,
        cell_type_key=args.cell_type_key,
        cell_type_query=args.cell_type_query,
        inplace=False,
    )

    if result is not None:
        result.to_csv(f"{args.prefix}_point_stats.csv")
    else:
        sdata.tables[table_key].obs.to_csv(f"{args.prefix}_point_stats.csv")


if __name__ == "__main__":
    main()
