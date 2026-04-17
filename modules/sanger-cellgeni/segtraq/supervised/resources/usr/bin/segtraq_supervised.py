#!/usr/bin/env python3
import argparse

import anndata as ad
import pandas as pd
import spatialdata
import segtraq


def parse_args():
    parser = argparse.ArgumentParser(description="SegTraQ supervised QC metrics for 10x Xenium")
    parser.add_argument("--zarr_dir", required=True)
    parser.add_argument("--labeled_h5ad", required=True, help="H5AD from segtraq/labeltransfer with transferred_cell_type in obs")
    parser.add_argument("--ref_h5ad", required=True, help="Reference scRNA-seq AnnData used to derive markers")
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--cell_type_key", default="transferred_cell_type")
    parser.add_argument("--ref_cell_type_key", default="celltype_major")
    parser.add_argument("--min_pos_frac", type=float, default=0.3)
    parser.add_argument("--n_jobs", type=int, default=1)
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

    # Inject transferred labels from labeltransfer output
    labeled_adata = ad.read_h5ad(args.labeled_h5ad)
    table_key = next(iter(sdata.tables.keys()))
    sdata.tables[table_key].obs[args.cell_type_key] = labeled_adata.obs[args.cell_type_key].values

    # Compute markers from reference
    adata_ref = ad.read_h5ad(args.ref_h5ad)
    markers = segtraq.markers_from_reference(
        adata_ref,
        cell_type_key=args.ref_cell_type_key,
        n_jobs=args.n_jobs,
        min_pos_frac=args.min_pos_frac,
    )

    st.sp.marker_purity(cell_type_key=args.cell_type_key, markers=markers)
    st.sp.neighbor_contamination(cell_type_key=args.cell_type_key)
    st.sp.mutually_exclusive_coexpression_rate(cell_type_key=args.cell_type_key, markers=markers)

    sdata.tables[table_key].obs.to_csv(f"{args.prefix}_supervised_obs.csv")

    # Save contamination matrix from uns if available
    uns = sdata.tables[table_key].uns
    contamination_key = next((k for k in uns if "contamination" in k.lower()), None)
    if contamination_key is not None:
        pd.DataFrame(uns[contamination_key]).to_csv(f"{args.prefix}_contamination_matrix.csv")
    else:
        pd.DataFrame().to_csv(f"{args.prefix}_contamination_matrix.csv")


if __name__ == "__main__":
    main()
