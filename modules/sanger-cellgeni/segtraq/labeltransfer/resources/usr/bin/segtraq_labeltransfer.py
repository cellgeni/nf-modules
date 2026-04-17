#!/usr/bin/env python3
import argparse

import anndata as ad
import spatialdata
import segtraq


def parse_args():
    parser = argparse.ArgumentParser(description="SegTraQ label transfer from reference scRNA-seq for 10x Xenium")
    parser.add_argument("--zarr_dir", required=True)
    parser.add_argument("--ref_h5ad", required=True, help="Reference scRNA-seq AnnData (H5AD)")
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--cell_type_key", default="celltype_major", help="Column in ref_h5ad.obs with cell type labels")
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

    adata_ref = ad.read_h5ad(args.ref_h5ad)
    st.run_label_transfer(adata_ref, ref_cell_type=args.cell_type_key)

    table_key = next(iter(sdata.tables.keys()))
    labeled = sdata.tables[table_key]
    labeled.write_h5ad(f"{args.prefix}_labeled.h5ad")


if __name__ == "__main__":
    main()
