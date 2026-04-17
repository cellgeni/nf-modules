#!/usr/bin/env python3
import argparse

import spatialdata
import segtraq


def parse_args():
    parser = argparse.ArgumentParser(description="SegTraQ 3D volume similarity metrics for 10x Xenium")
    parser.add_argument("--zarr_dir", required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--z_quantile", type=float, default=0.30,
                        help="Quantile splitting z-axis into top and bottom halves (default: 0.30)")
    parser.add_argument("--images_key", default="morphology_focus", help="Key for image in sdata.images (XeniumKeys.MORPHOLOGY_FOCUS_FILE)")
    parser.add_argument("--table_key", default="table", help="Key for cell table in sdata.tables")
    parser.add_argument("--centroid_x_key", default="x_centroid", help="XeniumKeys.CELL_X")
    parser.add_argument("--centroid_y_key", default="y_centroid", help="XeniumKeys.CELL_Y")
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
    st.vl.similarity_top_bottom(q=args.z_quantile)

    sdata.tables[args.table_key].obs.to_csv(f"{args.prefix}_volume_similarity_obs.csv")


if __name__ == "__main__":
    main()
