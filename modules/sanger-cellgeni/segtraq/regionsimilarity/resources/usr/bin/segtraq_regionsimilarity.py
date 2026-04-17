#!/usr/bin/env python3
import argparse

import spatialdata
import segtraq


def parse_args():
    parser = argparse.ArgumentParser(description="SegTraQ region similarity metrics for 10x Xenium")
    parser.add_argument("--zarr_dir", required=True)
    parser.add_argument("--prefix", required=True)
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

    st.rs.similarity_nucleus_cytoplasm()
    st.rs.similarity_border_neighborhood()

    table_key = next(iter(sdata.tables.keys()))
    sdata.tables[table_key].obs.to_csv(f"{args.prefix}_region_similarity_obs.csv")


if __name__ == "__main__":
    main()
