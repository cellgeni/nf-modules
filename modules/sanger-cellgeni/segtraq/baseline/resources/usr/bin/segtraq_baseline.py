#!/usr/bin/env python3
import argparse
import json

import spatialdata
import segtraq


def parse_args():
    parser = argparse.ArgumentParser(description="SegTraQ baseline QC metrics for 10x Xenium")
    parser.add_argument("--zarr_dir", required=True, help="Path to SpatialData zarr directory")
    parser.add_argument("--prefix", required=True, help="Output file prefix")
    parser.add_argument("--images_key", default="morphology_focus", help="Key for image in sdata.images (XeniumKeys.MORPHOLOGY_FOCUS_FILE)")
    parser.add_argument("--table_key", default="table", help="Key for cell table in sdata.tables")
    parser.add_argument("--centroid_x_key", default="x_centroid", help="Column name for centroid x in table obs (XeniumKeys.CELL_X)")
    parser.add_argument("--centroid_y_key", default="y_centroid", help="Column name for centroid y in table obs (XeniumKeys.CELL_Y)")
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

    num_cells = st.bl.num_cells()
    perc_unassigned = st.bl.perc_unassigned_transcripts()
    st.bl.transcripts_per_cell()
    st.bl.genes_per_cell()
    st.bl.mean_transcripts_per_gene_per_cell()
    st.bl.morphological_features()

    obs = sdata.tables[args.table_key].obs
    obs.to_csv(f"{args.prefix}_baseline_obs.csv")

    summary = {
        "num_cells": int(num_cells),
        "perc_unassigned_transcripts": float(perc_unassigned),
    }
    with open(f"{args.prefix}_baseline_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)


if __name__ == "__main__":
    main()
