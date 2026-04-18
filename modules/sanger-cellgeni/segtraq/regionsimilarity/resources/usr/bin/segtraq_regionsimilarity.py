#!/usr/bin/env python3
import argparse
import warnings

import spatialdata
import segtraq


def parse_args():
    parser = argparse.ArgumentParser(description="SegTraQ region similarity metrics for 10x Xenium")
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

    if nsk is not None:
        st.rs.similarity_nucleus_cytoplasm()
    try:
        st.rs.similarity_border_neighborhood()
    except (KeyError, ValueError) as e:
        warnings.warn(f"similarity_border_neighborhood skipped: {e}")

    sdata.tables[args.table_key].obs.to_csv(f"{args.prefix}_region_similarity_obs.csv")


if __name__ == "__main__":
    main()
