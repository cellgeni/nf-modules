#!/usr/bin/env python3
import argparse
import csv
import os
import omero2pandas


def parse_args():
    p = argparse.ArgumentParser(description="Attach a CSV file to an OMERO image")
    p.add_argument("--host", required=True)
    p.add_argument("--port", type=int, default=4064)
    p.add_argument("--user", required=True)
    p.add_argument("--password", required=True)
    p.add_argument("--image-id", type=int, required=True)
    p.add_argument("--roi-id", type=int, required=True)
    p.add_argument("--csv", required=True, help="Path to CSV file")
    p.add_argument("--name", default=None, help="Attachment name override (optional)")
    p.add_argument(
        "--ann-csv",
        default=None,
        help="Path to write annotation ids CSV (default: <input>.ann.csv)",
    )
    return p.parse_args()


def main():
    os.environ.setdefault("OMERO_TMPDIR", "/tmp")
    os.environ.setdefault("TMPDIR", "/tmp")
    args = parse_args()
    csv_path = os.path.abspath(args.csv)
    if not os.path.isfile(csv_path):
        raise FileNotFoundError(csv_path)

    table_name = args.name or os.path.basename(csv_path)
    ann_id = omero2pandas.upload_table(
        source=csv_path,
        table_name=table_name,
        parent_id=args.image_id,
        parent_type="Image",
        server=args.host,
        port=args.port,
        username=args.user,
        password=args.password,
        links=[("Image", args.image_id), ("Roi", args.roi_id)],
    )
    ann_csv_path = (
        os.path.abspath(args.ann_csv)
        if args.ann_csv
        else os.path.splitext(csv_path)[0] + ".ann.csv"
    )
    with open(ann_csv_path, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["ann_id", "image_id", "roi_id"])
        writer.writerow([ann_id, args.image_id, args.roi_id])


if __name__ == "__main__":
    main()
