#!/usr/bin/env python3
"""Download OMERO table annotations to local CSV files."""

from __future__ import annotations

import argparse
import getpass
import os
from pathlib import Path

os.environ.setdefault("OMERO_TMPDIR", "/tmp")
os.environ.setdefault("TMPDIR", "/tmp")

import omero2pandas
from omero.gateway import BlitzGateway

DEFAULT_HOST = "wsi-omero-prod-02.internal.sanger.ac.uk"


def parse_args() -> argparse.Namespace:
    tool_version = getattr(omero2pandas, "__version__", "unknown")
    p = argparse.ArgumentParser(
        description=(
            "Download OMERO table annotations containing spot locations. "
            "Prompts for credentials unless --user/--password or env vars are set."
        )
    )
    p.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {tool_version}",
    )
    p.add_argument("--host", default=DEFAULT_HOST)
    p.add_argument("--port", type=int, default=4064)
    p.add_argument(
        "-u",
        "--user",
        default=None,
        help="OMERO username (overrides OMERO_USER env var).",
    )
    p.add_argument(
        "-P",
        "--password",
        default=None,
        help="OMERO password (overrides OMERO_PASSWORD env var).",
    )
    p.add_argument(
        "--annotation-id",
        type=int,
        required=True,
        help="OMERO table annotation ID to download.",
    )
    p.add_argument(
        "--linking-index",
        type=int,
        required=True,
        help=(
            "Index to attach to all rows from this annotation as metadata "
            "(e.g. channel/time/plane linking key)."
        ),
    )
    p.add_argument(
        "--out-dir",
        default="downloaded_annotations",
        help="Directory for downloaded per-annotation CSVs",
    )
    return p.parse_args()


def ensure_credentials(user: str | None, password: str | None) -> tuple[str, str]:
    username = user or os.getenv("OMERO_USER")
    secret = password or os.getenv("OMERO_PASSWORD")
    if not username:
        username = input("OMERO username: ").strip()
    if not secret:
        secret = getpass.getpass("OMERO password: ").strip()
    if not username or not secret:
        raise ValueError("OMERO credentials are required.")
    return username, secret


def resolve_image_id(conn: BlitzGateway, annotation_id: int) -> int | None:
    """Resolve source image id for an annotation id via OMERO link APIs."""
    for link in conn.getAnnotationLinks("Image", ann_ids=[annotation_id]):
        parent = link.getParent()
        if parent is not None and getattr(parent, "id", None) is not None:
            return int(parent.id)

    # Fallback for annotations linked to ROI instead of Image.
    for link in conn.getAnnotationLinks("Roi", ann_ids=[annotation_id]):
        roi = link.getParent()
        if roi is None:
            continue
        image = roi.getImage()
        if image is not None and getattr(image, "id", None) is not None:
            return int(image.id)

    return None


def main() -> None:
    args = parse_args()
    username, password = ensure_credentials(args.user, args.password)
    annotation_id = int(args.annotation_id)
    linking_index = int(args.linking_index)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    conn = BlitzGateway(
        username,
        password,
        host=args.host,
        port=args.port,
    )
    if not conn.connect():
        raise RuntimeError(
            f"Failed to connect to OMERO server '{args.host}:{args.port}' for image-id resolution."
        )
    # Query across all groups accessible to the current user.
    conn.SERVICE_OPTS.setOmeroGroup("-1")

    try:
        image_id = resolve_image_id(conn, annotation_id)

        table = omero2pandas.read_table(
            annotation_id=annotation_id,
            server=args.host,
            port=args.port,
            username=username,
            password=password,
        )
        table = table.copy()
        table["source_annotation_id"] = annotation_id
        table["source_image_id"] = image_id
        table["linking_index"] = linking_index

        out_csv = out_dir / f"annotation_{annotation_id}.csv"
        table.to_csv(out_csv, index=False)
        print(f"Downloaded annotation {annotation_id} -> {out_csv.resolve()}")
    finally:
        conn.close()


if __name__ == "__main__":
    main()
