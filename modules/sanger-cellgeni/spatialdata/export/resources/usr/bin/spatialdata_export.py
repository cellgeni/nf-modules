#!/usr/bin/env python3
"""Export a spatialdata-io dataset to SpatialData .zarr format."""

import argparse
import inspect
import json
import logging
import sys
from pathlib import Path

try:
    import spatialdata_io as sdio
except ImportError as exc:
    raise SystemExit("ERROR: spatialdata-io is required.") from exc

try:
    import spatialdata as sd
except ImportError as exc:
    raise SystemExit("ERROR: spatialdata is required.") from exc

LOGGER = logging.getLogger(__name__)

READERS = {
    name: getattr(sdio, name)
    for name in (
        "codex", "cosmx", "curio", "dbit", "mcmicro", "merscope",
        "seqfish", "steinbock", "stereoseq", "visium", "visium_hd",
        "xenium", "generic", "image", "geojson",
    )
    if hasattr(sdio, name)
}
if hasattr(sdio, "experimental") and hasattr(sdio.experimental, "iss"):
    READERS["iss"] = sdio.experimental.iss

ALIASES = {k.lower(): v for k, v in {
    "phenocycler": "codex", "codex": "codex",
    "dbit-seq": "dbit", "dbit_seq": "dbit",
    "visiumhd": "visium_hd", "visium-hd": "visium_hd", "visium_hd": "visium_hd",
    "geneps": "seqfish", "genep": "seqfish", "seqfish": "seqfish",
    "stereo-seq": "stereoseq", "stereoseq": "stereoseq",
    "xenium": "xenium", "cosmx": "cosmx", "merscope": "merscope",
    "mcmicro": "mcmicro", "iss": "iss", "curio": "curio",
    "generic": "generic", "image": "image", "geojson": "geojson",
}.items() if v in READERS}

AUTO_ORDER = [
    "visium", "visium_hd", "xenium", "codex", "cosmx", "curio", "dbit",
    "iss", "mcmicro", "merscope", "seqfish", "steinbock", "stereoseq",
]


def _parse_val(s):
    if s.lower() in ("none", "null"):
        return None
    if s.lower() in ("true", "false"):
        return s.lower() == "true"
    for conv in (int, float):
        try:
            return conv(s)
        except ValueError:
            pass
    return s


def _call_reader(func, path, kwargs):
    sig = inspect.signature(func)
    has_varkw = any(p.kind == p.VAR_KEYWORD for p in sig.parameters.values())
    filtered = kwargs if has_varkw else {k: v for k, v in kwargs.items() if k in sig.parameters}
    return func(str(path), **filtered)


def _load_sdata(input_path, reader_name, kwargs):
    if reader_name != "auto":
        if reader_name not in READERS:
            raise ValueError(f"Unknown reader '{reader_name}'. Known: {', '.join(sorted(READERS))}")
        LOGGER.info("Using reader '%s'.", reader_name)
        return _call_reader(READERS[reader_name], input_path, kwargs)

    for name in [n for n in AUTO_ORDER if n in READERS]:
        try:
            LOGGER.info("Trying reader '%s'...", name)
            return _call_reader(READERS[name], input_path, kwargs)
        except Exception as exc:
            LOGGER.debug("Reader '%s' failed: %s", name, exc)
    raise RuntimeError("Auto-detection failed for all readers.")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Export spatialdata-io datasets to .zarr.")
    parser.add_argument("input_path")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--reader", default="auto",
                        help="Reader name (default: auto). Use --list-readers to see options.")
    parser.add_argument("--reader-arg", action="append", metavar="KEY=VALUE",
                        help="Reader kwargs as KEY=VALUE (repeatable).")
    parser.add_argument("--reader-kwargs-json",
                        help="JSON string or file path with reader kwargs.")
    parser.add_argument("--include-points", action="store_true",
                        help="Include points/transcripts in output (excluded by default).")
    parser.add_argument("--rename-image-to", default="raw_image",
                        help="Rename a single image element to this key (default: raw_image). "
                             "Pass empty string to skip renaming.")
    parser.add_argument("-@", "--threads", type=int, default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--list-readers", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument("--version", action="version", version="spatialdata_export.py 0.1.0")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    if args.list_readers:
        print("\n".join(sorted(READERS)))
        return 0

    input_path = Path(args.input_path)
    if not input_path.exists():
        LOGGER.error("Input path does not exist: %s", input_path)
        return 2

    output_path = Path(args.output)
    if output_path.exists() and not args.overwrite:
        LOGGER.error("Output already exists: %s (use --overwrite)", output_path)
        return 2

    reader_name = ALIASES.get(args.reader.strip().lower(), args.reader.strip())

    kwargs = {}
    if args.reader_kwargs_json:
        p = Path(args.reader_kwargs_json)
        raw = p.read_text() if p.exists() else args.reader_kwargs_json
        kwargs = json.loads(raw)
    for pair in args.reader_arg or []:
        k, _, v = pair.partition("=")
        kwargs[k] = _parse_val(v)
    if args.threads is not None:
        kwargs.setdefault("n_jobs", args.threads)

    try:
        sdata = _load_sdata(input_path, reader_name, kwargs)
    except Exception as exc:
        LOGGER.error("Failed to load: %s", exc)
        return 1

    images = dict(sdata.images)
    if args.rename_image_to and len(images) == 1:
        (old_key, img), = images.items()
        if old_key != args.rename_image_to:
            LOGGER.info("Renaming image '%s' -> '%s'.", old_key, args.rename_image_to)
            images = {args.rename_image_to: img}

    points = dict(sdata.points) if args.include_points else {}
    if sdata.points and not args.include_points:
        LOGGER.info("Dropping %d points element(s) (use --include-points to keep).", len(sdata.points))

    sdata = sd.SpatialData(
        images=images,
        labels=dict(sdata.labels),
        shapes=dict(sdata.shapes),
        tables=dict(sdata.tables),
        points=points,
    )

    LOGGER.info("Writing to %s", output_path)
    try:
        sdata.write(str(output_path), overwrite=args.overwrite)
    except Exception as exc:
        LOGGER.error("Failed to write zarr: %s", exc)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
