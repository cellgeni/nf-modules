#!/usr/bin/env python3
"""Export a spatialdata-io dataset to SpatialData .zarr format."""

import argparse
import inspect
import json
import logging
import os
import sys
import tempfile
from pathlib import Path

try:
    import spatialdata_io as sdio
except ImportError as exc:
    raise SystemExit("ERROR: spatialdata-io is required.") from exc


LOGGER = logging.getLogger(__name__)

READERS = {
    name: getattr(sdio, name)
    for name in (
        "codex",
        "cosmx",
        "curio",
        "dbit",
        "mcmicro",
        "merscope",
        "seqfish",
        "steinbock",
        "stereoseq",
        "visium",
        "visium_hd",
        "xenium",
        "generic",
        "image",
        "geojson",
    )
    if hasattr(sdio, name)
}
if hasattr(sdio, "experimental") and hasattr(sdio.experimental, "iss"):
    READERS["iss"] = sdio.experimental.iss

ALIASES = {
    k.lower(): v
    for k, v in {
        "phenocycler": "codex",
        "codex": "codex",
        "dbit-seq": "dbit",
        "dbit_seq": "dbit",
        "visiumhd": "visium_hd",
        "visium-hd": "visium_hd",
        "visium_hd": "visium_hd",
        "geneps": "seqfish",
        "genep": "seqfish",
        "seqfish": "seqfish",
        "stereo-seq": "stereoseq",
        "stereoseq": "stereoseq",
        "xenium": "xenium",
        "cosmx": "cosmx",
        "merscope": "merscope",
        "mcmicro": "mcmicro",
        "iss": "iss",
        "curio": "curio",
        "generic": "generic",
        "image": "image",
        "geojson": "geojson",
    }.items()
    if v in READERS
}

AUTO_ORDER = [
    "visium",
    "visium_hd",
    "xenium",
    "codex",
    "cosmx",
    "curio",
    "dbit",
    "iss",
    "mcmicro",
    "merscope",
    "seqfish",
    "steinbock",
    "stereoseq",
]

XENIUM_PARQUET_SCHEMAS = {
    "nucleus_boundaries.parquet": {
        "reader_flag": "nucleus_boundaries",
        "cell_id": "string",
        "vertex_x": "float32",
        "vertex_y": "float32",
        "label_id": "int64",
    },
    "cell_boundaries.parquet": {
        "reader_flag": "cells_boundaries",
        "cell_id": "string",
        "vertex_x": "float32",
        "vertex_y": "float32",
        "label_id": "int64",
    },
}


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
    filtered = (
        kwargs
        if has_varkw
        else {k: v for k, v in kwargs.items() if k in sig.parameters}
    )
    return func(str(path), **filtered)


def _prepare_xenium_bundle(input_path, kwargs):
    if not input_path.is_dir():
        return input_path, None

    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError:
        LOGGER.warning(
            "pyarrow is unavailable; skipping Xenium parquet schema preflight."
        )
        return input_path, None

    parquet_fixes = []
    for filename, expected_schema in XENIUM_PARQUET_SCHEMAS.items():
        parquet_path = input_path / filename
        if not parquet_path.exists():
            continue
        metadata = pq.ParquetFile(parquet_path).metadata
        reader_flag = expected_schema["reader_flag"]
        if metadata.num_rows == 0:
            if kwargs.get(reader_flag, True):
                kwargs[reader_flag] = False
                LOGGER.warning(
                    "Disabled Xenium reader option %s because %s is empty.",
                    reader_flag,
                    filename,
                )
            continue
        schema = pq.read_schema(parquet_path)
        null_columns = [
            field.name for field in schema if pa.types.is_null(field.type)
        ]
        if null_columns:
            parquet_fixes.append((filename, expected_schema, null_columns))

    if not parquet_fixes:
        return input_path, None

    temp_dir = tempfile.TemporaryDirectory(
        prefix=f".{input_path.name}.spatialdata_export_", dir="."
    )
    sanitized_path = Path(temp_dir.name)

    for child in input_path.iterdir():
        os.symlink(child.resolve(), sanitized_path / child.name)

    for filename, expected_schema, null_columns in parquet_fixes:
        source = input_path / filename
        target = sanitized_path / filename
        target.unlink()
        table = pq.read_table(source)
        schema = pa.schema(
            [
                pa.field(name, getattr(pa, type_name)())
                for name, type_name in expected_schema.items()
                if name != "reader_flag"
            ]
        )
        table = table.select([field.name for field in schema]).cast(schema)
        pq.write_table(table, target)
        LOGGER.warning(
            "Rewrote %s in a local sanitized Xenium bundle because columns %s "
            "had Arrow null type.",
            filename,
            ", ".join(null_columns),
        )

    return sanitized_path, temp_dir


def _load_sdata(input_path, reader_name, kwargs):
    if reader_name != "auto":
        if reader_name not in READERS:
            raise ValueError(
                f"Unknown reader '{reader_name}'. Known: {', '.join(sorted(READERS))}"
            )
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
    parser = argparse.ArgumentParser(
        description="Export spatialdata-io datasets to .zarr."
    )
    parser.add_argument("input_path")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument(
        "--reader",
        default="xenium",
        help="Reader name (default: xenium). Use --list-readers to see options.",
    )
    parser.add_argument(
        "--reader-arg",
        action="append",
        metavar="KEY=VALUE",
        help="Reader kwargs as KEY=VALUE (repeatable).",
    )
    parser.add_argument(
        "--reader-kwargs-json", help="JSON string or file path with reader kwargs."
    )
    parser.add_argument("-@", "--threads", type=int, default=None)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--list-readers", action="store_true")
    parser.add_argument("--verbose", action="store_true")
    parser.add_argument(
        "--version", action="version", version="spatialdata_export.py 0.1.0"
    )
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

    temp_bundle = None
    try:
        if reader_name == "xenium":
            input_path, temp_bundle = _prepare_xenium_bundle(input_path, kwargs)

        sdata = _load_sdata(input_path, reader_name, kwargs)
    except Exception as exc:
        LOGGER.error("Failed to load: %s", exc)
        return 1

    try:
        LOGGER.info("Writing to %s", output_path)
        sdata.write(str(output_path), overwrite=args.overwrite)
    except Exception as exc:
        LOGGER.error("Failed to write zarr: %s", exc)
        return 1
    finally:
        if temp_bundle is not None:
            temp_bundle.cleanup()

    return 0


if __name__ == "__main__":
    sys.exit(main())
