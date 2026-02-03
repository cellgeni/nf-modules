#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Export SpatialData-IO reader outputs to an AnnData .h5ad file."""

from __future__ import annotations

import argparse
import inspect
import json
import logging
import sys
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Mapping, Tuple

try:
    import spatialdata_io as sdio
    from spatialdata_io import experimental as sdio_experimental
except ImportError as exc:  # pragma: no cover - import guard for runtime
    raise SystemExit(
        "ERROR: spatialdata-io is required to run this script. "
        "Please install spatialdata-io in the environment."
    ) from exc
if not hasattr(sdio_experimental, "to_legacy_anndata"):
    raise SystemExit(
        "ERROR: spatialdata-io does not expose experimental.to_legacy_anndata. "
        "Please install a compatible spatialdata-io version."
    )


LOGGER = logging.getLogger("spatialdata_export")


def _build_readers() -> Dict[str, Callable[..., Any]]:
    readers: Dict[str, Callable[..., Any]] = {}
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
    ):
        if hasattr(sdio, name):
            readers[name] = getattr(sdio, name)
    if hasattr(sdio_experimental, "iss"):
        readers["iss"] = sdio_experimental.iss
    return readers


def _build_aliases(readers: Dict[str, Callable[..., Any]]) -> Dict[str, str]:
    aliases = {
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
        "experimental.iss": "iss",
        "curio": "curio",
        "generic": "generic",
        "image": "image",
        "geojson": "geojson",
    }
    return {key.lower(): value for key, value in aliases.items() if value in readers}


def _setup_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def _select_element(
    elements: Mapping[str, Any], key: str | None, kind: str
) -> Tuple[str, Any]:
    if not elements:
        raise ValueError(f"No {kind} available to export.")
    if key:
        if key not in elements:
            available = ", ".join(sorted(elements.keys()))
            raise ValueError(f"{kind.title()} key '{key}' not found. Available: {available}")
        return key, elements[key]
    first_key = sorted(elements.keys())[0]
    LOGGER.info("No %s key provided; using '%s'.", kind, first_key)
    return first_key, elements[first_key]


def _prepare_array(data: Any) -> Any:
    array = data
    dims = getattr(data, "dims", None)
    if dims:
        if "y" in dims and "x" in dims:
            channel_dims = [d for d in dims if d not in ("y", "x")]
            if channel_dims:
                try:
                    data = data.transpose("y", "x", *channel_dims)
                except Exception:  # pragma: no cover - best effort
                    data = data
            array = data
    array = getattr(array, "data", array)
    if hasattr(array, "compute"):
        array = array.compute()
    try:
        import numpy as np
    except ImportError as exc:  # pragma: no cover - import guard for runtime
        raise SystemExit("ERROR: numpy is required to export images.") from exc
    np_array = np.asarray(array)
    if np_array.ndim == 3 and np_array.shape[0] <= 10 and np_array.shape[1] > 32 and np_array.shape[2] > 32:
        np_array = np.transpose(np_array, (1, 2, 0))
    return np_array


def _write_tiff(output_path: Path, array: Any) -> None:
    try:
        import tifffile
    except ImportError:
        tifffile = None
    if tifffile is not None:
        tifffile.imwrite(str(output_path), array)
        return
    try:
        import imageio.v2 as imageio
    except ImportError as exc:  # pragma: no cover - import guard for runtime
        raise SystemExit(
            "ERROR: tifffile or imageio is required to export images."
        ) from exc
    imageio.imwrite(str(output_path), array)


def _parse_value(raw: str) -> Any:
    text = raw.strip()
    lower = text.lower()
    if lower in {"none", "null"}:
        return None
    if lower in {"true", "false"}:
        return lower == "true"
    if text and text[0] in "[{":
        try:
            return json.loads(text)
        except json.JSONDecodeError:
            return text
    try:
        return int(text)
    except ValueError:
        pass
    try:
        return float(text)
    except ValueError:
        return text


def _parse_kv_pairs(pairs: Iterable[str]) -> Dict[str, Any]:
    parsed: Dict[str, Any] = {}
    for pair in pairs:
        if "=" not in pair:
            raise ValueError(f"Invalid --reader-arg value '{pair}'. Expected KEY=VALUE.")
        key, value = pair.split("=", 1)
        parsed[key] = _parse_value(value)
    return parsed


def _load_reader_kwargs(args: argparse.Namespace) -> Dict[str, Any]:
    kwargs: Dict[str, Any] = {}
    if args.reader_kwargs_json:
        kwargs_path = Path(args.reader_kwargs_json)
        try:
            if kwargs_path.exists():
                kwargs = json.loads(kwargs_path.read_text())
            else:
                kwargs = json.loads(args.reader_kwargs_json)
        except json.JSONDecodeError as exc:
            raise ValueError("Failed to parse --reader-kwargs-json as JSON.") from exc
        if not isinstance(kwargs, dict):
            raise ValueError("--reader-kwargs-json must resolve to a JSON object.")
    kwargs.update(_parse_kv_pairs(args.reader_arg or []))
    if args.threads is not None and "n_jobs" not in kwargs:
        kwargs["n_jobs"] = args.threads
    return kwargs


def _filter_kwargs(func: Callable[..., Any], kwargs: Dict[str, Any]) -> Tuple[Dict[str, Any], List[str]]:
    signature = inspect.signature(func)
    accepts_kwargs = any(param.kind == param.VAR_KEYWORD for param in signature.parameters.values())
    if accepts_kwargs:
        return dict(kwargs), []
    accepted = {k: v for k, v in kwargs.items() if k in signature.parameters}
    rejected = [k for k in kwargs if k not in signature.parameters]
    return accepted, rejected


def _auto_readers(readers: Dict[str, Callable[..., Any]]) -> List[str]:
    preferred = [
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
    available = [name for name in preferred if name in readers]
    extras = [name for name in readers if name not in available]
    return available + extras


def _load_sdata(
    input_path: Path,
    reader_name: str,
    readers: Dict[str, Callable[..., Any]],
    reader_kwargs: Dict[str, Any],
) -> Tuple[Any, str]:
    if reader_name == "auto":
        errors: List[str] = []
        for name in _auto_readers(readers):
            func = readers[name]
            filtered_kwargs, rejected = _filter_kwargs(func, reader_kwargs)
            if rejected:
                LOGGER.debug("Ignoring unsupported kwargs for reader '%s': %s", name, ", ".join(rejected))
            try:
                LOGGER.info("Trying reader '%s'...", name)
                return func(str(input_path), **filtered_kwargs), name
            except Exception as exc:  # pragma: no cover - best-effort auto detection
                errors.append(f"{name}: {exc}")
                LOGGER.debug("Reader '%s' failed: %s", name, exc, exc_info=True)
        error_summary = "\n  - ".join(errors) if errors else "No readers available."
        raise RuntimeError(f"Auto-detection failed. Errors:\n  - {error_summary}")

    if reader_name not in readers:
        raise ValueError(f"Unknown reader '{reader_name}'. Use --list-readers to see options.")

    func = readers[reader_name]
    filtered_kwargs, rejected = _filter_kwargs(func, reader_kwargs)
    if rejected:
        LOGGER.warning("Ignoring unsupported kwargs for reader '%s': %s", reader_name, ", ".join(rejected))
    LOGGER.info("Using reader '%s'.", reader_name)
    return func(str(input_path), **filtered_kwargs), reader_name


def _build_parser(readers: Dict[str, Callable[..., Any]]) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Export SpatialData-IO supported datasets to a legacy AnnData .h5ad file."
    )
    parser.add_argument("input_path", help="Path to the spatial technology output folder or file.")
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output .h5ad file path.",
    )
    parser.add_argument(
        "--reader",
        default="auto",
        help=(
            "Spatialdata-io reader to use (default: auto). "
            "Use --list-readers to see options."
        ),
    )
    parser.add_argument(
        "--reader-arg",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help="Additional reader kwargs as key=value (can be repeated).",
    )
    parser.add_argument(
        "--reader-kwargs-json",
        help="JSON string or path to JSON file with reader kwargs.",
    )
    parser.add_argument(
        "--coordinate-system",
        default=None,
        help="Coordinate system passed to to_legacy_anndata.",
    )
    parser.add_argument(
        "--table-name",
        default='table',
        help="Table name passed to to_legacy_anndata.",
    )
    parser.add_argument(
        "--include-images",
        action="store_false",
        help="Include downscaled images in the legacy AnnData output.",
    )
    parser.add_argument(
        "--raw-image-output",
        default=None,
        help="Optional output path for exporting a raw image as TIFF.",
    )
    parser.add_argument(
        "--raw-image-key",
        default=None,
        help="Image key within the SpatialData object to export as raw image.",
    )
    parser.add_argument(
        "--label-image-output",
        default=None,
        help="Optional output path for exporting a label image as TIFF.",
    )
    parser.add_argument(
        "--label-image-key",
        default=None,
        help="Label key within the SpatialData object to export as label image.",
    )
    parser.add_argument(
        "-@",
        "--threads",
        type=int,
        default=None,
        help="Threads to pass to readers that support n_jobs.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite the output file if it exists.",
    )
    parser.add_argument(
        "--list-readers",
        action="store_true",
        help="List available readers and exit.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Enable verbose logging.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="spatialdata_export.py 0.1.0",
    )
    return parser


def main(argv: List[str] | None = None) -> int:
    readers = _build_readers()
    parser = _build_parser(readers)
    args = parser.parse_args(argv)
    _setup_logging(args.verbose)

    if args.list_readers:
        for name in sorted(readers.keys()):
            print(name)
        return 0

    input_path = Path(args.input_path)
    if not input_path.exists():
        LOGGER.error("Input path does not exist: %s", input_path)
        return 2

    output_path = Path(args.output)
    if output_path.exists() and not args.overwrite:
        LOGGER.error("Output file already exists: %s (use --overwrite to replace)", output_path)
        return 2
    output_path.parent.mkdir(parents=True, exist_ok=True)

    aliases = _build_aliases(readers)
    reader_key = args.reader.strip()
    reader_name = aliases.get(reader_key.lower(), reader_key)

    try:
        reader_kwargs = _load_reader_kwargs(args)
    except ValueError as exc:
        LOGGER.error("%s", exc)
        return 2

    try:
        sdata, used_reader = _load_sdata(input_path, reader_name, readers, reader_kwargs)
    except Exception as exc:
        LOGGER.error("Failed to load dataset with reader '%s': %s", reader_name, exc)
        return 1

    if args.raw_image_output:
        try:
            images = getattr(sdata, "images", {})
            key, raw_image = _select_element(images, args.raw_image_key, "raw image")
            raw_array = _prepare_array(raw_image)
            raw_output = Path(args.raw_image_output)
            raw_output.parent.mkdir(parents=True, exist_ok=True)
            LOGGER.info("Writing raw image '%s' to %s", key, raw_output)
            _write_tiff(raw_output, raw_array)
        except Exception as exc:
            LOGGER.error("Failed to export raw image: %s", exc)
            return 1

    if args.label_image_output:
        try:
            labels = getattr(sdata, "labels", {})
            key, label_image = _select_element(labels, args.label_image_key, "label image")
            label_array = _prepare_array(label_image)
            label_output = Path(args.label_image_output)
            label_output.parent.mkdir(parents=True, exist_ok=True)
            LOGGER.info("Writing label image '%s' to %s", key, label_output)
            _write_tiff(label_output, label_array)
        except Exception as exc:
            LOGGER.error("Failed to export label image: %s", exc)
            return 1

    LOGGER.info("Converting SpatialData to legacy AnnData...")
    try:
        adata = sdio_experimental.to_legacy_anndata(
            sdata,
            coordinate_system=args.coordinate_system,
            table_name=args.table_name,
            include_images=args.include_images,
        )
    except Exception as exc:
        LOGGER.error("Failed to convert to legacy AnnData: %s", exc)
        adata = sdata.tables[args.table_name]
        # return 1
    

    LOGGER.info("Writing AnnData to %s", output_path)
    try:
        adata.write_h5ad(output_path)
    except Exception as exc:
        LOGGER.error("Failed to write .h5ad file: %s", exc)
        return 1

    LOGGER.info("Done (reader: %s).", used_reader)
    return 0


if __name__ == "__main__":
    sys.exit(main())
