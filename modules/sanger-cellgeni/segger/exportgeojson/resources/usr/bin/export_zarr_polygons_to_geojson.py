#!/usr/bin/env python3
"""Export Xenium cells.zarr.zip polygon sets to 10x-compatible GeoJSON."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

VERSION = "0.2.0"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Export polygon_sets/<N> from a Xenium-style cells.zarr.zip file "
            "as one GeoJSON FeatureCollection with Polygon geometries. "
            "Geometries are validated before the output file is written."
        )
    )
    parser.add_argument("zarr_zip", type=Path, help="Path to cells.zarr.zip")
    parser.add_argument("output_geojson", type=Path, help="Output GeoJSON path")
    parser.add_argument(
        "--polygon-set",
        default="1",
        help="Polygon set to export. For Xenium/Segger cells.zarr.zip, 1 is cells and 0 is nuclei.",
    )
    parser.add_argument(
        "--allow-truncated",
        action="store_true",
        help=(
            "Allow rows where num_vertices is greater than the stored coordinate pairs. "
            "The stored coordinates are still validated and exported."
        ),
    )
    parser.add_argument(
        "--fail-on-skipped",
        action="store_true",
        help="Exit with an error if any invalid polygon rows are skipped.",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=VERSION,
        help="Print the tool version and exit.",
    )
    return parser.parse_args()


def finite_pair(x: Any, y: Any) -> bool:
    return math.isfinite(float(x)) and math.isfinite(float(y))


def signed_area(ring: list[list[float]]) -> float:
    return 0.5 * sum(
        ring[i][0] * ring[i + 1][1] - ring[i + 1][0] * ring[i][1]
        for i in range(len(ring) - 1)
    )


def ring_from_flat_vertices(row: Any) -> tuple[list[list[float]], int]:
    """Convert a flat [x0, y0, x1, y1, ...] row to a closed GeoJSON ring."""
    coords: list[list[float]] = []
    values = row.tolist()

    for x, y in zip(values[0::2], values[1::2]):
        if not finite_pair(x, y):
            continue
        x_float = float(x)
        y_float = float(y)
        if x_float == 0.0 and y_float == 0.0:
            continue
        coords.append([x_float, y_float])

    stored_pairs = len(coords)

    cleaned: list[list[float]] = []
    for point in coords:
        if not cleaned or point != cleaned[-1]:
            cleaned.append(point)

    if len(cleaned) >= 3 and cleaned[0] != cleaned[-1]:
        cleaned.append(cleaned[0])

    return cleaned, stored_pairs


def validate_ring(ring: list[list[float]]) -> str | None:
    if len(ring) < 4:
        return "ring has fewer than 4 coordinates"
    if ring[0] != ring[-1]:
        return "ring is not closed"
    if len({tuple(point) for point in ring[:-1]}) < 3:
        return "ring has fewer than 3 unique vertices"
    if abs(signed_area(ring)) == 0.0:
        return "ring has zero signed area"
    return None


def make_geometry(ring: list[list[float]]) -> dict[str, Any]:
    return {"type": "Polygon", "coordinates": [ring]}


def validate_geojson_feature(feature: dict[str, Any]) -> str | None:
    if feature.get("type") != "Feature":
        return "feature type is not Feature"

    geometry = feature.get("geometry")
    if not isinstance(geometry, dict):
        return "missing geometry"

    geometry_type = geometry.get("type")
    coordinates = geometry.get("coordinates")

    if geometry_type != "Polygon":
        return f"unsupported geometry type: {geometry_type!r}"

    if not coordinates or not coordinates[0]:
        return "empty polygon coordinates"
    return validate_ring(coordinates[0])


def feature_for_row(
    *,
    row_index: int,
    cell_ids: Any,
    cell_index: Any,
    method: Any,
    num_vertices: Any,
    vertices: Any,
    polygon_set: str,
) -> tuple[dict[str, Any] | None, str | None, bool]:
    ring, stored_pairs = ring_from_flat_vertices(vertices[row_index])
    reported_vertices = int(num_vertices[row_index])
    truncated = reported_vertices > stored_pairs

    reason = validate_ring(ring)
    if reason is not None:
        return None, reason, truncated

    cell_id_0 = int(cell_ids[row_index, 0])
    cell_id_1 = int(cell_ids[row_index, 1])
    feature = {
        "type": "Feature",
        "properties": {
            "cell_id": f"{cell_id_0}_{cell_id_1}",
            "cell_id_0": cell_id_0,
            "cell_id_1": cell_id_1,
            "cell_index": int(cell_index[row_index]),
            "polygon_set": polygon_set,
            "method": int(method[row_index]),
            "reported_num_vertices": reported_vertices,
            "stored_coordinate_pairs": stored_pairs,
            "truncated_reported_vertices": truncated,
        },
        "geometry": make_geometry(ring),
    }

    return feature, validate_geojson_feature(feature), truncated


def main() -> int:
    args = parse_args()
    if not args.zarr_zip.exists():
        print(f"ERROR: input does not exist: {args.zarr_zip}", file=sys.stderr)
        return 2

    import zarr
    from zarr.storage import ZipStore

    features: list[dict[str, Any]] = []
    skipped: list[tuple[int, str]] = []
    truncated_rows = 0

    # store = zarr.ZipStore(str(args.zarr_zip), mode="r")
    store = ZipStore(str(args.zarr_zip), mode="r")
    try:
        root = zarr.open(store, mode="r")
        if "cell_id" not in root:
            raise KeyError("missing required array: cell_id")
        if "polygon_sets" not in root or args.polygon_set not in root["polygon_sets"]:
            raise KeyError(f"missing polygon set: polygon_sets/{args.polygon_set}")

        polygon_group = root["polygon_sets"][args.polygon_set]
        for required in ("cell_index", "method", "num_vertices", "vertices"):
            if required not in polygon_group:
                raise KeyError(
                    f"missing required array: polygon_sets/{args.polygon_set}/{required}"
                )

        cell_ids = root["cell_id"][:]
        cell_index = polygon_group["cell_index"][:]
        method = polygon_group["method"][:]
        num_vertices = polygon_group["num_vertices"][:]
        vertices = polygon_group["vertices"][:]

        row_count = vertices.shape[0]
        if not (
            cell_ids.shape[0]
            == cell_index.shape[0]
            == method.shape[0]
            == num_vertices.shape[0]
            == row_count
        ):
            raise ValueError(
                "array length mismatch: "
                f"cell_id={cell_ids.shape[0]}, "
                f"cell_index={cell_index.shape[0]}, "
                f"method={method.shape[0]}, "
                f"num_vertices={num_vertices.shape[0]}, "
                f"vertices={row_count}"
            )

        for row_index in range(row_count):
            feature, reason, truncated = feature_for_row(
                row_index=row_index,
                cell_ids=cell_ids,
                cell_index=cell_index,
                method=method,
                num_vertices=num_vertices,
                vertices=vertices,
                polygon_set=args.polygon_set,
            )

            if truncated:
                truncated_rows += 1
                if not args.allow_truncated:
                    skipped.append(
                        (
                            row_index,
                            "reported num_vertices is greater than stored coordinate pairs",
                        )
                    )
                    continue

            if reason is not None:
                skipped.append((row_index, reason))
                continue
            if feature is not None:
                features.append(feature)
    finally:
        store.close()

    if not features:
        print("ERROR: no valid polygons to write", file=sys.stderr)
        return 1

    if skipped and args.fail_on_skipped:
        for row_index, reason in skipped[:20]:
            print(f"SKIPPED row={row_index}: {reason}", file=sys.stderr)
        if len(skipped) > 20:
            print(f"... {len(skipped) - 20} more skipped rows", file=sys.stderr)
        print(
            "ERROR: skipped rows found and --fail-on-skipped was set", file=sys.stderr
        )
        return 1

    geojson = {"type": "FeatureCollection", "features": features}
    args.output_geojson.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = args.output_geojson.with_suffix(args.output_geojson.suffix + ".tmp")
    with tmp_path.open("w") as handle:
        json.dump(geojson, handle, separators=(",", ":"))
        handle.write("\n")
    tmp_path.replace(args.output_geojson)

    print(f"wrote: {args.output_geojson}")
    print(f"features written: {len(features)}")
    print(f"rows skipped: {len(skipped)}")
    print(
        f"rows with reported num_vertices > stored coordinate pairs: {truncated_rows}"
    )
    if skipped:
        print("first skipped rows:", file=sys.stderr)
        for row_index, reason in skipped[:10]:
            print(f"  row={row_index}: {reason}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
