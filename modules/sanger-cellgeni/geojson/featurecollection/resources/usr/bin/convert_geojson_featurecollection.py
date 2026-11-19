#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

VERSION = "1.0.0"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Normalize GeoJSON to minimal Polygon FeatureCollection format"
    )
    parser.add_argument("--input", required=False, help="Input GeoJSON file path")
    parser.add_argument("--output", required=False, help="Output GeoJSON file path")
    parser.add_argument(
        "--from-properties-id",
        action="store_true",
        help="Fallback to feature.properties.id if feature.id is missing",
    )
    parser.add_argument(
        "--keep-missing-id",
        action="store_true",
        help="Keep Polygon features without an id",
    )
    parser.add_argument("--version", action="store_true", help="Print version and exit")
    args = parser.parse_args()

    if args.version:
        return args
    if not args.input or not args.output:
        parser.error("--input and --output are required unless --version is used")
    return args


def normalize_feature(feature: dict, from_properties_id: bool, keep_missing_id: bool) -> dict | None:
    if not isinstance(feature, dict):
        return None
    geometry = feature.get("geometry")
    if not isinstance(geometry, dict) or geometry.get("type") != "Polygon":
        return None
    coordinates = geometry.get("coordinates")
    if coordinates is None:
        return None

    feature_id = feature.get("id")
    if feature_id is None and from_properties_id:
        properties = feature.get("properties")
        if isinstance(properties, dict):
            feature_id = properties.get("id")
    if feature_id is None and not keep_missing_id:
        return None

    out = {
        "type": "Feature",
        "geometry": {"type": "Polygon", "coordinates": coordinates},
    }
    if feature_id is not None:
        out["id"] = feature_id
    return out


def coerce_to_feature_list(payload: object) -> list[dict]:
    if isinstance(payload, dict):
        payload_type = payload.get("type")

        if payload_type == "FeatureCollection" and isinstance(payload.get("features"), list):
            return payload["features"]

        if payload_type == "GeometryCollection" and isinstance(payload.get("geometries"), list):
            return [
                {"type": "Feature", "id": str(idx), "geometry": geometry}
                for idx, geometry in enumerate(payload["geometries"])
            ]

        if isinstance(payload.get("geometries"), list):
            return [
                {"type": "Feature", "id": str(idx), "geometry": geometry}
                for idx, geometry in enumerate(payload["geometries"])
            ]

        if isinstance(payload.get("geometry"), dict):
            return [payload]

        if payload_type in {"Polygon", "MultiPolygon", "Point", "MultiPoint", "LineString", "MultiLineString"}:
            return [{"type": "Feature", "id": "0", "geometry": payload}]

    if isinstance(payload, list):
        return [
            {"type": "Feature", "id": str(idx), "geometry": geometry}
            for idx, geometry in enumerate(payload)
            if isinstance(geometry, dict)
        ]

    return []


def main() -> None:
    args = parse_args()
    if args.version:
        print(f"convert_geojson_featurecollection.py\t{VERSION}")
        return

    with Path(args.input).open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    features = coerce_to_feature_list(payload)

    normalized_features = []
    for feature in features:
        normalized = normalize_feature(feature, args.from_properties_id, args.keep_missing_id)
        if normalized is not None:
            normalized_features.append(normalized)

    output = {"type": "FeatureCollection", "features": normalized_features}
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(output, handle, ensure_ascii=False, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
