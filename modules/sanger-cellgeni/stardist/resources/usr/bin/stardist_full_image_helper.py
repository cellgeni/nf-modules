#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2026 Wellcome Sanger Institute

import json
import logging
from pathlib import Path

import fire
import numpy as np
from csbdeep.utils import normalize
from shapely.geometry import Polygon, mapping
from stardist.models import StarDist2D


logging.basicConfig(level="INFO", format="[%(asctime)s][%(levelname)s] %(message)s")


def _select_plane(array: np.ndarray, channel: int, z: int) -> np.ndarray:
    image = np.asarray(array)
    image = np.squeeze(image)

    if image.ndim == 2:
        return image

    if image.ndim == 3:
        if image.shape[-1] <= 16:
            return image[..., channel]
        if image.shape[0] <= 16:
            return image[channel]
        return image[z]

    while image.ndim > 2:
        small_axes = [axis for axis, size in enumerate(image.shape) if size <= 16]
        axis = small_axes[0] if small_axes else 0
        index = channel if axis == small_axes[-1] else z
        image = np.take(image, index, axis=axis)
        image = np.squeeze(image)

    return image


def load_image_plane(
    image_path: str,
    channel: int = 0,
    z: int = 0,
    resolution_level: int = 0,
) -> np.ndarray:
    path = Path(image_path)

    if path.suffix.lower() == ".pgm":
        return _read_pgm(path)

    try:
        from aicsimageio import AICSImage

        image = AICSImage(path)
        return np.squeeze(
            image.get_image_data("YX", C=channel, Z=z, T=0, S=resolution_level)
        )
    except Exception as error:
        logging.info("Falling back to tifffile reader: %s", error)

    import tifffile

    with tifffile.TiffFile(path) as tif:
        series = tif.series[min(resolution_level, len(tif.series) - 1)]
        return _select_plane(series.asarray(), channel=channel, z=z)


def _read_pgm(path: Path) -> np.ndarray:
    tokens = []
    with open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            tokens.extend(line.split())

    if tokens[0] != "P2":
        raise ValueError(f"Unsupported PGM format in {path}")

    width = int(tokens[1])
    height = int(tokens[2])
    max_value = int(tokens[3])
    data = np.array(tokens[4:], dtype=np.float32).reshape(height, width)
    if max_value != 255:
        data = data / max_value * 255
    return data.astype(np.float32)


def segment(
    image_path: str,
    output_name: str,
    channel: int = 0,
    z: int = 0,
    resolution_level: int = 0,
    model_name: str = "2D_versatile_fluo",
):
    logging.info("Loading full image plane from '%s'", image_path)
    image = load_image_plane(
        image_path,
        channel=channel,
        z=z,
        resolution_level=resolution_level,
    )

    logging.info("Loading StarDist2D model '%s'", model_name)
    model = StarDist2D.from_pretrained(model_name)

    labels, details = model.predict_instances(normalize(image, 1, 99.8, axis=(0, 1)))
    coord = details["coord"]

    logging.info("Converting outlines to GeoJSON format")
    features = []
    if coord.shape[0] != 0:
        for object_id, polygon in enumerate(coord, start=1):
            flat_coords = [(xy[1], xy[0]) for xy in polygon.reshape(-1, 2)]
            features.append(
                {
                    "type": "Feature",
                    "properties": {"object_id": object_id},
                    "geometry": mapping(Polygon(flat_coords)),
                }
            )
    else:
        logging.info("No outlines found")

    geojson = {"type": "FeatureCollection", "features": features}
    with open(output_name, "wt") as handle:
        json.dump(geojson, handle, separators=(",", ":"))


if __name__ == "__main__":
    options = {
        "run": segment,
        "version": "0.0.1",
    }
    fire.Fire(options)
