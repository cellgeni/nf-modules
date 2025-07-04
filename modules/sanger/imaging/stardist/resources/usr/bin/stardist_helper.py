#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2025 Wellcome Sanger Institute

import numpy as np
from csbdeep.utils import normalize
from stardist.models import StarDist2D
from csbdeep.data import Normalizer, normalize_mi_ma
import fire
from imagetileprocessor import slice_and_crop_image
from shapely import wkt, MultiPolygon

import logging


class MyNormalizer(Normalizer):
    def __init__(self, mi, ma):
        self.mi, self.ma = mi, ma

    def before(self, x, axes):
        return normalize_mi_ma(x, self.mi, self.ma, dtype=np.float32)

    def after(*args, **kwargs):
        assert False

    @property
    def do_after(self):
        return False


# configure logging
logging.basicConfig(level="INFO", format="[%(asctime)s][%(levelname)s] %(message)s")


def segment(
    image_path: str,
    x_min: int,
    x_max: int,
    y_min: int,
    y_max: int,
    resolution_level: int = 0,
    model_name: str = "2D_versatile_fluo",
    output_name: str = None,
    DAPI_index: int = 0,
    z: int = 0,
):
    """
    Main function to perform image segmentation using tiles.

    Parameters:
        image_path (str): The path to the image.
        x_min (int): The minimum x-coordinate of the tile.
        x_max (int): The maximum x-coordinate of the tile.
        y_min (int): The minimum y-coordinate of the tile.
        y_max (int): The maximum y-coordinate of the tile.
        resolution_level (str): The pyramid level to use for segmentation.
        model_name (str): The StarDist model type to use for segmentation.
        output (str): The output directory to save the segmentation results.
        image_id (str): The image identification name.
        DAIP_index (int): The index of the DAPI channel.
    """
    img = slice_and_crop_image(
        image_path,
        x_min,
        x_max,
        y_min,
        y_max,
        channel=np.array([DAPI_index]),
        zs=np.array([z]),
        resolution_level=resolution_level,
    )

    logging.info(f"Loading StarDist2D model '{model_name}'")
    model = StarDist2D.from_pretrained(model_name)

    logging.info(f"Loading full image")

    # mi, ma = np.percentile(crop[::8], [1,99.8])                      # compute percentiles from low-resolution image
    # normalizer = MyNormalizer(mi, ma)
    labels, details = model.predict_instances(normalize(img, 1, 99.8, axis=(0, 1)))
    coord, points, prob = details["coord"], details["points"], details["prob"]

    # logging.info(f"Normalize image")
    # norm_image = normalize(image, 1, 99.8, axis=(0, 1))

    # convert cellpose outlines to WTK
    logging.info(f"Converting outlines to WKT format")
    wkts = []
    if coord.shape[0] != 0:
        for polygon in coord:
            flat_coords = [
                (xy[1] + x_min, xy[0] + y_min) for xy in polygon.reshape(-1, 2)
            ]
            wkts.append(
                wkt.loads(
                    "POLYGON (("
                    + ", ".join(f"{x} {y}" for x, y in flat_coords + [flat_coords[0]])
                    + "))"
                )
            )
        with open(output_name, "wt") as f:
            f.write(wkt.dumps(MultiPolygon(wkts)))
    else:
        logging.info("No outlines file found")
        with open(output_name, "wt") as f:
            f.write("")


if __name__ == "__main__":
    options = {"run": segment, "version": "0.0.2"}
    fire.Fire(options)
