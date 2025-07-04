#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (c) 2024 Wellcome Sanger Institute

"""
This script will slice the image in XY dimension and save the slices coordinates in json files
"""
import fire
from scipy import ndimage
from instanseg import InstanSeg
import cv2
import numpy as np
from shapely.geometry import Polygon, MultiPolygon
from shapely import to_wkt
from shapely.affinity import translate
from imagetileprocessor import slice_and_crop_image

import logging

logger = logging.getLogger(__name__)


def get_largest_polygon(multi_polygon: MultiPolygon):
    # Initialize variables to store the largest polygon and its area
    largest_polygon = None
    largest_area = 0

    # Iterate through each polygon in the MultiPolygon
    for polygon in multi_polygon.geoms:
        # Calculate the area of the current polygon
        area = polygon.area
        # Update the largest polygon if the current one is larger
        if area > largest_area:
            largest_area = area
        largest_polygon = polygon
    if largest_polygon is None:
        logger.warning("No polygon found")
        return multi_polygon
    else:
        return largest_polygon


def get_shapely(label):
    """
    get outlines of masks as a list to loop over for plotting
    """
    polygons = {}
    simpler_polys = {}
    slices = ndimage.find_objects(label)
    for i, bbox in enumerate(slices):
        if not bbox:
            continue
        cur_cell_label = i + 1
        msk = (label[bbox[0], bbox[1]] == cur_cell_label).astype(np.uint8).copy()
        cnts, _ = cv2.findContours(msk, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
        # if len(cnts) > 1:
        # print(len(cnts), cur_cell_label)
        current_polygons = [
            Polygon((cnt + [bbox[1].start, bbox[0].start]).squeeze())
            #             else Point((cnt + [bbox[1].start, bbox[0].start]).squeeze())
            for cnt in cnts
            if len(cnt) > 2
        ]
        multipoly_obj = MultiPolygon(current_polygons)
        polygons[cur_cell_label] = multipoly_obj
        simpler_polys[cur_cell_label] = get_largest_polygon(multipoly_obj).simplify(
            1, preserve_topology=True
        )
    return polygons, simpler_polys


def main(
    image_path: str,
    x_min: int,
    x_max: int,
    y_min: int,
    y_max: int,
    output_name: str,
    model: str = "fluorescence_nuclei_and_cells",
    C: list = [0],
    Z: list = [0],
):
    instanseg_fluorescence = InstanSeg(model, verbosity=1)

    crop = slice_and_crop_image(
        image_path,
        x_min,
        x_max,
        y_min,
        y_max,
        channel=np.array([C]),
        zs=np.array([Z]),
        resolution_level=0,
    )

    labeled_output = instanseg_fluorescence.eval_small_image(
        np.array(crop).astype(np.uint16),
        None,
        return_image_tensor=False,
        target="nuclei",
        resolve_cell_and_nucleus=False,
        cleanup_fragments=True,
    )
    polys = get_shapely(np.squeeze(np.array(labeled_output)).astype(np.uint16))
    with open(output_name, "wt") as f:
        f.write(
            to_wkt(
                translate(MultiPolygon(list(polys[1].values())), xoff=x_min, yoff=y_min)
            )
        )


if __name__ == "__main__":
    options = {
        "run": main,
        "version": "0.0.1",
    }
    fire.Fire(options)
