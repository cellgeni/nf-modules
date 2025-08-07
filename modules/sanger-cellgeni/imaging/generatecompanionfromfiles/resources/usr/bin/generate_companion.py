#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import fire
from ome_types.model import OME
from ome_types.model import TiffData
from ome_types.model import Plate, Well, WellSample, ImageRef
from aicsimageio import AICSImage
from ome_types import to_xml
import uuid
from glob import glob

import os
import re


def get_alphabet_index(alphabet_str):
    """
    Get the index of a character in the alphabet string.
    """
    return "ABCDEFGHIJKLMNOPQRSTUVWXYZ".index(alphabet_str.upper())


def main(input_folder, regex, out_companion_xml, out_folder="./", channel_names=None):
    ome = OME(uuid=f"urn:uuid:{uuid.uuid4()}")
    images = glob(f"{input_folder}/{regex}")
    print(
        f"Found {len(images)} images matching the regex '{regex}' in folder '{input_folder}'."
    )
    plate_ome = Plate(
        id=f"Plate_{uuid.uuid4()}",
        name="Plate 1",
        well_count=len(images),
    )
    ome.plates.append(plate_ome)

    for ind, img in enumerate(images):
        img_basename = os.path.basename(img)
        row = re.search(r"([A-Z])", os.path.basename(img)).group(
            1
        )  # Extract row name using regex
        column = re.search(r"(\d{2})", os.path.basename(img)).group(
            1
        )  # Extract column name using regex
        img_obj = AICSImage(img)

        image_id = f"Image:{ind}"
        new_img_obj = img_obj.ome_metadata.images[0]
        if channel_names:
            assert len(new_img_obj.pixels.channels) == len(channel_names)
        else:
            channel_names = [
                f"Channel {i + 1}" for i in range(len(new_img_obj.pixels.channels))
            ]
        for i, channel in enumerate(new_img_obj.pixels.channels):
            # print(f"Channel {i}: {channel.name} ({channel.id})")
            new_img_obj.pixels.channels[i].name = channel_names[i]
        new_img_obj.id = image_id
        new_img_obj.name = img_basename
        # print("Rows:", row, get_alphabet_index(row))
        # print("Columns:", column, int(column))
        well = Well(
            id=f"Well:{ind}",
            row=get_alphabet_index(row),
            column=int(column),
        )
        well_sample = WellSample(
            index=0,
            image_ref=ImageRef(id=image_id),
        )
        well.well_samples.append(well_sample)
        new_img_obj.pixels.tiff_data_blocks[0] = TiffData(
            first_c=(
                new_img_obj.pixels.tiff_data_blocks[0].first_c
                if new_img_obj.pixels.tiff_data_blocks
                else 0
            ),
            first_t=(
                new_img_obj.pixels.tiff_data_blocks[0].first_t
                if new_img_obj.pixels.tiff_data_blocks
                else 0
            ),
            first_z=(
                new_img_obj.pixels.tiff_data_blocks[0].first_z
                if new_img_obj.pixels.tiff_data_blocks
                else 0
            ),
            uuid=TiffData.UUID(
                value=f"urn:uuid:{uuid.uuid4()}", file_name=img_basename
            ),
        )
        ome.plates[0].wells.append(well)
        ome.images.append(new_img_obj)
    with open(f"{out_folder}/{out_companion_xml}", "w") as xml_file:
        xml_file.write(to_xml(ome))
    return ome


if __name__ == "__main__":
    options = {"version": "0.1.0", "run": main}
    fire.Fire(options)
