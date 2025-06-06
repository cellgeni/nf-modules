#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import fire
from aicsimageio import AICSImage
import numpy as np
from ome_types import to_xml
from ome_types.model import TiffData
import pandas as pd
import os
import uuid


def generate_tiles_to_process_csv(
    img_obj: AICSImage,
):
    """
    Generate a CSV file with the tiles to process.
    """
    df = pd.DataFrame(
        {
            "image_id": [
                img_obj.metadata.images[i].id for i in np.arange(len(img_obj.scenes))
            ]
        }
    )
    df["index"] = np.arange(len(img_obj.scenes))
    return df


def append_tiffdata_to_ome_metadata(
    img_obj: AICSImage,
    index: int,
):
    """
    Append TiffData to the OME metadata.
    """
    current_image = img_obj.ome_metadata.images[index]
    filename = current_image.id + ".tif"
    tiff_uuid = TiffData.UUID(value=f"urn:uuid:{uuid.uuid4()}", file_name=filename)
    tiff = TiffData(first_c=0, first_t=0, first_z=0, uuid=tiff_uuid)
    img_obj.ome_metadata.images[index].pixels.tiff_data_blocks.append(tiff)


def main(
    image_root_folder, tiles_csv, companion_xml, out_folder, master_file="Index.idx.xml"
):
    """
    Generate a companion file for a given image file.
    """
    img = AICSImage(f"{image_root_folder}/{master_file}")
    df = generate_tiles_to_process_csv(img)
    df["root_xml"] = os.path.realpath(image_root_folder)
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
    df.to_csv(f"{out_folder}/{tiles_csv}", index=False)

    for i in np.arange(len(img.scenes)):
        append_tiffdata_to_ome_metadata(img, i)

    # Save the OME metadata as an XML file
    with open(f"{out_folder}/{companion_xml}", "w") as xml_file:
        xml_file.write(to_xml(img.ome_metadata))


def version():
    """
    Print the version of the script.
    """
    return "0.2.0"


if __name__ == "__main__":
    options = {
        "run": main,
        "version": version,
    }
    fire.Fire(options)
