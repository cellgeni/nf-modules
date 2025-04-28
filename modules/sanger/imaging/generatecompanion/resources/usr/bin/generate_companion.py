#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import fire
from aicsimageio import AICSImage
import numpy as np
from ome_types import to_xml
import pandas as pd
import os


def main(file_with_ome_md, tiles_csv, companion_xml, master_file="Index.idx.xml"):
    """
    Generate a companion file for a given image file.
    """
    img = AICSImage(f"{file_with_ome_md}/{master_file}")
    absolute_path = os.path.realpath(file_with_ome_md)
    ome_str = to_xml(img.metadata)
    df = pd.DataFrame(
        {
            "image_id" : [img.metadata.images[i].id for i in np.arange(len(img.scenes))]
        }
    )
    df["root_xml"] = absolute_path
    df["index"] = np.arange(len(img.scenes))

    df.to_csv(tiles_csv, index=False)

    # Save the OME metadata as an XML file
    with open(companion_xml, "w") as xml_file:
        xml_file.write(ome_str)

def version():
    """
    Print the version of the script.
    """
    return "0.1.0"

if __name__ == "__main__":
    options = {
        "run" : main,
        "version" : version,
    }
    fire.Fire(options)