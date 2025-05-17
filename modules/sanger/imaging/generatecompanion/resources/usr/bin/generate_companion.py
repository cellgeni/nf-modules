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
    md = img.metadata

    wells = []
    if len(md.plates) > 0:
        # Plate format
        for well in md.plates[0].wells:
            print(well)
            for fov in well.well_samples:
                row = fov.model_dump()
                fov_md = {}
                fov_md["index"] = fov.index
                fov_md['well_id'] = well.id
                fov_md['fov_id'] = row['image_ref']['id']
                wells.append(fov_md)
    else:
        raise ValueError("No plates found in the metadata.")
    df = pd.DataFrame(wells)
    df["root_folder"] = os.path.realpath(file_with_ome_md)

    df.to_csv(tiles_csv, index=False)

    # Save the OME metadata as an XML file
    ome_str = to_xml(md)
    with open(companion_xml, "w") as xml_file:
        xml_file.write(ome_str)

def version():
    """
    Print the version of the script.
    """
    return "0.2.0"

if __name__ == "__main__":
    options = {
        "run" : main,
        "version" : version,
    }
    fire.Fire(options)