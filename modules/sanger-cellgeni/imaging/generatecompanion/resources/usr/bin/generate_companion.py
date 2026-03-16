#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import fire
from bioio import BioImage
import numpy as np
from ome_types import from_xml, to_xml
from ome_types.model import TiffData
import pandas as pd
import os
import uuid
import re


def generate_image_paths_from_ome(ome):
    """
    Build {image_id_int: "A/1/0"} mapping from an OME object.
    Uses image.name (e.g. "Well A01, Field 1") when available,
    and falls back to plate/well_samples metadata otherwise.
    """
    out = {}

    # 1) Try to parse directly from image.name
    pat = re.compile(r"Well\s+([A-Za-z]+)(\d+)\s*,\s*Field\s+(\d+)")
    for img in ome.images:
        image_id = int(str(img.id).split(":")[-1])
        name = img.name or ""
        m = pat.search(name)
        if m:
            row_label, col_str, field_str = m.groups()
            col_num = int(col_str)
            fov_idx = int(field_str) - 1
            out[image_id] = f"{row_label.upper()}/{col_num}/{fov_idx}"

    # 2) Fallback from plate metadata for any image not mapped above
    if len(out) < len(ome.images) and ome.plates:
        for plate in ome.plates:
            for well in plate.wells:
                row_label = chr(ord("A") + well.row)
                col_num = well.column + 1
                samples = sorted(well.well_samples, key=lambda s: s.index)
                for local_fov_idx, sample in enumerate(samples):
                    image_id = int(str(sample.image_ref.id).split(":")[-1])
                    out.setdefault(image_id, f"{row_label}/{col_num}/{local_fov_idx}")

    return dict(sorted(out.items()))


def generate_tiles_to_process_csv(img_obj: BioImage, prefix: str):
    """
    Generate a CSV file with the tiles to process.
    """
    df = pd.DataFrame(
        {
            "image_id": [
                prefix + "_" + img_obj.metadata.images[i].id
                for i in np.arange(len(img_obj.scenes))
            ]
        }
    )
    df["index"] = np.arange(len(img_obj.scenes))
    return df


def generate_tiles_to_process_csv_from_ome(ome, prefix: str):
    """
    Generate a CSV file with the tiles to process from an OME object.
    """
    df = pd.DataFrame({"image_id": [prefix + "_" + image.id for image in ome.images]})
    df["index"] = np.arange(len(ome.images))
    return df


def append_tiffdata_to_ome_metadata(
    img_obj: BioImage,
    index: int,
    prefix: str,
):
    """
    Append TiffData to the OME metadata.
    """
    current_image = img_obj.ome_metadata.images[index]
    filename = prefix + "_" + current_image.id + ".tif"
    tiff_uuid = TiffData.UUID(value=f"urn:uuid:{uuid.uuid4()}", file_name=filename)
    tiff = TiffData(first_c=0, first_t=0, first_z=0, uuid=tiff_uuid)
    img_obj.ome_metadata.images[index].pixels.tiff_data_blocks.append(tiff)


def append_tiffdata_to_ome_image(current_image, prefix: str):
    """
    Append a TiffData block to a single OME image object.
    """
    filename = prefix + "_" + current_image.id + ".tif"
    tiff_uuid = TiffData.UUID(value=f"urn:uuid:{uuid.uuid4()}", file_name=filename)
    tiff = TiffData(first_c=0, first_t=0, first_z=0, uuid=tiff_uuid)
    current_image.pixels.tiff_data_blocks.append(tiff)


def main(
    image_root_folder, tiles_csv, companion_xml, prefix, master_file="Index.idx.xml"
):
    """
    Generate a companion file for a given image file.
    """
    ome_xml_relative_path = os.path.join("OME", "METADATA.ome.xml")
    ome_xml_path = os.path.join(image_root_folder, ome_xml_relative_path)
    is_omezarr_plate = os.path.isfile(ome_xml_path)

    if is_omezarr_plate:
        ome = from_xml(ome_xml_path)
        df = generate_tiles_to_process_csv_from_ome(ome, prefix)
        image_paths = generate_image_paths_from_ome(ome)
        df["path"] = [
            image_paths.get(int(str(image.id).split(":")[-1]), "")
            for image in ome.images
        ]
        effective_master_file = ome_xml_relative_path
    else:
        img = BioImage(f"{image_root_folder}/{master_file}")
        ome = img.ome_metadata
        df = generate_tiles_to_process_csv(img, prefix)
        image_paths = generate_image_paths_from_ome(ome)
        df["path"] = [
            image_paths.get(
                int(str(img.metadata.images[i].id).split(":")[-1]),
                "",
            )
            for i in np.arange(len(img.scenes))
        ]
        effective_master_file = master_file

    df["root_xml"] = os.path.realpath(image_root_folder)
    df["master_file"] = effective_master_file
    df.to_csv(f"{tiles_csv}", index=False)

    for image in ome.images:
        append_tiffdata_to_ome_image(image, prefix)

    # Save the OME metadata as an XML file
    with open(f"{companion_xml}", "w") as xml_file:
        xml_file.write(to_xml(ome))


def version():
    """
    Print the version of the script.
    """
    return "0.2.1"


if __name__ == "__main__":
    options = {
        "run": main,
        "version": version,
    }
    fire.Fire(options)
