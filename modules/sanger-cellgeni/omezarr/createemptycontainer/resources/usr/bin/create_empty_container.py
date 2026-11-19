#!/usr/bin/env python3

import json
import logging

import fire
from ngio import (
    create_empty_ome_zarr,
    create_empty_plate,
    open_ome_zarr_container,
    open_ome_zarr_plate,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def _load_tiles(tiles_json):
    """
    fire pre-parses CLI args that look like Python literals, so tiles_json may
    already be a list by the time it reaches here.
    """
    if isinstance(tiles_json, str):
        return json.loads(tiles_json)
    return list(tiles_json)


def _expected_shape_and_axes(img, z_project):
    """
    Mirror process.py's own read/projection logic so the placeholder created here
    has exactly the shape process.py will later write into.
    """
    stack = img.get_as_dask(axes_order=("t", "c", "z", "y", "x"))
    t, c, z, y, x = stack.shape
    if z_project:
        return (t, c, y, x), ("t", "c", "y", "x")
    return (t, c, z, y, x), ("t", "c", "z", "y", "x")


def main(root_folder, tiles_json, out_zarr, z_project=True, overwrite=True):
    """
    Pre-create one empty, correctly-shaped image slot per tile in a fresh output
    OME-Zarr store (an HCS plate when any tile carries an hcs_path, otherwise a
    standalone whole-slide image), so that OMEZARR_PREPROCESS tasks can later
    write into it concurrently without racing on shared zarr group metadata.
    """
    tiles = _load_tiles(tiles_json)
    non_hcs_tiles = [t for t in tiles if not t[1]]
    if len(non_hcs_tiles) > 1:
        raise ValueError(
            "Non-HCS OME-Zarr targets support exactly one tile per out_zarr; "
            f"got {len(non_hcs_tiles)} tile(s) without an hcs_path."
        )

    plate_in = None
    plate_out = None

    for image_id, hcs_path in tiles:
        if hcs_path:
            if plate_in is None:
                plate_in = open_ome_zarr_plate(root_folder)
            if plate_out is None:
                plate_out = create_empty_plate(
                    store=out_zarr, name="preprocessed", overwrite=overwrite
                )
            row, column, fov = hcs_path.split("/")
            container_in = plate_in.get_image(row=row, column=column, image_path=fov)
            img_in = container_in.get_image()
            shape, axes_names = _expected_shape_and_axes(img_in, z_project)

            plate_out.add_image(row, column, fov)
            target_store = plate_out.get_image_store(row, column, fov)
        else:
            container_in = open_ome_zarr_container(root_folder)
            img_in = container_in.get_image()
            shape, axes_names = _expected_shape_and_axes(img_in, z_project)
            target_store = out_zarr

        pixelsize = img_in.pixel_size
        create_empty_ome_zarr(
            store=target_store,
            shape=shape,
            pixelsize=(pixelsize.x, pixelsize.y),
            z_spacing=1.0 if z_project else pixelsize.z,
            axes_names=axes_names,
            channels_meta=list(container_in.channel_labels),
            dtype="uint16",
            name=image_id,
            overwrite=overwrite,
        )
        logger.info(
            f"Pre-created empty image {image_id!r} (hcs_path={hcs_path!r}) "
            f"with shape {shape}."
        )


def version():
    """
    Print the version of the script.
    """
    return "0.1.0"


if __name__ == "__main__":
    options = {
        "run": main,
        "version": version,
    }
    fire.Fire(options)
