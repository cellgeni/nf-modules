#!/usr/bin/env python3

import json
import logging

import fire
from ngio import open_ome_zarr_container, open_ome_zarr_plate

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


def main(out_zarr, tiles_json):
    """
    Verify that every tile expected in this OME-Zarr store actually has data
    written, once all parallel OMEZARR_PREPROCESS tasks for the dataset have
    completed. Meant as a final gate before the store is considered done.
    """
    tiles = _load_tiles(tiles_json)

    plate = None
    missing = []
    for image_id, hcs_path in tiles:
        if hcs_path:
            if plate is None:
                plate = open_ome_zarr_plate(out_zarr)
            row, column, fov = hcs_path.split("/")
            img = plate.get_image(row=row, column=column, image_path=fov).get_image()
        else:
            img = open_ome_zarr_container(out_zarr).get_image()

        n_written = img.zarr_array.nchunks_initialized
        if n_written == 0:
            missing.append(f"{image_id} ({hcs_path or 'n/a'})")
        else:
            logger.info(
                f"Validated {image_id!r} (hcs_path={hcs_path!r}): "
                f"{n_written} chunk(s) written."
            )

    if missing:
        raise RuntimeError(
            f"{len(missing)}/{len(tiles)} expected tile(s) have no data written "
            f"in {out_zarr}: {', '.join(missing)}"
        )
    logger.info(f"All {len(tiles)} tile(s) validated in {out_zarr}.")


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
