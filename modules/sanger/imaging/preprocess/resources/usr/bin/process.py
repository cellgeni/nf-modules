#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import fire
from aicsimageio import AICSImage
from aicsimageio.writers.ome_tiff_writer import OmeTiffWriter
import numpy as np
import tifffile as tf
import time  # Import the time module
from clij2fft.richardson_lucy import richardson_lucy_nc
from pathlib import Path
import pyopencl as cl
import logging

# Set up logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
)

logger = logging.getLogger(__name__)


def main(root_folder, index, out_img_name,
         iterations=100,
         master_file="Index.idx.xml",
         psf_folder="psfs",
         z_project=True):
    """
    Generate a companion file for a given image file.
    """
    try:
        platforms = cl.get_platforms()
        if len(platforms) <= 0:
            raise RuntimeError("Could not find a valid open cl platform. Check your enviroment.")
        devices=platforms[0].get_devices()

        for device in devices:
            logger.info(f"Found open CL device: {device}")
            logger.info(f"Device has {device.get_info(cl.device_info.GLOBAL_MEM_SIZE)} mem available.")
    except:
        logger.warning("Could not find a valid open cl platform. Fall back to CPU.")

    # Start the timer
    start_time = time.time()
    cursor = start_time
    img = AICSImage(f"{root_folder}/{master_file}")
    # Print the time taken to load the image
    print(f"Time taken to load the image: {time.time() - cursor} seconds")
    cursor = time.time()
    img.set_scene(img.scenes[index])
    print(f"Elapsed time for scene setting: {time.time() - cursor} seconds")
    cursor = time.time()
    processed_hyper_stack = []
    for t in range(img.dims.T):
        c_stack = []
        cursor = time.time()
        cz_stack = img.get_image_data("CZYX", T=t)
        for c, c_name in enumerate(img.channel_names):
            print(c_name)
            current_psf = f"{root_folder}/{psf_folder}/{c_name}.tif"
            if Path(current_psf).exists():
                logger.info(f"PSF file found: {current_psf}")
                psf = tf.imread(current_psf)
                before_decon = time.time()
                z_stack = richardson_lucy_nc(cz_stack[c], psf, iterations)
                after_decon = time.time()
                logger.info(f"Deconvolution time: {after_decon - before_decon} seconds")
            else:
                logger.info(f"PSF file not found: {current_psf}. Using original image data.")
                z_stack = cz_stack[c]
            if z_project:
                z_stack = np.max(z_stack, axis=0)
            c_stack.append(z_stack)
        processed_hyper_stack.append(c_stack)
    processed_hyper_stack = np.array(processed_hyper_stack)
    print(processed_hyper_stack.shape)
    cursor = time.time()
    new_dim_order = "TCYX" if z_project else "TCZYX"
    OmeTiffWriter.save(
        processed_hyper_stack, out_img_name,
        dim_order=new_dim_order,
        channel_names=img.channel_names,
        image_names=out_img_name.replace(".ome.tif", ""),
        physical_pixel_sizes=img.physical_pixel_sizes,
    )
    print(f"Elapsed time for saving the image: {time.time() - cursor} seconds")
    

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