#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import fire
from aicsimageio import AICSImage
from aicsimageio.writers.ome_tiff_writer import OmeTiffWriter
import numpy as np
import tifffile as tf
import time  # Import the time module
from clij2fft.richardson_lucy import richardson_lucy_nc, richardson_lucy
# from cucim.skimage.restoration import richardson_lucy
# from skimage.restoration import richardson_lucy
from pathlib import Path
import pyopencl as cl
import logging
import os

# Set up logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
)
logger = logging.getLogger(__name__)


def load_and_process_psf(
        z_stack,
        psf_file,
        original_z_step,
        psf_z_step=0.1
    ):
    """
    Load the PSF file and process it to match the dimensions of the input stack.
    """
    psf = tf.imread(psf_file)
    psf_shape = psf.shape
    if len(psf_shape) != 3:
        raise ValueError("PSF file must have 3 dimensions.")
    Z = z_stack.shape[0]
    # Get the shape of the PSF
    if Z > psf_shape[0]:
        raise ValueError("Z planes in the input stack is greater than the PSF Z planes.")
    step = original_z_step // psf_z_step
    print(step)
    indices_to_keep = np.array([psf_shape[0]//2 + int(step) * (i - Z // 2) for i in range(Z)])
    indices_to_keep = indices_to_keep[(indices_to_keep >= 0) & (indices_to_keep < psf_shape[0])]
    print(indices_to_keep)
    # Subsample the PSF to match the number of Z planes in the input stack
    return psf[indices_to_keep, :, :]


def main(root_folder, index, out_img_name, well,
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
    physical_pixel_sizes = img.physical_pixel_sizes
    print(f"Physical pixel sizes: {physical_pixel_sizes}")
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
            current_psf = f"{psf_folder}/{c_name}.tif"
            if Path(current_psf).exists():
                logger.info(f"PSF file found: {current_psf}")
                psf = load_and_process_psf(cz_stack[c], current_psf, physical_pixel_sizes.Z)
                before_decon = time.time()
                z_stack = richardson_lucy_nc(cz_stack[c], psf, iterations)
                # z_stack = richardson_lucy(cz_stack[c], psf, iterations)
                # z_stack = richardson_lucy(cp.asarray(cz_stack[c]), cp.asarray(psf), iterations, clip=False).get()
                # z_stack = richardson_lucy(cz_stack[c], psf, iterations, filter_epsilon=0.05)
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
    if not os.path.exists(f"{well}"):
        os.makedirs(f"{well}")
    OmeTiffWriter.save(
        processed_hyper_stack, f"{well}/{out_img_name}",
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