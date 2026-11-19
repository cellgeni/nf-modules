#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import argparse
# from pathlib import Path
# import scanpy as sc
import anndata as ad
import numpy as np
import sys

# NEMO
from terra import harmonize_tokenize_embed_pipeline

__version__ = "0.1.1"

def add_terra_features(args):

    # === 1. Load h5ad ===
    h5ad_path = args.h5ad_path
    assert os.path.exists(h5ad_path), f"h5ad file {h5ad_path} does not exist"
    if h5ad_path.endswith('.h5ad'):
        adata = ad.read_h5ad(h5ad_path)
    elif h5ad_path.endswith('.zarr'):
        adata = ad.read_zarr(h5ad_path)
    else:
        raise ValueError(f"h5ad file {h5ad_path} is not in h5ad or zarr format")

    assert 'spatial' in adata.obsm.keys(), f"h5ad file {h5ad_path} does not have spatial coordinates in obsm['spatial']"


    # === 2. Add NEMO embedding ===
    model_folder_path = '/opt/terra_model'
    assert os.path.exists(model_folder_path), "terra model folder path does not exists"

    adata = harmonize_tokenize_embed_pipeline(
        adata=adata,
        sample_key=None,          # a single sample
        batch_key="sample",
        model_folder_path=model_folder_path,
        cache_directory_path=args.tokenized_data_path,
        batch_size=128,           # lower this if you run out of GPU memory
    )


    # Save embeddings as CSV with high precision.
    output_dir = os.path.dirname(args.output_file) or "."
    output_stem = os.path.splitext(os.path.basename(args.output_file))[0]
    cell_csv_path = os.path.join(output_dir, f"{output_stem}_cell_emb.csv")
    neighborhood_csv_path = os.path.join(output_dir, f"{output_stem}_neighborhood_emb.csv")
    np.savetxt(
        cell_csv_path,
        np.asarray(adata.obsm['cell_emb'], dtype=np.float64),
        delimiter=",",
        fmt="%.18e",
    )
    np.savetxt(
        neighborhood_csv_path,
        np.asarray(adata.obsm['neighborhood_emb'], dtype=np.float64),
        delimiter=",",
        fmt="%.18e",
    )

    # === 3. Save anndata ===
    adata.write(args.output_file)
    


if __name__ == "__main__":
    if len(sys.argv) == 2 and sys.argv[1] == "version":
        print(__version__)
        raise SystemExit(0)

    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad_path", type=str)
    parser.add_argument("--tokenized_data_path", type=str)
    parser.add_argument("--output_file", type=str)
    parser.add_argument("--version", action='store_true', help="Print version and exit")
    args = parser.parse_args()
    if args.version:
        print(__version__)
        raise SystemExit(0)
    add_terra_features(args)
