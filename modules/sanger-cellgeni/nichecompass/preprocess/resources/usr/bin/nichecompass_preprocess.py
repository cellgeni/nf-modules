#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nichecompass_preprocess.py

Preprocess spatial H5AD batches for NicheCompass training:
  1. Copy NicheCompass reference data from the container image
  2. Load H5AD batches, falling back to .X if the counts layer is absent
  3. Compute spatial neighbour graphs with Squidpy
  4. Concatenate batches with a block-diagonal adjacency matrix
  5. Sanitize data (NaN / Inf handling)
  6. Build combined prior GP mask (OmniPath + NicheNet v2 + MEBOCOST)
  7. Embed GP masks into AnnData

Outputs written to the working directory:
  <prefix>_preprocessed.h5ad  — preprocessed AnnData ready for training
  data/                       — NicheCompass reference data for downstream analysis
"""

import argparse
import json
import logging
import random
import shutil
import sys
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any, cast, Literal, TypeAlias

import anndata as ad
import numpy as np
import scanpy as sc
import squidpy as sq
import scipy.sparse as sp
import torch

from nichecompass.models import NicheCompass
from nichecompass.utils import (
    add_gps_from_gp_dict_to_adata,
    filter_and_combine_gp_dict_gps_v2,
    extract_gp_dict_from_mebocost_ms_interactions,
    extract_gp_dict_from_nichenet_lrt_interactions,
    extract_gp_dict_from_omnipath_lr_interactions,
)

Species: TypeAlias = Literal["human", "mouse"]

_ADJ_KEY: str = "spatial_connectivities"


def setup_logging(log_path: Path, debug: bool) -> None:
    level = logging.DEBUG if debug else logging.INFO
    fmt = "%(asctime)s | %(levelname)s | %(message)s"
    datefmt = "%Y%m%d %H:%M:%S"
    logging.basicConfig(
        level=level,
        format=fmt,
        datefmt=datefmt,
        handlers=[
            logging.FileHandler(log_path, encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
        force=True,
    )


def fixed_seeds(seed: int = 0) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if getattr(torch, "cuda", None) and torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    import os
    os.environ["PYTHONHASHSEED"] = str(seed)


#### Parameter dataclass and parsing ####

@dataclass
class PreprocessParams:
    # MAIN
    batches: list[Path]
    prefix: str
    species: Species = "human"
    debug: bool = False

    # DATASET / GRAPH
    spatial_key: str = "spatial"
    n_neighbors: int = 4
    sample_key: str = "batch"

    # AnnData keys
    counts_key: str = "counts"
    gp_names_key: str = "nichecompass_gp_names"
    active_gp_names_key: str = "nichecompass_active_gp_names"
    gp_targets_mask_key: str = "nichecompass_gp_targets"
    gp_targets_categories_mask_key: str = "nichecompass_gp_targets_categories"
    gp_sources_mask_key: str = "nichecompass_gp_sources"
    gp_sources_categories_mask_key: str = "nichecompass_gp_sources_categories"
    latent_key: str = "nichecompass_latent"

    # Metadata for downstream use (training, analysis)
    cell_type_key: str = "Main_molecular_cell_type"
    cat_covariates_keys: list[str] | None = None


def load_config_json(path: Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        cfg = json.load(f)
    if not isinstance(cfg, dict):
        raise ValueError("--config JSON must be an object/dict at top level")
    return cfg


def validate_known_keys(config: dict[str, Any], allowed: list[str]) -> None:
    unknown = sorted(set(config.keys()) - set(allowed))
    if unknown:
        raise ValueError(f"Unknown keys in --config: {unknown}")


def _validate_and_cast_species(val: str) -> Species:
    if val not in ("human", "mouse"):
        raise ValueError(f"Invalid species {val!r}; must be 'human' or 'mouse'.")
    return cast(Species, val)


def str2bool(v: str) -> bool:
    if isinstance(v, bool):
        return v
    val = v.strip().lower()
    if val in {"true", "t", "1", "yes", "y"}:
        return True
    if val in {"false", "f", "0", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError(f"Invalid bool: {v!r}")


def build_parser() -> tuple[argparse.ArgumentParser, argparse.ArgumentParser]:
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config", type=Path, default=None, help="Path to JSON config file.")

    parser = argparse.ArgumentParser(
        prog="nichecompass_preprocess.py",
        description="Preprocess spatial H5AD batches for NicheCompass training.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[pre],
    )

    g_main = parser.add_argument_group("MAIN")
    g_main.add_argument("--batches", nargs="+", type=Path, default=argparse.SUPPRESS,
                        help="Paths to input .h5ad files (>=1).")
    g_main.add_argument("--prefix", type=str, required=True,
                        help="Output file prefix (e.g. sample1 → sample1_preprocessed.h5ad).")
    g_main.add_argument("--species", type=str, choices=["human", "mouse"], default="human",
                        help="Species tag for prior knowledge lookup.")
    g_main.add_argument("--debug", action="store_true", help="Enable DEBUG logging.")

    g_dataset = parser.add_argument_group("DATASET / GRAPH")
    g_dataset.add_argument("--sample_key", type=str, default="batch",
                           help="obs key for sample/batch.")
    g_dataset.add_argument("--cell_type_key", type=str, default="Main_molecular_cell_type",
                           help="obs key for cell type labels.")
    g_dataset.add_argument("--spatial_key", type=str, default="spatial",
                           help="obsm key for spatial coordinates.")
    g_dataset.add_argument("--n_neighbors", type=int, default=4,
                           help="Number of spatial neighbours per node.")

    g_ad = parser.add_argument_group("AnnData keys")
    g_ad.add_argument("--counts_key", type=str, default="counts",
                      help="Layer name for counts (falls back to X if absent).")
    g_ad.add_argument("--gp_names_key", type=str, default="nichecompass_gp_names")
    g_ad.add_argument("--active_gp_names_key", type=str, default="nichecompass_active_gp_names")
    g_ad.add_argument("--gp_targets_mask_key", type=str, default="nichecompass_gp_targets")
    g_ad.add_argument("--gp_targets_categories_mask_key", type=str,
                      default="nichecompass_gp_targets_categories")
    g_ad.add_argument("--gp_sources_mask_key", type=str, default="nichecompass_gp_sources")
    g_ad.add_argument("--gp_sources_categories_mask_key", type=str,
                      default="nichecompass_gp_sources_categories")
    g_ad.add_argument("--latent_key", type=str, default="nichecompass_latent")
    g_ad.add_argument("--cat_covariates_keys", nargs="+", type=str, default=argparse.SUPPRESS,
                      help="obs keys for categorical covariates (default: [sample_key]).")

    return parser, pre


def merge_config_and_args(args: argparse.Namespace, cfg: dict[str, Any]) -> PreprocessParams:
    all_params = set(f.name for f in fields(PreprocessParams))
    allowed_for_config = sorted(list(all_params))
    if cfg:
        validate_known_keys(cfg, allowed_for_config)

    merged: dict[str, Any] = {}
    merged.update(cfg or {})

    for k, v in vars(args).items():
        if k == "config":
            continue
        if v is not None:
            merged[k] = v

    if "batches" not in merged or not merged["batches"]:
        raise ValueError("You must provide --batches PATH [PATH ...] or set 'batches' in --config JSON.")

    merged["batches"] = [Path(p) for p in merged["batches"]]

    species_val = merged.get("species", "human")
    if isinstance(species_val, str):
        merged["species"] = _validate_and_cast_species(species_val)
    else:
        raise TypeError("species must be a string.")

    pp = PreprocessParams(**merged)
    if pp.cat_covariates_keys is None:
        pp.cat_covariates_keys = [pp.sample_key]

    return pp


##########################################################
#### Prior knowledge GP mask ####

def create_prior_gp_mask(
    nichecompass_data_dir: Path,
    species: Species,
    figure_folder_path: Path,
) -> dict[str, Any]:
    logging.info(f"Preparing GP reference paths from {nichecompass_data_dir}")

    ga_data_folder_path = nichecompass_data_dir / "gene_annotations"
    gp_data_folder_path = nichecompass_data_dir / "gene_programs"

    omnipath_lr_network_file_path = gp_data_folder_path / "omnipath_lr_network.csv"
    nichenet_lr_network_file_path = gp_data_folder_path / f"nichenet_lr_network_v2_{species}.csv"
    nichenet_ligand_target_matrix_file_path = (
        gp_data_folder_path / f"nichenet_ligand_target_matrix_v2_{species}.csv"
    )
    mebocost_enzyme_sensor_interactions_folder_path = (
        gp_data_folder_path / "metabolite_enzyme_sensor_gps"
    )
    gene_orthologs_mapping_file_path = ga_data_folder_path / "human_mouse_gene_orthologs.csv"

    logging.info("Extracting OmniPath GP dict...")
    omnipath_gp_dict = extract_gp_dict_from_omnipath_lr_interactions(
        species=species,
        load_from_disk=True,
        save_to_disk=False,
        lr_network_file_path=str(omnipath_lr_network_file_path),
        gene_orthologs_mapping_file_path=str(gene_orthologs_mapping_file_path),
        plot_gp_gene_count_distributions=False,
        gp_gene_count_distributions_save_path=str(
            figure_folder_path / "omnipath_gp_gene_count_distributions.svg"
        ),
    )
    logging.info(f"OmniPath GP count: {len(omnipath_gp_dict)}")

    logging.info("Extracting NicheNet GP dict...")
    nichenet_gp_dict = extract_gp_dict_from_nichenet_lrt_interactions(
        species=species,
        version="v2",
        keep_target_genes_ratio=1.0,
        max_n_target_genes_per_gp=250,
        load_from_disk=False,
        save_to_disk=True,
        lr_network_file_path=str(nichenet_lr_network_file_path),
        ligand_target_matrix_file_path=str(nichenet_ligand_target_matrix_file_path),
        gene_orthologs_mapping_file_path=str(gene_orthologs_mapping_file_path),
        plot_gp_gene_count_distributions=False,
        gp_gene_count_distributions_save_path=str(
            figure_folder_path / "nichenet_gp_gene_count_distributions.svg"
        ),
    )
    logging.info(f"NicheNet GP count: {len(nichenet_gp_dict)}")

    logging.info("Extracting MEBOCOST GP dict...")
    mebocost_gp_dict = extract_gp_dict_from_mebocost_ms_interactions(
        dir_path=str(mebocost_enzyme_sensor_interactions_folder_path),
        species=species,
        plot_gp_gene_count_distributions=False,
        gp_gene_count_distributions_save_path=str(
            figure_folder_path / "mebocost_gp_gene_count_distributions.svg"
        ),
    )
    logging.info(f"MEBOCOST GP count: {len(mebocost_gp_dict)}")

    logging.info("Combining GP dicts...")
    combined_gp_dict = filter_and_combine_gp_dict_gps_v2(
        [omnipath_gp_dict, nichenet_gp_dict, mebocost_gp_dict],
        verbose=True,
    )
    logging.info(f"Combined GP count: {len(combined_gp_dict)}")
    return combined_gp_dict


##########################################################
#### Data loading and preparation ####

def load_batches(
    batch_paths: list[Path],
    counts_key: str,
) -> tuple[list[ad.AnnData], str]:
    """Load .h5ad/.zarr batches. Falls back to X for all if any batch is missing the counts layer."""
    adata_batch_list: list[ad.AnnData] = []
    use_x = False

    for p in batch_paths:
        logging.info(f"Loading batch: {p}")
        try:
            a = ad.read_zarr(p) if p.suffix == ".zarr" else sc.read_h5ad(p)
        except Exception as e:
            raise RuntimeError(f"Failed to read H5AD: {p}") from e

        if counts_key not in a.layers.keys():
            logging.warning(
                f"Layer '{counts_key}' not found in {p}; falling back to X for all batches."
            )
            use_x = True

        adata_batch_list.append(a)

    return adata_batch_list, ("X" if use_x else counts_key)


def compute_spatial_neighbors_for_batches(
    adata_batch_list: list[ad.AnnData],
    *,
    spatial_key: str,
    n_neighbors: int,
    adj_key: str,
) -> None:
    """Compute spatial neighbours and symmetrize adjacency in-place for each batch."""
    for idx, a in enumerate(adata_batch_list):
        logging.info(f"Computing spatial neighbours for batch {idx} (n_neighbors={n_neighbors})...")

        if spatial_key not in a.obsm_keys():
            raise KeyError(
                f"Missing obsm['{spatial_key}'] in batch #{idx}. "
                f"Provide --spatial_key or precompute spatial coordinates."
            )

        effective_n = min(n_neighbors, a.n_obs - 1)
        if effective_n < 1:
            logging.warning(f"Batch {idx} has only {a.n_obs} cell(s); skipping spatial neighbours.")
            import scipy.sparse as sp
            a.obsp[adj_key] = sp.csr_matrix((a.n_obs, a.n_obs))
            a.obsp[adj_key.replace("connectivities", "distances")] = sp.csr_matrix((a.n_obs, a.n_obs))
        else:
            if effective_n < n_neighbors:
                logging.warning(f"Batch {idx}: clamping n_neighbors {n_neighbors} → {effective_n} (only {a.n_obs} cells).")
            sq.gr.spatial_neighbors(a, coord_type="generic", spatial_key=spatial_key, n_neighs=effective_n)

        if adj_key not in a.obsp:
            raise KeyError(
                f"Expected obsp['{adj_key}'] after spatial_neighbors in batch #{idx}."
            )

        adj = a.obsp[adj_key]
        a.obsp[adj_key] = adj.maximum(adj.T)


def concat_batches_with_block_adj(
    adata_batch_list: list[ad.AnnData],
    *,
    adj_key: str,
    batch_key: str,
) -> ad.AnnData:
    """Concatenate batches (inner join on vars) with block-diagonal adjacency."""
    logging.info(f"Concatenating {len(adata_batch_list)} batches...")
    adata = ad.concat(adata_batch_list, join="inner", label=batch_key)

    logging.info("Assembling block-diagonal adjacency...")
    blocks = [a.obsp[adj_key] for a in adata_batch_list]
    adata.obsp[adj_key] = sp.block_diag(blocks, format="csr")
    return adata


def _replace_nonfinite_in_matrix(mat: Any) -> int:
    if sp.issparse(mat):
        data = mat.data
        mask = ~np.isfinite(data)
        count = int(mask.sum())
        if count:
            data[mask] = 0.0
            mat.eliminate_zeros()
        return count
    arr = np.asarray(mat)
    mask = ~np.isfinite(arr)
    count = int(mask.sum())
    if count:
        arr[mask] = 0.0
    return count


def sanitize_adata_for_training(
    adata: ad.AnnData,
    *,
    counts_key_effective: str,
    adj_key: str,
    spatial_key: str,
) -> None:
    """Replace NaN/Inf in expression and adjacency matrices; ensure float32 dtype."""
    if not adata.obs_names.is_unique:
        logging.warning("Observation names are not unique; making them unique.")
        adata.obs_names_make_unique()

    if spatial_key in adata.obsm:
        spatial = np.asarray(adata.obsm[spatial_key])
        if not np.isfinite(spatial).all():
            bad = int((~np.isfinite(spatial)).sum())
            raise ValueError(
                f"Non-finite values in obsm['{spatial_key}'] ({bad} entries). "
                "Please fix spatial coordinates."
            )

    if counts_key_effective == "X":
        mat = adata.X
        count = _replace_nonfinite_in_matrix(mat)
        if mat.dtype != np.float32:
            logging.warning("Casting expression matrix (X) from %s to float32.", mat.dtype)
            adata.X = mat.astype(np.float32)
    else:
        mat = adata.layers[counts_key_effective]
        count = _replace_nonfinite_in_matrix(mat)
        if mat.dtype != np.float32:
            logging.warning("Casting expression matrix (%s) from %s to float32.",
                            counts_key_effective, mat.dtype)
            adata.layers[counts_key_effective] = mat.astype(np.float32)
    if count:
        logging.warning("Replaced %d non-finite values in expression matrix (%s) with 0.",
                        count, counts_key_effective)

    if adj_key in adata.obsp:
        adj = adata.obsp[adj_key]
        count = _replace_nonfinite_in_matrix(adj)
        if count:
            logging.warning("Replaced %d non-finite values in adjacency (%s) with 0.",
                            count, adj_key)


def add_gp_masks_to_adata(
    adata: ad.AnnData,
    combined_gp_dict: dict[str, Any],
    *,
    gp_targets_mask_key: str,
    gp_targets_categories_mask_key: str,
    gp_sources_mask_key: str,
    gp_sources_categories_mask_key: str,
    gp_names_key: str,
) -> None:
    logging.info("Adding GP masks to AnnData...")
    add_gps_from_gp_dict_to_adata(
        gp_dict=combined_gp_dict,
        adata=adata,
        gp_targets_mask_key=gp_targets_mask_key,
        gp_targets_categories_mask_key=gp_targets_categories_mask_key,
        gp_sources_mask_key=gp_sources_mask_key,
        gp_sources_categories_mask_key=gp_sources_categories_mask_key,
        gp_names_key=gp_names_key,
        min_genes_per_gp=2,
        min_source_genes_per_gp=1,
        min_target_genes_per_gp=1,
        max_genes_per_gp=None,
        max_source_genes_per_gp=None,
        max_target_genes_per_gp=None,
    )
    logging.info("GP masks added.")


##########################################################

def main(argv: list[str] | None = None) -> None:
    parser, pre = build_parser()
    pre_args, remaining = pre.parse_known_args(argv)

    cfg_params: dict[str, Any] = {}
    if pre_args.config:
        cfg_params = load_config_json(pre_args.config)

    args = parser.parse_args(remaining)

    try:
        params = merge_config_and_args(args, cfg_params)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(2)

    fixed_seeds(0)

    setup_logging(Path("preprocess.log"), params.debug)
    logging.info("=== NicheCompass Preprocess: START ===")
    logging.info("Parameters: %s", json.dumps(asdict(params), indent=2, default=str))

    # Copy reference data from container image
    nichecompass_data_dir = Path("data")
    logging.info("Copying reference data from container image...")
    shutil.copytree("/app/data", nichecompass_data_dir)

    # Figures from GP mask creation go alongside the data dir
    gp_figure_dir = nichecompass_data_dir / "gp_figures"
    gp_figure_dir.mkdir(parents=True, exist_ok=True)

    # Build prior GP mask
    logging.info("Creating prior gene program mask...")
    combined_gp_dict = create_prior_gp_mask(
        nichecompass_data_dir=nichecompass_data_dir,
        species=params.species,
        figure_folder_path=gp_figure_dir,
    )

    # Load batches and build spatial graphs
    logging.info("Loading data batches...")
    adata_batch_list, counts_key_effective = load_batches(
        batch_paths=params.batches,
        counts_key=params.counts_key,
    )

    compute_spatial_neighbors_for_batches(
        adata_batch_list,
        spatial_key=params.spatial_key,
        n_neighbors=params.n_neighbors,
        adj_key=_ADJ_KEY,
    )

    adata = concat_batches_with_block_adj(
        adata_batch_list,
        adj_key=_ADJ_KEY,
        batch_key=params.sample_key,
    )

    sanitize_adata_for_training(
        adata,
        counts_key_effective=counts_key_effective,
        adj_key=_ADJ_KEY,
        spatial_key=params.spatial_key,
    )

    # Add GP masks
    add_gp_masks_to_adata(
        adata=adata,
        combined_gp_dict=combined_gp_dict,
        gp_targets_mask_key=params.gp_targets_mask_key,
        gp_targets_categories_mask_key=params.gp_targets_categories_mask_key,
        gp_sources_mask_key=params.gp_sources_mask_key,
        gp_sources_categories_mask_key=params.gp_sources_categories_mask_key,
        gp_names_key=params.gp_names_key,
    )

    # Store preprocessing context for the training step
    adata.uns["nichecompass_preprocess_params"] = {
        "counts_key_effective": counts_key_effective,
        "adj_key": _ADJ_KEY,
        "gp_names_key": params.gp_names_key,
        "active_gp_names_key": params.active_gp_names_key,
        "gp_targets_mask_key": params.gp_targets_mask_key,
        "gp_targets_categories_mask_key": params.gp_targets_categories_mask_key,
        "gp_sources_mask_key": params.gp_sources_mask_key,
        "gp_sources_categories_mask_key": params.gp_sources_categories_mask_key,
        "latent_key": params.latent_key,
        "sample_key": params.sample_key,
        "cell_type_key": params.cell_type_key,
        "cat_covariates_keys": params.cat_covariates_keys,
        "species": params.species,
        "n_neighbors": params.n_neighbors,
        "spatial_key": params.spatial_key,
    }

    # Write preprocessed AnnData
    out_path = Path(f"{params.prefix}_preprocessed.h5ad")
    logging.info(f"Writing preprocessed AnnData to {out_path}")
    adata.write_h5ad(out_path)

    logging.info("=== NicheCompass Preprocess: DONE ===")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
