#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nichecompass_train_sample_integration.py

Train NicheCompass (sample integration) on multiple H5AD batches and save the model
and artifacts into a timestamped folder. At the end, prints a single line:
  TIMESTAMP=YYYYMMDD_HHMMSS
Use this timestamp for the analysis notebook.

Notes:
  - Seed is hard-fixed to 0.
  - Training script does not display figures, but silently saves GP gene-count
    distribution SVGs to the figures folder.
"""

import argparse
import json
import logging
import random
import shutil
import sys
from dataclasses import asdict, dataclass, field, fields
from datetime import datetime
from pathlib import Path
from typing import Any, cast, Literal, TypeAlias

import anndata as ad
import numpy as np
import scanpy as sc
import squidpy as sq
import scipy.sparse as sp
import torch

#TODO: Need to fix downloading omnipathdb every single run
from nichecompass.models import NicheCompass
from nichecompass.utils import (
    add_gps_from_gp_dict_to_adata,
    filter_and_combine_gp_dict_gps_v2,
    extract_gp_dict_from_mebocost_ms_interactions,
    extract_gp_dict_from_nichenet_lrt_interactions,
    extract_gp_dict_from_omnipath_lr_interactions,
)

_UMAP_RANDOM_STATE: int = 0
Species: TypeAlias = Literal["human", "mouse"]
Touch_geometric: TypeAlias = Literal["gcnconv", "gatv2conv"]


def setup_logging(run_root: Path, debug: bool) -> None:
    log_path = run_root / "train.log"
    level = logging.DEBUG if debug else logging.INFO
    fmt = "%(asctime)s | %(levelname)s | %(message)s"
    datefmt = "%Y%m%d %H:%M:%S"
    logging.basicConfig(
        level=level,
        format=fmt,
        datefmt=datefmt,
        handlers=[logging.FileHandler(log_path, encoding="utf-8"),
                  logging.StreamHandler(sys.stdout)],
        force=True,
    )


def fixed_seeds(seed: int = 0) -> None:
    """
    Use fixed seeds for reproducibility.
    """
    global _UMAP_RANDOM_STATE
    _UMAP_RANDOM_STATE = seed

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if getattr(torch, "cuda", None) and torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


#### Dataclass and functions to parse parameters ####

@dataclass
class RunParams:
    """
    dataclass for storing parameters (user-defined & runtime)
    """
    # MAIN
    batches: list[Path]
    outdir: Path
    species: Species = "human"
    nichecompass_version: str = "0.3.0"
    debug: bool = False

    # DATASET / GRAPH
    spatial_key: str = "spatial"
    n_neighbors: int = 4
    sample_key: str = "batch"
    cell_type_key: str = "Main_molecular_cell_type"

    # AnnData keys
    counts_key: str = "counts"
    adj_key: str = field(default="spatial_connectivities", init=False)
    gp_names_key: str = "nichecompass_gp_names"
    active_gp_names_key: str = "nichecompass_active_gp_names"
    gp_targets_mask_key: str = "nichecompass_gp_targets"
    gp_targets_categories_mask_key: str = "nichecompass_gp_targets_categories"
    gp_sources_mask_key: str = "nichecompass_gp_sources"
    gp_sources_categories_mask_key: str = "nichecompass_gp_sources_categories"
    latent_key: str = "nichecompass_latent"

    # MODEL / ARCH
    cat_covariates_keys: list[str] | None = None  # default to [sample_key] later
    cat_covariates_embeds_injection: list[str] = field(default_factory=lambda: ["gene_expr_decoder"])
    cat_covariates_embeds_nums: list[int] = field(default_factory=lambda: [3])
    cat_covariates_no_edges: list[bool] = field(default_factory=lambda: [True])
    conv_layer_encoder: Touch_geometric = "gatv2conv"
    active_gp_thresh_ratio: float = 0.01

    # TRAINER
    n_epochs: int = 400
    n_epochs_all_gps: int = 25
    lr: float = 1e-3
    lambda_edge_recon: float = 500000.0
    lambda_gene_expr_recon: float = 300.0
    lambda_l1_masked: float = 0.0
    lambda_l1_addon: float = 30.0
    edge_batch_size: int = 16384
    n_sampled_neighbors: int = 4
    use_cuda_if_available: bool = True

    # Derived at runtime
    timestamp: str | None = None
    run_root: Path | None = None
    artifacts_folder_path: Path | None = None
    model_folder_path: Path | None = None
    figure_folder_path: Path | None = None
    nichecompass_data_dir: Path | None = None

    def finalize_paths(self) -> None:
        self.timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.run_root = self.outdir / f"nichecompass_{self.timestamp}"
        self.nichecompass_data_dir = self.run_root / "data"
        self.artifacts_folder_path = self.run_root / "artifacts"
        self.model_folder_path = self.artifacts_folder_path / "model"
        self.figure_folder_path = self.artifacts_folder_path / "figures"


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


def normalise_list_arg(val: list[Any] | None, *, expected_len: int, default_item: Any) -> list[Any]:
    """
    Normalize a possibly missing or length-mismatched list argument to a fixed length.

    This helper is used for options that conceptually take one value per category
    (e.g., one setting per covariate). It ensures the output list has exactly
    `expected_len` items using the following rules:

    - If ``val is None``: return ``[default_item] * expected_len``.
    - If ``expected_len <= 0``: return ``[]``.
    - If ``val`` is an empty list: warn and return ``[default_item] * expected_len``.
    - If ``len(val) == expected_len``: return ``val`` unchanged.
    - If ``len(val) == 1``: repeat the single value to reach ``expected_len``.
    - Otherwise (length mismatch): log a warning and
      *truncate* to ``expected_len`` if too long, or *pad* by repeating the last
      element if too short.

    Parameters
    ----------
    val
        The user-provided list (from CLI/config) or ``None``.
    expected_len
        Desired length of the normalized list (e.g., number of covariates).
    default_item
        Filler value used when ``val`` is ``None`` or requires padding.

    Returns
    -------
    list[Any]
        A list of length exactly ``expected_len``.
    """
    if val is None:
        return [default_item] * expected_len
    if expected_len <= 0:
        return []
    if len(val) == expected_len:
        return val
    if len(val) == 1:
        return [val[0]] * expected_len
    logging.warning(f"List length {len(val)} does not match expected {expected_len}; adjusting.")
    out = list(val)[:expected_len]
    while len(out) < expected_len:
        out.append(val[-1])
    return out


def build_parser() -> tuple[argparse.ArgumentParser, argparse.ArgumentParser]:
    """
    Argument parser
    """
    # 'pre' parser for parsing params from json config
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config", type=Path, default=None, help="Path to JSON config file.")

    # 'parser' for parsing params from CLI
    parser = argparse.ArgumentParser(
        prog="nichecompass_train_sample_integration.py",
        description=(
            "Train NicheCompass (sample integration) on multiple H5AD batches and save the model + artifacts "
            "into a timestamped folder. Prints TIMESTAMP=YYYYMMDD_HHMMSS at the end."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[pre],
    )

    # Grouping parameters
    g_main = parser.add_argument_group("MAIN (I/O & run identity)")
    g_main.add_argument("--batches", nargs="+", type=Path, default=argparse.SUPPRESS, help="Paths to input .h5ad files (≥1).")
    g_main.add_argument("--outdir", type=Path, default=argparse.SUPPRESS, help="Base output directory (default: current working directory).")
    g_main.add_argument("--species", type=str, choices=["human", "mouse"], default="human", help="Species tag for prior knowledge lookup.")
    g_main.add_argument("--debug", action="store_true", help="Enable DEBUG logging.")

    g_dataset = parser.add_argument_group("DATASET / GRAPH")
    g_dataset.add_argument("--sample_key", type=str, default="batch", help="obs key for sample/batch.")
    g_dataset.add_argument("--cell_type_key", type=str, default="Main_molecular_cell_type",
                           help="obs key for cell type labels.")
    g_dataset.add_argument("--spatial_key", type=str, default="spatial", help="obsm key for spatial coordinates.")
    g_dataset.add_argument("--n_neighbors", type=int, default=4, help="Number of spatial neighbors per node.")

    g_ad = parser.add_argument_group("AnnData keys")
    g_ad.add_argument("--counts_key", type=str, default="counts", help="layer name for counts (falls back to X).")
    # g_ad.add_argument("--adj_key", type=str, default="spatial_connectivities", help="obsp key for spatial adjacency matrix.")
    g_ad.add_argument("--gp_names_key", type=str, default="nichecompass_gp_names", help="uns key for gene program names.")
    g_ad.add_argument("--active_gp_names_key", type=str, default="nichecompass_active_gp_names", help="uns key for active gene program names.")
    g_ad.add_argument("--gp_targets_mask_key", type=str, default="nichecompass_gp_targets", help="varm key for gene program targets mask.")
    g_ad.add_argument("--gp_targets_categories_mask_key", type=str, default="nichecompass_gp_targets_categories", help="varm key for gene program targets categories mask.")
    g_ad.add_argument("--gp_sources_mask_key", type=str, default="nichecompass_gp_sources", help="varm key for gene program sources mask.")
    g_ad.add_argument("--gp_sources_categories_mask_key", type=str, default="nichecompass_gp_sources_categories", help="varm key for gene program sources categories mask.")
    g_ad.add_argument("--latent_key", type=str, default="nichecompass_latent", help="obsm key for latent / gene program representation of active gene programs after model training.")

    g_model = parser.add_argument_group("MODEL / ARCHITECTURE")
    g_model.add_argument("--cat_covariates_keys", nargs="+", type=str, default=argparse.SUPPRESS,
                         help="obs key for categorical covariates (default: [sample_key]).")
    g_model.add_argument("--cat_covariates_embeds_injection", nargs="+", type=str,
                         default=["gene_expr_decoder"], help="List of VGPGAE modules in which the categorical covariates embeddings are injected.")
    g_model.add_argument("--cat_covariates_embeds_nums", nargs="+", type=int, default=[3],
                         help="List of number of embedding nodes for all categorical covariates.")
    g_model.add_argument("--cat_covariates_no_edges", nargs="+", type=str, default=["true"],
                         help="List of booleans that indicate whether there can be edges between different categories of the categorical covariates.")
    g_model.add_argument("--conv_layer_encoder", type=str, choices=["gcnconv", "gatv2conv"], default="gatv2conv", help="Encoder conv layer type.")
    g_model.add_argument("--active_gp_thresh_ratio", type=float, default=0.01,
                         help="Threshold ratio for active GP selection.")

    g_tr = parser.add_argument_group("TRAINER")
    g_tr.add_argument("--n_epochs", type=int, default=400, help="Total training epochs.")
    g_tr.add_argument("--n_epochs_all_gps", type=int, default=25, help="Warmup epochs training all GPs.")
    g_tr.add_argument("--lr", type=float, default=1e-3, help="Learning rate.")
    g_tr.add_argument("--lambda_edge_recon", type=float, default=500000.0, help="Edge reconstruction weight.")
    g_tr.add_argument("--lambda_gene_expr_recon", type=float, default=300.0, help="Gene expression reconstruction weight.")
    g_tr.add_argument("--lambda_l1_masked", type=float, default=0.0, help="L1 regularization on masked GPs.")
    g_tr.add_argument("--lambda_l1_addon", type=float, default=30.0, help="L1 regularization on addon GPs.")
    g_tr.add_argument("--edge_batch_size", type=int, default=16384, help="Edge batch size.")
    g_tr.add_argument("--n_sampled_neighbors", type=int, default=4, help="Number of sampled neighbors.")
    g_tr.add_argument("--use_cuda_if_available", type=str, choices=["true", "false"], default="true",
                      help="Use CUDA if available.")
    return parser, pre


def merge_config_and_args(args: argparse.Namespace, cfg: dict[str, Any]) -> RunParams:
    """
    Merge params from config.json and CLI

    Args:
        args (argparse.Namespace): params from CLI
        cfg (dict[str, Any]): params from config

    Raises:
        ValueError: if required params (batches) is not provided

    Returns:
        RunParams: dataclass with all params
    """

    all_params = set(f.name for f in fields(RunParams))
    runtime_params = {"timestamp", "run_root", "artifacts_folder_path", "model_folder_path",
               "figure_folder_path", "nichecompass_data_dir"}
    sealed_params = {"adj_key"}
    
    # Extract the user-defined params and validate the params from config
    allowed_for_config = sorted(list(all_params - runtime_params - sealed_params))
    if cfg:
        validate_known_keys(cfg, allowed_for_config)

    # Merge params from config and CLI
    merged: dict[str, Any] = {}
    merged.update(cfg or {})

    # Overwrite config params with params from CLI 
    for k, v in vars(args).items():
        # Ignore --config as it's already used
        if k == "config":
            continue
        if v is not None:
            if k == "use_cuda_if_available" and isinstance(v, str):
                merged[k] = (v.lower() == "true")
            elif k == "cat_covariates_no_edges" and isinstance(v, list):
                merged[k] = [str2bool(x) if isinstance(x, str) else bool(x) for x in v]
            else:
                merged[k] = v

    # Check required batches params are provided
    if "batches" not in merged or not merged["batches"]:
        raise ValueError("You must provide --batches PATH [PATH ...] or set 'batches' in --config JSON.")

    # Normalise to Path types
    merged["batches"] = [Path(p) for p in merged["batches"]]
    merged["outdir"] = Path(merged["outdir"]) if isinstance(merged.get("outdir"), (str, Path)) else Path.cwd()

    # Validate  species
    species_val = merged.get("species", "human")
    if isinstance(species_val, str):
        merged["species"] = _validate_and_cast_species(species_val)
    else:
        raise TypeError("species must be a string.")

    # Construct typed dataclass
    rp = RunParams(**merged)

    # Set [sample_key] as default for cat_covariates_keys
    if rp.cat_covariates_keys is None:
        rp.cat_covariates_keys = [rp.sample_key]

    return rp


#####################################################

#### Functions to create Prior Knowledge Gene Program (GP) Mask ####
def create_prior_gp_mask(
    nichecompass_data_dir: Path,
    species: Species,
    figure_folder_path: Path,
) -> dict[str, Any]:
    """
    Create the prior knowledge GP mask (combined) from OmniPath, NicheNet, and MEBOCOST.

    - Always reads reference files from `nichecompass_data_dir` (downloaded earlier).
    - Disables plotting but still writes the SVG summaries into `figure_folder_path`.
    - Returns the combined GP dictionary to be passed into `add_gps_from_gp_dict_to_adata`.

    Parameters
    ----------
    nichecompass_data_dir : Path
        Root path containing the 'data/' folder fetched from the NicheCompass repo.
    species : Species
        Species tag ('human' by default) to select correct reference files.
    figure_folder_path : Path
        Where to save *_gp_gene_count_distributions.svg.

    Returns
    -------
    dict[str, Any]
        Combined GP dictionary after filtering/merging.
    """
    logging.info(f"Preparing GP reference paths from {nichecompass_data_dir}")

    ga_data_folder_path = nichecompass_data_dir / "gene_annotations"
    gp_data_folder_path = nichecompass_data_dir / "gene_programs"

    omnipath_lr_network_file_path = gp_data_folder_path / "omnipath_lr_network.csv"
    nichenet_lr_network_file_path = gp_data_folder_path / f"nichenet_lr_network_v2_{species}.csv"
    nichenet_ligand_target_matrix_file_path = gp_data_folder_path / f"nichenet_ligand_target_matrix_v2_{species}.csv"
    mebocost_enzyme_sensor_interactions_folder_path = gp_data_folder_path / "metabolite_enzyme_sensor_gps"
    gene_orthologs_mapping_file_path = ga_data_folder_path / "human_mouse_gene_orthologs.csv"

    logging.info("Extracting OmniPath GP dict…")
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

    logging.info("Extracting NicheNet GP dict…")
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

    logging.info("Extracting MEBOCOST GP dict…")
    mebocost_gp_dict = extract_gp_dict_from_mebocost_ms_interactions(
        dir_path=str(mebocost_enzyme_sensor_interactions_folder_path),
        species=species,
        plot_gp_gene_count_distributions=False,
        gp_gene_count_distributions_save_path=str(
            figure_folder_path / "mebocost_gp_gene_count_distributions.svg"
        ),
    )
    logging.info(f"MEBOCOST GP count: {len(mebocost_gp_dict)}")

    logging.info("Combining GP dicts…")
    combined_gp_dict = filter_and_combine_gp_dict_gps_v2(
        [omnipath_gp_dict, nichenet_gp_dict, mebocost_gp_dict],
        verbose=True,
    )
    logging.info(
        f"Number of gene programs after filtering and combining: {len(combined_gp_dict)}",
    )
    return combined_gp_dict

#####################################################

#### Functions to prepare training data ####
def load_batches(
    batch_paths: list[Path],
    counts_key: str,
) -> tuple[list[ad.AnnData], str]:
    """
    Load .h5ad batches. 
    If any batch is missing the requested counts layer, fall back to 'X' for all
        (graceful, consistent with earlier behavior).
    Returns the list of AnnData and the effective counts key.
    """
    adata_batch_list: list[ad.AnnData] = []
    use_x = False

    for p in batch_paths:
        logging.info(f"Loading batch: {p}")
        try:
            a = sc.read_h5ad(p)
        except Exception as e:
            raise RuntimeError(f"Failed to read H5AD: {p}") from e


        if counts_key not in a.layers.keys():
            logging.warning(f"Layer '{counts_key}' not found in {p}; falling back to X for all batches.")
            use_x = True

        adata_batch_list.append(a)

    counts_key_effective = "X" if use_x else counts_key

    return adata_batch_list, counts_key_effective


def compute_spatial_neighbors_for_batches(
    adata_batch_list: list[ad.AnnData],
    *,
    spatial_key: str,
    n_neighbors: int,
    adj_key: str,
) -> None:
    """
    For each batch:
      - verify spatial coords exist,
      - compute spatial neighbors with squidpy,
      - ensure obsp[adj_key] exists,
      - symmetrize adjacency.
    Modifies AnnData objects in-place.
    """
    for idx, a in enumerate(adata_batch_list):
        logging.info(f"Computing spatial neighbors for batch {idx} (n_neighbors={n_neighbors})...")

        if spatial_key not in a.obsm_keys():
            raise KeyError(
                f"Missing obsm['{spatial_key}'] in batch #{idx}. Provide --spatial_key or precompute spatial coordinates."
            )

        sq.gr.spatial_neighbors(
            a,
            coord_type="generic",
            spatial_key=spatial_key,
            n_neighs=n_neighbors,
        )

        if adj_key not in a.obsp:
            raise KeyError(
                f"Expected obsp['{adj_key}'] after spatial_neighbors in batch #{idx}. "
                f"Check that squidpy writes to this key."
            )

        # Make adjacency symmetric
        adj = a.obsp[adj_key]
        a.obsp[adj_key] = adj.maximum(adj.T)


def concat_batches_with_block_adj(
    adata_batch_list: list[ad.AnnData],
    *,
    adj_key: str,
) -> ad.AnnData:
    """
    Concatenate batches (inner join on vars) and set a block-diagonal adjacency in obsp[adj_key].
    """
    logging.info(f"Concatenating {len(adata_batch_list)} batches...")
    adata = ad.concat(adata_batch_list, join="inner")

    logging.info("Assembling block-diagonal adjacency…")
    blocks = [a.obsp[adj_key] for a in adata_batch_list]
    # sp.block_diag returns coo by default; ask for csr
    adata.obsp[adj_key] = sp.block_diag(blocks, format="csr")
    return adata


def add_gp_masks_to_adata(
    adata: ad.AnnData,
    combined_gp_dict: dict[str, Any],
    params: RunParams,
) -> None:
    """
    Add the combined prior GP masks into AnnData in-place.
    """
    logging.info("Adding GP masks to AnnData...")
    add_gps_from_gp_dict_to_adata(
        gp_dict=combined_gp_dict,
        adata=adata,
        gp_targets_mask_key=params.gp_targets_mask_key,
        gp_targets_categories_mask_key=params.gp_targets_categories_mask_key,
        gp_sources_mask_key=params.gp_sources_mask_key,
        gp_sources_categories_mask_key=params.gp_sources_categories_mask_key,
        gp_names_key=params.gp_names_key,
        min_genes_per_gp=2,
        min_source_genes_per_gp=1,
        min_target_genes_per_gp=1,
        max_genes_per_gp=None,
        max_source_genes_per_gp=None,
        max_target_genes_per_gp=None,
    )
    logging.info("GP masks added.")


#####################################################

#### Functions to train nichecompass model ####
def build_model(
    adata: ad.AnnData,
    params: RunParams,
    counts_key_effective: str,
) -> NicheCompass:
    """
    Construct a NicheCompass model from AnnData and run parameters.
    """
    logging.info("Initializing NicheCompass model...")
    cat_keys = params.cat_covariates_keys or [params.sample_key]
    n_cov = len(cat_keys)

    inj = normalise_list_arg(
        params.cat_covariates_embeds_injection,
        expected_len=n_cov,
        default_item="gene_expr_decoder",
    )
    emb_dims = normalise_list_arg(
        params.cat_covariates_embeds_nums,
        expected_len=n_cov,
        default_item=3,
    )
    no_edges = normalise_list_arg(
        params.cat_covariates_no_edges,
        expected_len=n_cov,
        default_item=True,
    )

    model = NicheCompass(
        adata,
        counts_key=counts_key_effective,
        adj_key=params.adj_key,
        cat_covariates_embeds_injection=inj,
        cat_covariates_keys=cat_keys,
        cat_covariates_no_edges=no_edges,
        cat_covariates_embeds_nums=emb_dims,
        gp_names_key=params.gp_names_key,
        active_gp_names_key=params.active_gp_names_key,
        gp_targets_mask_key=params.gp_targets_mask_key,
        gp_targets_categories_mask_key=params.gp_targets_categories_mask_key,
        gp_sources_mask_key=params.gp_sources_mask_key,
        gp_sources_categories_mask_key=params.gp_sources_categories_mask_key,
        latent_key=params.latent_key,
        conv_layer_encoder=params.conv_layer_encoder,
        active_gp_thresh_ratio=params.active_gp_thresh_ratio,
    )
    logging.info("Model initialized.")
    return model


def train_and_embed(
    model: NicheCompass,
    params: RunParams,
) -> None:
    """
    Train the model, then compute neighbors and UMAP in latent space.
    """
    logging.info("Training model...")
    model.train(
        n_epochs=params.n_epochs,
        n_epochs_all_gps=params.n_epochs_all_gps,
        lr=params.lr,
        lambda_edge_recon=params.lambda_edge_recon,
        lambda_gene_expr_recon=params.lambda_gene_expr_recon,
        lambda_l1_masked=params.lambda_l1_masked,
        lambda_l1_addon=params.lambda_l1_addon,
        edge_batch_size=params.edge_batch_size,
        n_sampled_neighbors=params.n_sampled_neighbors,
        use_cuda_if_available=params.use_cuda_if_available,
        verbose=False,
    )
    logging.info("Training done. Computing neighbors and UMAP in latent space...")
    sc.pp.neighbors(model.adata, use_rep=params.latent_key, key_added=params.latent_key)
    sc.tl.umap(model.adata, neighbors_key=params.latent_key, random_state=_UMAP_RANDOM_STATE)
    logging.info("Neighbors/UMAP computed.")


def save_model_and_config(
    model: NicheCompass,
    params: RunParams,
    counts_key_effective: str,
) -> None:
    """
    Save trained model and run_config.json.
    """
    logging.info(f"Saving trained model to {params.model_folder_path}")
    model.save(
        dir_path=str(params.model_folder_path),
        overwrite=True,
        save_adata=True,
        adata_file_name="model.h5ad",
    )

    run_config = asdict(params)
    run_config["counts_key_effective"] = counts_key_effective
    cfg_path = params.run_root / "run_config.json"
    with open(cfg_path, "w", encoding="utf-8") as f:
        json.dump(run_config, f, indent=2, default=str)
    logging.info(f"Saved run_config.json at {cfg_path}")


def main(argv: list[str] | None = None) -> None:
    ### 1. Parse parameters ###
    parser, pre = build_parser()
    pre_args, remaining = pre.parse_known_args(argv)

    # if --config is provided, load params from json file
    cfg_params: dict[str, Any] = {}
    if pre_args.config:
        cfg_params = load_config_json(pre_args.config)

    # parse the params provided from other options
    args = parser.parse_args(remaining)

    # Merge the params provided from --config and other options
    try:
        params = merge_config_and_args(args, cfg_params)
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(2)

    ###########################

    # Fix seed for reproducibility
    fixed_seeds(0)

    ###########################

    ### 2. Create output directories and set up loggers ###
    # Set output directories paths
    params.finalize_paths()

    params.run_root.mkdir(parents=True, exist_ok=True)
    setup_logging(params.run_root, params.debug)
    logging.info("=== NicheCompass Sample Integration: START ===")
    logging.info("Resolved parameters: %s", json.dumps(asdict(params), indent=2, default=str))

    # Create artifacts directories
    params.figure_folder_path.mkdir(parents=True, exist_ok=True)
    params.model_folder_path.mkdir(parents=True, exist_ok=True)


    ### 3. Prepare prior knowledge gene program (GP) mask
    # Copy reference data dir from image
    logging.info("Copy reference data from container image...")
    shutil.copytree("/app/data", params.nichecompass_data_dir)

    logging.info("Creating prior gene program mask...")
    combined_gp_dict = create_prior_gp_mask(
        nichecompass_data_dir=params.nichecompass_data_dir,
        species=params.species,
        figure_folder_path=params.figure_folder_path,
    )
    logging.info(f"Number of gene programs after filtering and combining: {len(combined_gp_dict)}")


    ### 4. Load data & build spatial graphs
    logging.info("Loading data batches...")
    adata_batch_list, counts_key_effective = load_batches(
        batch_paths=params.batches,
        counts_key=params.counts_key,
    )

    compute_spatial_neighbors_for_batches(
        adata_batch_list,
        spatial_key=params.spatial_key,
        n_neighbors=params.n_neighbors,
        adj_key=params.adj_key,
    )

    adata = concat_batches_with_block_adj(
        adata_batch_list,
        adj_key=params.adj_key,
    )

    ### 5. Add GP masks to data
    add_gp_masks_to_adata(
        adata=adata,
        combined_gp_dict=combined_gp_dict,
        params=params,
    )

    ### 6. Build, train, and embed model
    model = build_model(
        adata=adata,
        params=params,
        counts_key_effective=counts_key_effective,
    )
    train_and_embed(model, params)

    # 7. Save outputs
    save_model_and_config(
        model=model,
        params=params,
        counts_key_effective=counts_key_effective,
    )

    logging.info("=== NicheCompass Sample Integration: DONE ===")
    timestamp_path = params.run_root / "timestamp.txt"
    with open(timestamp_path, 'w', encoding="utf-8") as f:
        f.write(f"TIMESTAMP={params.timestamp}")
    print(f"TIMESTAMP={params.timestamp}", file=sys.stdout)
    sys.stdout.flush()


if __name__ == "__main__":
    main()
