#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
nichecompass_train.py

Train a NicheCompass model on a preprocessed H5AD produced by nichecompass_preprocess.py.

Steps:
  1. Load the preprocessed AnnData and read preprocessing context from adata.uns
  2. Build the NicheCompass model
  3. Train the model
  4. Compute neighbours and UMAP in latent space
  5. Save the trained model and a run_config.json compatible with the analysis step

Usage:
  nichecompass_train.py --preprocessed_h5ad sample_preprocessed.h5ad \\
                        --model_dir sample_nichecompass_dir [options]

Outputs written under ./<model_dir>/:
  artifacts/model/   — trained model files including trained.h5ad
  artifacts/figures/ — optional training figures
  run_config.json    — full parameter record for the analysis step
"""

import argparse
import json
import logging
import random
import sys
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import Any, Literal, TypeAlias

import anndata as ad
import numpy as np
import scanpy as sc
import torch

from nichecompass.models import NicheCompass

_UMAP_RANDOM_STATE: int = 0

Touch_geometric: TypeAlias = Literal["gcnconv", "gatv2conv"]


def _json_default(obj: Any) -> Any:
    """JSON serializer that converts numpy scalars/arrays to Python natives.

    anndata loads uns values from h5ad as numpy types; without this, json.dump
    with default=str would turn array(['batch']) into the string "['batch']",
    which breaks downstream consumers expecting a JSON array.
    """
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    return str(obj)


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
    global _UMAP_RANDOM_STATE
    _UMAP_RANDOM_STATE = seed
    import os
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if getattr(torch, "cuda", None) and torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)


def last_n_levels(path_str: Any, n: int) -> Path:
    p = Path(path_str)
    return Path(*p.parts[-n:])


#### Parameter dataclass and parsing ####

@dataclass
class TrainParams:
    # MAIN
    preprocessed_h5ad: Path
    model_dir: str
    outdir: Path = field(default_factory=Path.cwd)
    debug: bool = False

    # MODEL / ARCHITECTURE
    cat_covariates_keys: list[str] | None = None       # default: read from preprocess uns
    cat_covariates_embeds_injection: list[str] = field(
        default_factory=lambda: ["gene_expr_decoder"]
    )
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

    # Runtime paths (set by finalize_paths)
    run_root: Path | None = None
    artifacts_folder_path: Path | None = None
    model_folder_path: Path | None = None
    figure_folder_path: Path | None = None

    def finalize_paths(self) -> None:
        self.run_root = self.outdir / self.model_dir
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


def str2bool(v: str) -> bool:
    if isinstance(v, bool):
        return v
    val = v.strip().lower()
    if val in {"true", "t", "1", "yes", "y"}:
        return True
    if val in {"false", "f", "0", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError(f"Invalid bool: {v!r}")


def normalise_list_arg(
    val: list[Any] | None, *, expected_len: int, default_item: Any
) -> list[Any]:
    if val is None:
        return [default_item] * expected_len
    if expected_len <= 0:
        return []
    if len(val) == expected_len:
        return val
    if len(val) == 1:
        return [val[0]] * expected_len
    logging.warning("List length %d does not match expected %d; adjusting.", len(val), expected_len)
    out = list(val)[:expected_len]
    while len(out) < expected_len:
        out.append(val[-1])
    return out


def build_parser() -> tuple[argparse.ArgumentParser, argparse.ArgumentParser]:
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("--config", type=Path, default=None, help="Path to JSON config file.")

    parser = argparse.ArgumentParser(
        prog="nichecompass_train.py",
        description=(
            "Train a NicheCompass model on a preprocessed H5AD. "
            "Preprocessing context is read from adata.uns['nichecompass_preprocess_params']."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        parents=[pre],
    )

    g_main = parser.add_argument_group("MAIN")
    g_main.add_argument("--preprocessed_h5ad", type=Path, required=True,
                        help="Path to the preprocessed H5AD from nichecompass_preprocess.py.")
    g_main.add_argument("--model_dir", type=str, required=True,
                        help="Output directory name (created under --outdir).")
    g_main.add_argument("--outdir", type=Path, default=argparse.SUPPRESS,
                        help="Base output directory (default: current working directory).")
    g_main.add_argument("--debug", action="store_true", help="Enable DEBUG logging.")

    g_model = parser.add_argument_group("MODEL / ARCHITECTURE")
    g_model.add_argument("--cat_covariates_keys", nargs="+", type=str, default=argparse.SUPPRESS,
                         help="obs keys for categorical covariates (default: [sample_key] from preprocess).")
    g_model.add_argument("--cat_covariates_embeds_injection", nargs="+", type=str,
                         default=["gene_expr_decoder"])
    g_model.add_argument("--cat_covariates_embeds_nums", nargs="+", type=int, default=[3])
    g_model.add_argument("--cat_covariates_no_edges", nargs="+", type=str, default=["true"],
                         help="Per-covariate bool: no edges between different categories.")
    g_model.add_argument("--conv_layer_encoder", type=str, choices=["gcnconv", "gatv2conv"],
                         default="gatv2conv")
    g_model.add_argument("--active_gp_thresh_ratio", type=float, default=0.01)

    g_tr = parser.add_argument_group("TRAINER")
    g_tr.add_argument("--n_epochs", type=int, default=400)
    g_tr.add_argument("--n_epochs_all_gps", type=int, default=25)
    g_tr.add_argument("--lr", type=float, default=1e-3)
    g_tr.add_argument("--lambda_edge_recon", type=float, default=500000.0)
    g_tr.add_argument("--lambda_gene_expr_recon", type=float, default=300.0)
    g_tr.add_argument("--lambda_l1_masked", type=float, default=0.0)
    g_tr.add_argument("--lambda_l1_addon", type=float, default=30.0)
    g_tr.add_argument("--edge_batch_size", type=int, default=16384)
    g_tr.add_argument("--n_sampled_neighbors", type=int, default=4)
    g_tr.add_argument("--use_cuda_if_available", type=str, choices=["true", "false"],
                       default="true")

    return parser, pre


def merge_config_and_args(args: argparse.Namespace, cfg: dict[str, Any]) -> TrainParams:
    all_params = set(f.name for f in fields(TrainParams))
    runtime_params = {"run_root", "artifacts_folder_path", "model_folder_path", "figure_folder_path"}
    allowed_for_config = sorted(list(all_params - runtime_params))
    if cfg:
        validate_known_keys(cfg, allowed_for_config)

    merged: dict[str, Any] = {}
    merged.update(cfg or {})

    for k, v in vars(args).items():
        if k == "config":
            continue
        if v is not None:
            if k == "use_cuda_if_available" and isinstance(v, str):
                merged[k] = v.lower() == "true"
            elif k == "cat_covariates_no_edges" and isinstance(v, list):
                merged[k] = [str2bool(x) if isinstance(x, str) else bool(x) for x in v]
            else:
                merged[k] = v

    if "outdir" not in merged:
        merged["outdir"] = Path.cwd()
    else:
        merged["outdir"] = Path(merged["outdir"])

    return TrainParams(**merged)


##########################################################
#### Model construction, training, saving ####

def build_model(
    adata: ad.AnnData,
    params: TrainParams,
    preprocess_ctx: dict[str, Any],
    counts_key_effective: str,
) -> NicheCompass:
    logging.info("Initialising NicheCompass model...")

    cat_keys: list[str] = [
        str(k) for k in (
            params.cat_covariates_keys
            or preprocess_ctx.get("cat_covariates_keys")
            or [preprocess_ctx["sample_key"]]
        )
    ]
    n_cov = len(cat_keys)

    inj = normalise_list_arg(params.cat_covariates_embeds_injection,
                             expected_len=n_cov, default_item="gene_expr_decoder")
    emb_dims = normalise_list_arg(params.cat_covariates_embeds_nums,
                                  expected_len=n_cov, default_item=3)
    no_edges = normalise_list_arg(params.cat_covariates_no_edges,
                                  expected_len=n_cov, default_item=True)

    model = NicheCompass(
        adata,
        counts_key=counts_key_effective,
        adj_key=preprocess_ctx["adj_key"],
        cat_covariates_embeds_injection=inj,
        cat_covariates_keys=cat_keys,
        cat_covariates_no_edges=no_edges,
        cat_covariates_embeds_nums=emb_dims,
        gp_names_key=preprocess_ctx["gp_names_key"],
        active_gp_names_key=preprocess_ctx["active_gp_names_key"],
        gp_targets_mask_key=preprocess_ctx["gp_targets_mask_key"],
        gp_targets_categories_mask_key=preprocess_ctx["gp_targets_categories_mask_key"],
        gp_sources_mask_key=preprocess_ctx["gp_sources_mask_key"],
        gp_sources_categories_mask_key=preprocess_ctx["gp_sources_categories_mask_key"],
        latent_key=preprocess_ctx["latent_key"],
        conv_layer_encoder=params.conv_layer_encoder,
        active_gp_thresh_ratio=params.active_gp_thresh_ratio,
    )
    logging.info("Model initialised.")
    return model


def train_and_embed(model: NicheCompass, params: TrainParams, latent_key: str) -> None:
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
    logging.info("Training done. Computing neighbours and UMAP in latent space...")
    sc.pp.neighbors(model.adata, use_rep=latent_key, key_added=latent_key)
    sc.tl.umap(model.adata, neighbors_key=latent_key, random_state=_UMAP_RANDOM_STATE)
    logging.info("Neighbours/UMAP computed.")


def save_model_and_config(
    model: NicheCompass,
    params: TrainParams,
    preprocess_ctx: dict[str, Any],
    counts_key_effective: str,
    cat_covariates_keys: list[str],
) -> None:
    logging.info(f"Saving trained model to {params.model_folder_path}")
    model.save(
        dir_path=str(params.model_folder_path),
        overwrite=True,
        save_adata=True,
        adata_file_name="trained.h5ad",
    )

    cfg_path = params.run_root / "run_config.json"

    # Build relative paths for portability (analysis notebook resolves these
    # relative to its own working directory where Nextflow stages the outputs)
    run_config = {
        # Preprocessing context (needed by analysis notebook)
        "counts_key_effective": counts_key_effective,
        "adj_key": preprocess_ctx["adj_key"],
        "gp_names_key": preprocess_ctx["gp_names_key"],
        "active_gp_names_key": preprocess_ctx["active_gp_names_key"],
        "gp_targets_mask_key": preprocess_ctx["gp_targets_mask_key"],
        "gp_targets_categories_mask_key": preprocess_ctx["gp_targets_categories_mask_key"],
        "gp_sources_mask_key": preprocess_ctx["gp_sources_mask_key"],
        "gp_sources_categories_mask_key": preprocess_ctx["gp_sources_categories_mask_key"],
        "latent_key": preprocess_ctx["latent_key"],
        "sample_key": preprocess_ctx["sample_key"],
        "cell_type_key": preprocess_ctx["cell_type_key"],
        "cat_covariates_keys": cat_covariates_keys,
        "species": preprocess_ctx["species"],
        "n_neighbors": preprocess_ctx["n_neighbors"],
        "spatial_key": preprocess_ctx["spatial_key"],
        # Training params
        "n_epochs": params.n_epochs,
        "n_epochs_all_gps": params.n_epochs_all_gps,
        "lr": params.lr,
        "lambda_edge_recon": params.lambda_edge_recon,
        "lambda_gene_expr_recon": params.lambda_gene_expr_recon,
        "lambda_l1_masked": params.lambda_l1_masked,
        "lambda_l1_addon": params.lambda_l1_addon,
        "edge_batch_size": params.edge_batch_size,
        "n_sampled_neighbors": params.n_sampled_neighbors,
        "use_cuda_if_available": params.use_cuda_if_available,
        "conv_layer_encoder": params.conv_layer_encoder,
        "active_gp_thresh_ratio": params.active_gp_thresh_ratio,
        # Relative output paths (for reference)
        "run_root": params.run_root.name,
        "model_folder_path": str(last_n_levels(params.model_folder_path, 3)),
        "figure_folder_path": str(last_n_levels(params.figure_folder_path, 3)),
        "nichecompass_data_dir": f"{params.model_dir}/data",
        "outdir": "./",
    }

    with open(cfg_path, "w", encoding="utf-8") as f:
        json.dump(run_config, f, indent=2, default=_json_default)
    logging.info(f"Saved run_config.json at {cfg_path}")


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

    params.finalize_paths()
    params.run_root.mkdir(parents=True, exist_ok=True)
    params.figure_folder_path.mkdir(parents=True, exist_ok=True)
    params.model_folder_path.mkdir(parents=True, exist_ok=True)

    setup_logging(params.run_root / "train.log", params.debug)
    logging.info("=== NicheCompass Train: START ===")
    logging.info("Parameters: %s", json.dumps(asdict(params), indent=2, default=str))

    # Load preprocessed AnnData
    logging.info(f"Loading preprocessed AnnData from {params.preprocessed_h5ad}")
    adata = sc.read_h5ad(params.preprocessed_h5ad)

    preprocess_ctx: dict[str, Any] = adata.uns.get("nichecompass_preprocess_params", {})
    if not preprocess_ctx:
        raise KeyError(
            "adata.uns['nichecompass_preprocess_params'] is missing. "
            "Ensure the H5AD was produced by nichecompass_preprocess.py."
        )

    counts_key_effective: str = preprocess_ctx["counts_key_effective"]
    logging.info(f"Effective counts key: {counts_key_effective!r}")

    # Build, train, embed
    model = build_model(adata, params, preprocess_ctx, counts_key_effective)

    cat_keys: list[str] = [
        str(k) for k in (
            params.cat_covariates_keys
            or preprocess_ctx.get("cat_covariates_keys")
            or [preprocess_ctx["sample_key"]]
        )
    ]

    train_and_embed(model, params, latent_key=preprocess_ctx["latent_key"])

    save_model_and_config(
        model=model,
        params=params,
        preprocess_ctx=preprocess_ctx,
        counts_key_effective=counts_key_effective,
        cat_covariates_keys=cat_keys,
    )

    logging.info("=== NicheCompass Train: DONE ===")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
