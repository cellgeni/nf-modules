#!/usr/bin/env python3
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import sparse

import scanpy as sc
import tangram as tg
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


def _read_markers(markers_path: str | None) -> list[str] | None:
    """
    Read marker genes from csv.
    Expected format: either a 1-column csv, or any matrix where values are gene names.
    Returns a de-duplicated list of strings, preserving order as best as possible.
    """
    if not markers_path:
        return None
    if not os.path.exists(markers_path):
        raise FileNotFoundError(f"Markers file not found: {markers_path}")

    df = pd.read_csv(markers_path, index_col=0)
    vals = df.values.ravel()
    vals = [str(x) for x in vals if pd.notnull(x)]
    # de-duplicate preserving order
    seen = set()
    out = []
    for g in vals:
        if g not in seen:
            out.append(g)
            seen.add(g)
    return out or None


def _pp_adatas_compat(ad_sc, ad_sp, genes=None):
    """
    Tangram has moved preprocessing helper across versions.
    """
    try:
        preprocess = tg.pp_adatas  # older tangram
    except AttributeError:
        preprocess = tg.utils.pp_adatas  # tangram-sc >= 1.0.4
    preprocess(ad_sc, ad_sp, genes=genes)


def _normalize_rows(mat):
    """
    Row-normalize a (dense or sparse) matrix so each row sums to 1 (if sum>0).
    This is useful for interpreting proportions per spatial cell/location.
    """
    if sparse.issparse(mat):
        row_sums = np.asarray(mat.sum(axis=1)).ravel()
        inv = np.zeros_like(row_sums, dtype=float)
        nz = row_sums > 0
        inv[nz] = 1.0 / row_sums[nz]
        D = sparse.diags(inv)
        return D @ mat
    else:
        row_sums = mat.sum(axis=1)
        inv = np.zeros_like(row_sums, dtype=float)
        nz = row_sums > 0
        inv[nz] = 1.0 / row_sums[nz]
        return (mat.T * inv).T


def _subsample_adata(
    ad_sc, max_cells: int, cell_type_label: str | None = None, seed: int = 42
):
    """
    Subsample AnnData to max_cells, optionally preserving cell type distribution.

    Args:
        ad_sc: AnnData object to subsample
        max_cells: Maximum number of cells to retain
        cell_type_label: Optional column name for stratified sampling
        seed: Random seed for reproducibility

    Returns:
        Subsampled AnnData object
    """
    n_cells = ad_sc.n_obs

    if n_cells <= max_cells:
        print(
            f"[Tangram] scRNA-seq has {n_cells:,} cells (under {max_cells:,} threshold, no subsampling needed)"
        )
        return ad_sc

    print(
        f"[Tangram] WARNING: scRNA-seq has {n_cells:,} cells, which may cause memory issues."
    )
    print(f"[Tangram] Subsampling to {max_cells:,} representative cells...")

    np.random.seed(seed)

    # Stratified sampling if cell type label is available
    if cell_type_label and cell_type_label in ad_sc.obs.columns:
        cell_types = ad_sc.obs[cell_type_label]
        sample_indices = []

        for ct in cell_types.unique():
            ct_indices = np.where(cell_types == ct)[0]
            n_ct = len(ct_indices)
            # Sample proportionally, but at least 10 cells per type (or all if less than 10)
            n_sample = max(min(10, n_ct), int(max_cells * n_ct / n_cells))
            n_sample = min(n_sample, n_ct)  # Don't sample more than available

            sampled = np.random.choice(ct_indices, size=n_sample, replace=False)
            sample_indices.extend(sampled)

        # If we sampled too many, trim to max_cells
        if len(sample_indices) > max_cells:
            sample_indices = np.random.choice(
                sample_indices, size=max_cells, replace=False
            )

        ad_sc_sub = ad_sc[sample_indices, :].copy()
        print(
            f"[Tangram] Stratified subsampling complete: {ad_sc_sub.n_obs:,} cells retained"
        )
        print(
            f"[Tangram] Cell type distribution preserved across {len(cell_types.unique())} cell types"
        )
    else:
        # Random sampling
        sample_indices = np.random.choice(n_cells, size=max_cells, replace=False)
        ad_sc_sub = ad_sc[sample_indices, :].copy()
        print(
            f"[Tangram] Random subsampling complete: {ad_sc_sub.n_obs:,} cells retained"
        )

    return ad_sc_sub


def _plot_training_scores(adata_map, bins=10, alpha=0.7, save_path=None):
    """
    Plots the 4-panel training diagnosis plot. Adapted from Tangram's plot_utils.py.
    """
    if "train_genes_df" not in adata_map.uns:
        raise ValueError("No training info found in adata_map.uns['train_genes_df'].")

    df = adata_map.uns["train_genes_df"]

    fig, axs = plt.subplots(1, 4, figsize=(14, 3), sharey=True)
    axs_f = axs.flatten()

    axs_f[0].set_ylim([0.0, 1.0])
    for i in range(1, len(axs_f)):
        axs_f[i].set_xlim([0.0, 1.0])
        axs_f[i].set_ylim([0.0, 1.0])

    sns.histplot(data=df, y="train_score", bins=bins, ax=axs_f[0], color="coral")
    axs_f[0].set_title("Distribution of training scores")

    sns.scatterplot(
        data=df,
        y="train_score",
        x="sparsity_sc",
        ax=axs_f[1],
        alpha=alpha,
        color="coral",
    )
    axs_f[1].set_title("score vs sparsity (single cells)")

    sns.scatterplot(
        data=df,
        y="train_score",
        x="sparsity_sp",
        ax=axs_f[2],
        alpha=alpha,
        color="coral",
    )
    axs_f[2].set_title("score vs sparsity (spatial)")

    sns.scatterplot(
        data=df,
        y="train_score",
        x="sparsity_diff",
        ax=axs_f[3],
        alpha=alpha,
        color="coral",
    )
    axs_f[3].set_title("score vs sparsity (sp - sc)")

    plt.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[Tangram] Saved training scores plot to: {save_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Tangram wrapper: map scRNA-seq data to spatial transcriptomics and derive spatial cell-type calls."
    )
    parser.add_argument("-s", "--scrna", required=True, help="Path to scRNA .h5ad")
    parser.add_argument("-p", "--spatial", required=True, help="Path to spatial .h5ad")
    parser.add_argument(
        "-m",
        "--markers",
        required=False,
        default=None,
        help="Optional markers CSV (genes to use for training).",
    )
    parser.add_argument(
        "-t",
        "--tech",
        default="h5ad",
        help="Spatial tech (only 'h5ad' supported in this script).",
    )
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    parser.add_argument(
        "--mode",
        default="cells",
        choices=["cells", "clusters"],
        help="Mapping mode (default: cells)",
    )
    parser.add_argument(
        "-c",
        "--cluster-label",
        default=None,
        help="Obs column with cluster IDs (required if --mode clusters)",
    )
    parser.add_argument(
        "-d",
        "--device",
        default="cpu",
        help="Compute device: 'cpu', 'cuda', 'cuda:0', etc. Default: cpu",
    )

    # Cell type label
    parser.add_argument(
        "--cell-type-label",
        required=True,
        help="Obs column in scRNA AnnData containing cell types (e.g. 'class_label').",
    )

    # Training parameters
    parser.add_argument(
        "--num-epochs",
        type=int,
        default=1000,
        help="Number of training epochs (default: 1000)",
    )
    parser.add_argument(
        "--learning-rate",
        type=float,
        default=0.1,
        help="Learning rate for training (default: 0.1)",
    )

    # Subsampling parameters
    parser.add_argument(
        "--subsample",
        action="store_true",
        help="Enable subsampling of scRNA-seq data to reduce memory usage",
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=10000,
        help="Maximum number of cells to use from scRNA-seq when subsampling (default: 10000)",
    )
    parser.add_argument(
        "--subsample-seed",
        type=int,
        default=42,
        help="Random seed for subsampling reproducibility (default: 42)",
    )

    # Outputs
    parser.add_argument(
        "--save-props-csv",
        action="store_true",
        help="Also save tangram cell-type proportions as CSV.",
    )
    parser.add_argument(
        "--ct-props-key",
        default="tangram_ct_props",
        help="Key for storing cell-type proportions in adata_sp.obsm (default: tangram_ct_props).",
    )
    parser.add_argument(
        "--ct-max-key",
        default="tangram_ct_max",
        help="Column name for hard label in adata_sp.obs (default: tangram_ct_max).",
    )
    parser.add_argument(
        "--ct-max-score-key",
        default="tangram_ct_max_score",
        help="Column name for the max proportion score in adata_sp.obs.",
    )

    # Training diagnostics plot
    parser.add_argument(
        "--plot-bins",
        type=int,
        default=10,
        help="Number of bins for the training score histogram (default=10).",
    )
    parser.add_argument(
        "--plot-alpha",
        type=float,
        default=0.7,
        help="Opacity for the training score scatterplots (default=0.7).",
    )

    args = parser.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    # Guard against gzipped inputs
    for p, label in [(args.scrna, "scRNA"), (args.spatial, "spatial")]:
        if p.endswith(".h5ad.gz"):
            sys.exit(
                f"[Tangram] ERROR: {label} file '{p}' is gzipped (.h5ad.gz). "
                f"Please gunzip to a plain .h5ad first."
            )

    print(f"[Tangram] Reading scRNA:   {args.scrna}")
    ad_sc = sc.read_h5ad(args.scrna)

    print(f"[Tangram] Reading spatial: {args.spatial}")
    if args.tech.lower() != "h5ad":
        print(
            "[Tangram] Warning: only 'h5ad' loading is supported here; reading as h5ad."
        )
    ad_sp = sc.read_h5ad(args.spatial)

    if args.cell_type_label not in ad_sc.obs.columns:
        avail = ", ".join(ad_sc.obs.columns[:30]) + (
            " ..." if ad_sc.obs.shape[1] > 30 else ""
        )
        sys.exit(
            f"[Tangram] ERROR: --cell-type-label '{args.cell_type_label}' not found in ad_sc.obs. "
            f"Available columns (truncated): {avail}"
        )

    markers = _read_markers(args.markers)
    if markers:
        print(f"[Tangram] Loaded {len(markers)} marker genes from {args.markers}")
    else:
        print(
            "[Tangram] No markers provided; Tangram will use intersecting genes (default preprocessing)."
        )

    print("[Tangram] Preprocessing (harmonize gene space, normalize/log)...")
    _pp_adatas_compat(ad_sc, ad_sp, genes=markers)

    # Apply subsampling if requested
    if args.subsample:
        ad_sc = _subsample_adata(
            ad_sc,
            max_cells=args.max_cells,
            cell_type_label=args.cell_type_label,
            seed=args.subsample_seed,
        )

        # ===== FIX: Filter genes with zero expression after subsampling =====
        print("[Tangram] Filtering genes with zero expression after subsampling...")

        # Count genes before filtering
        n_genes_before = ad_sc.n_vars

        # Filter genes with at least 1 cell expressing in scRNA
        sc.pp.filter_genes(ad_sc, min_cells=1)

        # Filter genes with at least 1 cell expressing in spatial
        sc.pp.filter_genes(ad_sp, min_cells=1)

        # Re-harmonize gene space after filtering
        common_genes = sorted(list(set(ad_sc.var_names) & set(ad_sp.var_names)))
        print(
            f"[Tangram] Gene filtering: {n_genes_before} → {len(common_genes)} genes (removed {n_genes_before - len(common_genes)} zero-expression genes)"
        )

        ad_sc = ad_sc[:, common_genes].copy()
        ad_sp = ad_sp[:, common_genes].copy()

        # CRITICAL: Re-run Tangram preprocessing after filtering to update internal gene lists
        print("[Tangram] Re-running preprocessing to synchronize gene spaces...")
        _pp_adatas_compat(ad_sc, ad_sp, genes=markers)
        # ===== END FIX =====
    else:
        print(
            f"[Tangram] Subsampling disabled. Using all {ad_sc.n_obs:,} cells from scRNA-seq."
        )

    # Device sanity (Tangram uses torch under the hood)
    dev = args.device.strip()
    dev_l = dev.lower()
    if dev_l != "cpu" and not dev_l.startswith("cuda"):
        print(
            f"[Tangram] Warning: device '{args.device}' not recognized; falling back to CPU."
        )
        dev = "cpu"

    # Clusters mode checks
    if args.mode == "clusters":
        if not args.cluster_label:
            sys.exit(
                "[Tangram] ERROR: --cluster-label is required when --mode clusters"
            )
        if args.cluster_label not in ad_sc.obs.columns:
            avail = ", ".join(ad_sc.obs.columns[:30]) + (
                " ..." if ad_sc.obs.shape[1] > 30 else ""
            )
            sys.exit(
                f"[Tangram] ERROR: --cluster-label '{args.cluster_label}' not found in ad_sc.obs. "
                f"Available columns (truncated): {avail}"
            )

    print(f"[Tangram] Running mapping (mode={args.mode}, device={dev}) ...")
    print(
        f"[Tangram] Training parameters: num_epochs={args.num_epochs}, learning_rate={args.learning_rate}"
    )

    ad_map = tg.map_cells_to_space(
        adata_sc=ad_sc,
        adata_sp=ad_sp,
        mode=args.mode,
        cluster_label=args.cluster_label,
        device=dev,
        num_epochs=args.num_epochs,
        learning_rate=args.learning_rate,
        density_prior="rna_count_based",
    )

    print("[Tangram] Plotting training diagnostics ...")
    try:
        figures_dir = os.path.join(args.outdir, "figures")
        os.makedirs(figures_dir, exist_ok=True)
        _plot_training_scores(
            ad_map,
            bins=args.plot_bins,
            alpha=args.plot_alpha,
            save_path=os.path.join(figures_dir, "training_scores.png"),
        )
    except ValueError as e:
        print(f"[Tangram] Skipping training diagnostics plot: {e}")

    print("[Tangram] Projecting genes ...")
    # NOTE: for gene projection, cluster_label is only relevant if you project from clusters;
    # here we keep your previous behavior, but do not force it.
    ad_ge = tg.project_genes(
        adata_map=ad_map,
        adata_sc=ad_sc,
        cluster_label=None,  # keep it simple; mapping is cell-level unless user uses cluster mode + special handling
    )

    # -------------------------
    # NEW: derive spatial cell-type proportions + hard labels
    # -------------------------
    print("[Tangram] Deriving spatial cell-type proportions from mapping matrix ...")

    # mapping: cells(sc) x spatial
    M = ad_map.X
    if sparse.issparse(M):
        M = M.tocsr()

    # one-hot encode sc cell types (cells x types)
    ct = ad_sc.obs[args.cell_type_label].astype("category")
    ct_categories = list(ct.cat.categories)
    onehot = pd.get_dummies(ct, dtype=float)  # dense (cells x types)
    onehot_mat = onehot.values

    # compute (spatial x types) = (spatial x cells) @ (cells x types)
    # Note: ad_map is (cells x spatial), so transpose first
    if sparse.issparse(M):
        props = M.T @ onehot_mat  # yields dense np.ndarray (spatial x types)
    else:
        props = M.T.dot(onehot_mat)

    # normalize each spatial row to sum to 1 => proportions
    props = _normalize_rows(props)

    # store in spatial AnnData with obs exactly preserved
    ad_sp_ct = ad_sp.copy()
    ad_sp_ct.obsm[args.ct_props_key] = props
    ad_sp_ct.uns[f"{args.ct_props_key}_columns"] = ct_categories

    # hard label + score
    if sparse.issparse(props):
        props_dense = props.toarray()
    else:
        props_dense = np.asarray(props)

    max_idx = props_dense.argmax(axis=1)
    max_score = props_dense[np.arange(props_dense.shape[0]), max_idx]
    pred_ct = [ct_categories[i] for i in max_idx]

    ad_sp_ct.obs[args.ct_max_key] = pd.Categorical(pred_ct, categories=ct_categories)
    ad_sp_ct.obs[args.ct_max_score_key] = max_score

    # -------------------------
    # Save outputs
    # -------------------------
    out_map_h5ad = os.path.join(args.outdir, "tangram_aligned.h5ad")
    out_ge_h5ad = os.path.join(args.outdir, "tangram_gene_proj.h5ad")
    out_sp_ct_h5ad = os.path.join(args.outdir, "spatial_with_tangram_celltypes.h5ad")
    out_log = os.path.join(args.outdir, "tangram.log")

    ad_map.write_h5ad(out_map_h5ad)
    print(f"[Tangram] Saved mapping matrix → {out_map_h5ad}")

    ad_ge.write_h5ad(out_ge_h5ad)
    print(f"[Tangram] Saved gene projection → {out_ge_h5ad}")

    ad_sp_ct.write_h5ad(out_sp_ct_h5ad)
    print(f"[Tangram] Saved spatial w/ cell types → {out_sp_ct_h5ad}")

    if args.save_props_csv:
        csv_path = os.path.join(args.outdir, "tangram_celltype_props.csv")
        pd.DataFrame(
            props_dense, index=ad_sp_ct.obs_names, columns=ct_categories
        ).to_csv(csv_path)
        print(f"[Tangram] Saved cell-type proportions CSV → {csv_path}")

    with open(out_log, "w") as f:
        f.write(f"scRNA: {args.scrna}\nSpatial: {args.spatial}\n")
        f.write(f"Mode: {args.mode}\nDevice: {dev}\n")
        f.write(
            f"Training: num_epochs={args.num_epochs}, learning_rate={args.learning_rate}\n"
        )
        f.write(f"Cell type label (scRNA obs): {args.cell_type_label}\n")
        if args.subsample:
            f.write(
                f"Subsampling: enabled (max_cells={args.max_cells}, seed={args.subsample_seed})\n"
            )
        else:
            f.write(f"Subsampling: disabled\n")
        if args.mode == "clusters":
            f.write(f"Cluster label: {args.cluster_label}\n")
        f.write(f"Markers: {args.markers if args.markers else 'None'}\n")
        f.write(
            f"Saved: {out_map_h5ad}\nSaved: {out_ge_h5ad}\nSaved: {out_sp_ct_h5ad}\n"
        )
        f.write("Done.\n")
    print(f"[Tangram] Log written → {out_log}")

    print("[Tangram] Done.")


if __name__ == "__main__":
    main()
