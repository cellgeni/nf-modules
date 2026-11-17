#!/usr/bin/env python3
"""Transfer annotations from a single-cell reference onto spatial data with TACCO.

Runs ``tacco.tools.annotate`` (https://github.com/simonwm/tacco) and writes:

* ``<prefix>.h5ad``: the spatial AnnData with the annotation fractions in
  ``.obsm[<result_key>]``, the top annotation in ``.obs[<result_key>]`` and its
  fraction in ``.obs[<result_key>_fraction]``.
* ``<prefix>_annotation.csv.gz``: the same per-observation table (fractions per
  category, top annotation and its fraction).

Both AnnData objects must hold raw counts, in ``.X`` or in the layer given by
``--spatial-layer`` / ``--reference-layer``.
"""

import argparse
import json
import logging
import sys
from importlib.metadata import version

import anndata as ad
import pandas as pd
import tacco as tc

log = logging.getLogger("tacco_annotate")


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    io = p.add_argument_group("input/output")
    io.add_argument("--spatial", required=True, help="Spatial AnnData (.h5ad) to annotate.")
    io.add_argument("--reference", required=True, help="Annotated reference AnnData (.h5ad).")
    io.add_argument("--prefix", required=True, help="Prefix of the output files.")
    io.add_argument(
        "--spatial-layer",
        default=None,
        help="Layer of --spatial holding raw counts (default: .X).",
    )
    io.add_argument(
        "--reference-layer",
        default=None,
        help="Layer of --reference holding raw counts (default: .X).",
    )

    tl = p.add_argument_group("tacco.tl.annotate")
    tl.add_argument(
        "--annotation-key",
        default=None,
        help="Column of reference .obs with the annotation to transfer, e.g. cell_type "
        "(default: inferred by TACCO, which only works if it is the only .obs column).",
    )
    tl.add_argument(
        "--result-key",
        default=None,
        help="Key for the result in the spatial .obsm/.obs (default: --annotation-key).",
    )
    tl.add_argument(
        "--method",
        default="OT",
        help="Annotation method: OT, nnls, projection, svm, ... (default: OT).",
    )
    tl.add_argument("--bisections", type=int, default=None)
    tl.add_argument("--bisection-divisor", type=int, default=3)
    tl.add_argument("--platform-iterations", type=int, default=None)
    tl.add_argument("--normalize-to", choices=["adata", "reference"], default="adata")
    tl.add_argument("--multi-center", type=int, default=None)
    tl.add_argument("--max-annotation", type=int, default=None)
    tl.add_argument("--min-counts-per-gene", type=int, default=None)
    tl.add_argument("--min-counts-per-cell", type=int, default=None)
    tl.add_argument("--min-cells-per-gene", type=int, default=None)
    tl.add_argument("--min-genes-per-cell", type=int, default=None)
    tl.add_argument("--min-log2foldchange", type=float, default=None)
    tl.add_argument("--min-expression", type=float, default=None)
    tl.add_argument("--n-hvg", type=int, default=None)
    tl.add_argument("--remove-mito", action="store_true")
    tl.add_argument(
        "--assume-valid-counts",
        action="store_true",
        help="Skip TACCO's check for non-integer or negative counts.",
    )
    tl.add_argument(
        "--method-kwargs",
        type=json.loads,
        default={},
        help='JSON object of extra keyword arguments for the method, e.g. \'{"epsilon": 0.01}\'.',
    )
    tl.add_argument("--verbose", type=int, default=1)
    return p.parse_args(argv)


def use_counts_layer(adata, layer, name):
    if layer is None:
        return adata
    if layer not in adata.layers:
        sys.exit(f"ERROR: layer '{layer}' not found in {name} (layers: {list(adata.layers)})")
    adata.X = adata.layers[layer]
    return adata


def main(argv=None):
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    spatial = use_counts_layer(ad.read_h5ad(args.spatial), args.spatial_layer, "--spatial")
    reference = use_counts_layer(ad.read_h5ad(args.reference), args.reference_layer, "--reference")
    log.info("spatial: %d obs x %d vars", spatial.n_obs, spatial.n_vars)
    log.info("reference: %d obs x %d vars", reference.n_obs, reference.n_vars)

    if args.annotation_key is not None and args.annotation_key not in reference.obs:
        sys.exit(
            f"ERROR: --annotation-key '{args.annotation_key}' not in reference .obs "
            f"(columns: {list(reference.obs.columns)})"
        )
    shared = spatial.var_names.intersection(reference.var_names)
    log.info("%d genes shared between spatial and reference", len(shared))
    if len(shared) == 0:
        sys.exit("ERROR: spatial and reference share no var_names (genes)")

    fractions = tc.tl.annotate(
        spatial,
        reference,
        annotation_key=args.annotation_key,
        result_key=None,
        method=args.method,
        bisections=args.bisections,
        bisection_divisor=args.bisection_divisor,
        platform_iterations=args.platform_iterations,
        normalize_to=args.normalize_to,
        multi_center=args.multi_center,
        max_annotation=args.max_annotation,
        min_counts_per_gene=args.min_counts_per_gene,
        min_counts_per_cell=args.min_counts_per_cell,
        min_cells_per_gene=args.min_cells_per_gene,
        min_genes_per_cell=args.min_genes_per_cell,
        min_log2foldchange=args.min_log2foldchange,
        min_expression=args.min_expression,
        n_hvg=args.n_hvg,
        remove_mito=args.remove_mito,
        assume_valid_counts=args.assume_valid_counts,
        verbose=args.verbose,
        **args.method_kwargs,
    )

    # TACCO may drop observations (e.g. without counts): keep all, with NaN for those.
    fractions = fractions.reindex(spatial.obs_names)
    fractions.columns = fractions.columns.astype(str)
    dropped = int(fractions.isna().all(axis=1).sum())
    if dropped:
        log.warning("%d observations were not annotated by TACCO (NaN)", dropped)

    result_key = args.result_key or args.annotation_key or "tacco"
    annotated = fractions.notna().any(axis=1)
    top = pd.Series(pd.NA, index=fractions.index, dtype=object)
    top[annotated] = fractions[annotated].idxmax(axis=1)
    top = pd.Categorical(top, categories=list(fractions.columns))

    spatial.obsm[result_key] = fractions
    spatial.obs[result_key] = top
    spatial.obs[f"{result_key}_fraction"] = fractions.max(axis=1)
    spatial.uns["tacco"] = {
        "method": args.method,
        "annotation_key": str(args.annotation_key),
        "result_key": result_key,
        "n_shared_genes": len(shared),
        "tacco_version": version("tacco"),
    }
    spatial.write_h5ad(f"{args.prefix}.h5ad", compression="gzip")

    table = fractions.copy()
    table[result_key] = top
    table[f"{result_key}_fraction"] = spatial.obs[f"{result_key}_fraction"]
    table.index.name = "obs_name"
    table.to_csv(
        f"{args.prefix}_annotation.csv.gz",
        compression={"method": "gzip", "mtime": 0},
        float_format="%.6g",
    )
    log.info("annotation counts:\n%s", top.value_counts().to_string())


if __name__ == "__main__":
    main()
