# dotpy

Nextflow DSL2 module that runs [DOTpy](https://github.com/earmingol/DOTpy) — Deconvolution by Optimal
Transport — to transfer cell-type annotations from a single-cell reference (h5ad) onto spatial
transcriptomics data (h5ad).

## Contents

```
main.nf                        DOTPY process definition (runs run_dot_cli.py)
meta.yml                       module input/output/ext.args documentation
resources/usr/bin/             run_dot_cli.py — CLI wrapper around the dotpy Python package,
                                staged onto PATH via Nextflow module binaries
tests/                         nf-test suite and small synthetic h5ad fixtures
```

## Processes

| Process | Label                                                         | Container                       |
| ------- | ------------------------------------------------------------- | ------------------------------- |
| `DOTPY` | `process_gpu` if `params.dotpy_use_gpu` else `process_normal` | `quay.io/cellgeni/dotpy:latest` |

Set `--dotpy_use_gpu true` to schedule the job onto a GPU queue; the script itself
auto-detects CUDA (`--device auto`) regardless of which queue it lands on.

## Input

| Name       | Type | Description                                                           |
| ---------- | ---- | --------------------------------------------------------------------- |
| `meta`     | map  | Sample metadata. Must contain at least `id`, e.g. `[ id: 'sample1' ]` |
| `sc_adata` | file | Single-cell reference AnnData file (`.h5ad`)                          |
| `sp_adata` | file | Spatial transcriptomics AnnData file (`.h5ad`)                        |

## Output

| Emit             | Pattern             | Description                                                               |
| ---------------- | ------------------- | ------------------------------------------------------------------------- |
| `weights`        | `*_weights.csv`     | Per-spot cell-type deconvolution weights (one column per cell type)       |
| `annotations`    | `*_annotations.csv` | Per-spot winning cell-type assignment (argmax of weights)                 |
| `versions_dotpy` | topic: `versions`   | `dotpy` package version, emitted via the nf-core topic-channel convention |

See [`meta.yml`](meta.yml) for full descriptions.

### `ext.args`

All `run_dot_cli.py` flags are passed via `task.ext.args` in the including pipeline's
`nextflow.config`. The module always sets `--no-h5ad --no-plots` (CSV output only).

Example `ext.args` per platform:

- **Visium** (`lowres` mode; spots cover multiple cells, so abundance-ratio matching helps):

  ```groovy
  ext.args = '--cell-type-key cell_type --mode lowres --ratios-weight 0.3'
  ```

- **Xenium** (`highres` mode; single-cell resolution, ratio matching disabled):

  ```groovy
  ext.args = '--cell-type-key cell_type --mode highres'
  ```

### `--ext_args` flags

| Flag                 | Default           | Description                                                                                                                   |
| -------------------- | ----------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| `--no-plots`         | `off`             | Skip plot generation                                                                                                          |
| `--no-h5ad`          | `off`             | Skip saving per-sample h5ad files                                                                                             |
| `--cell-type-key`    | `cell_type`       | Column in `sc_adata.obs` containing cell-type labels                                                                          |
| `--mode`             | `highres`         | Resolution mode: `highres` (Xenium/MERFISH/CosMx) or `lowres` (Visium/ST)                                                     |
| `--counts-layer`     | uses `sp_adata.X` | Layer in `sp_adata` with raw counts — set this if `sp_adata.X` is normalized. Use `X` to explicitly force use of `sp_adata.X` |
| `--ref-counts-layer` | uses `sc_adata.X` | Layer in `sc_adata` with raw counts — set this if `sc_adata.X` is normalized                                                  |
| `--device`           | `auto`            | Compute device: `auto`, `cuda`, or `cpu`                                                                                      |
| `--iterations`       | `100`             | Number of Frank-Wolfe optimisation iterations                                                                                 |
| `--batch-size`       | `5000`            | GPU batch size (reduce if CUDA OOM)                                                                                           |
| `--ratios-weight`    | `0.0`             | Weight for cell-type abundance matching (0 disables; ~0.3 for lowres)                                                         |
| `--subcluster-size`  | `10`              | Maximum subclusters per cell type in the reference                                                                            |
| `--max-genes`        | `5000`            | Maximum number of genes used during optimisation                                                                              |
| `--th-spatial`       | `0.84`            | Threshold on spatial similarity of adjacent spots                                                                             |
| `--th-gene-low`      | `0.01`            | Minimum fraction of spots a spatial gene must be expressed in (0 disables)                                                    |
| `--th-gene-high`     | `0.99`            | Maximum fraction of spots a spatial gene must be expressed in (1 disables)                                                    |
| `--mixed-precision`  | `false`           | Boolean flag (no value); enables float16 on GPU                                                                               |
| `--sample-key`       | `''`              | Column in `sp_adata.obs` for multi-sample splitting within one h5ad (omit for single)                                         |
| `--lineage-key`      | `''`              | Column in `sc_adata.obs` with lineage annotations (optional)                                                                  |
| `--save-combined`    | `false`           | Boolean flag (no value); saves a combined results file at the end                                                             |
| `--checkpoint-dir`   | unset (disabled)  | Directory for saving optimisation checkpoints                                                                                 |
| `--checkpoint-freq`  | `10`              | Save a checkpoint every N iterations (only used if `--checkpoint-dir` is set)                                                 |
| `--resume-from`      | unset             | Path to a checkpoint file to resume from                                                                                      |
| `--verbose`          | `false`           | Boolean flag (no value); prints detailed progress messages                                                                    |

See [`meta.yml`](meta.yml) for full descriptions.
