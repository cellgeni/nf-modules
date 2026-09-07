# baysor

Nextflow DSL2 modules wrapping [Baysor](https://github.com/kharchenkolab/Baysor) — Bayesian cell
segmentation of imaging-based spatial transcriptomics data (Xenium, MERFISH, CosMx, MERSCOPE).

Baysor segments cell bodies directly from the molecule table using a Bayesian mixture model over
transcript positions and expression, optionally guided by a prior (e.g. nuclear) segmentation.

## Contents

```
run/main.nf     BAYSOR_RUN process definition (runs `baysor run`)
run/meta.yml    module input/output documentation
run/tests/      nf-test suite and small synthetic molecule-table fixtures
```

Upstream also ships `baysor preview` and `baysor segfree` subcommands; only `run` is wrapped so far.
Add them as `preview/` and `segfree/` alongside `run/` if needed.

## Container

`quay.io/cellgeni/baysor:cpp-0.8.2` — the **C++ port** of Baysor
([`cpp-0.8.2`](https://github.com/kharchenkolab/Baysor/tree/cpp-0.8.2) upstream), not the original
Julia implementation. This matters: the CLI and the output file names differ from Julia Baysor
≤ 0.7.1. See [Migrating from the Julia-era module](#migrating-from-the-julia-era-module).

The `baysor` binary exposes **no version flag**, so `BAYSOR_RUN` reports a hardcoded `0.8.2` on the
`versions` topic. Update that literal in `run/main.nf` whenever the container tag is bumped.

## Threading

Baysor is OpenMP-parallel but has no `--threads` option, so the module exports
`OMP_NUM_THREADS=${task.cpus}`. Without it Baysor would grab every core on the host regardless of
what the scheduler allocated. Measured on a 120k-molecule synthetic dataset, 200 iterations:

| `OMP_NUM_THREADS` | wall  | cpu   |
| ----------------- | ----- | ----- |
| 1                 | 39.9s | 39.8s |
| 8                 | 24.3s | 70.3s |

Set `cpus` through the `process_high` label or a `withName: BAYSOR_RUN` selector.

## Input

| Position | Name                 | Type      | Description                                                                                                  |
| -------- | -------------------- | --------- | ------------------------------------------------------------------------------------------------------------ |
| `[0]`    | `meta`               | map       | Sample metadata; must contain at least `id`, e.g. `[ id: 'sample1' ]`                                        |
| `[0]`    | `coordinates`        | file      | Molecule table (CSV/Parquet) or a Xenium `experiment.xenium` manifest                                        |
| `[0]`    | `prior_segmentation` | file      | Optional prior segmentation — label image mask, or cell boundary CSV/Parquet. Pass `[]` for none              |
| `[1]`    | `config`             | file      | Optional TOML config (`--config`). Pass `[]` to configure entirely through `ext.args`                        |

Column names are **not** module inputs — pass them through `ext.args` (`-x`, `-y`, `-z`, `-g`,
`--qv-column`). Defaults are `x`, `y`, `z`, `gene`, `qv`, which rarely match a vendor export.

### Prior segmentation

Baysor accepts the prior as a second positional argument in three forms. Two are files and go
through the `prior_segmentation` input; the third is a reference to a column already present in the
molecule table, which is not a file and therefore goes through `ext.prior` instead:

```groovy
process {
    withName: BAYSOR_RUN {
        // use the `nucleus_id` column of the molecule table as the prior
        ext.prior = ':nucleus_id'
        ext.args  = '-x x_location -y y_location -g feature_name --prior-segmentation-confidence 0.5'
    }
}
```

`prior_segmentation` takes precedence over `ext.prior` when both are set.

## Output

Output file names depend on `--output-style` (`legacy`, the default, or `parquet`). Each emit
covers both styles under one semantic name:

| Emit              | `legacy`                          | `parquet`                   | Notes                                     |
| ----------------- | --------------------------------- | --------------------------- | ----------------------------------------- |
| `outdir`          | `*_baysor/`                       | `*_baysor/`                 | The whole bundle, as a directory          |
| `molecules`       | `segmentation.csv`                | `molecules.parquet`         | Molecules + assigned cell id / confidence |
| `cell_stats`      | `segmentation_cell_stats.csv`     | `cells.parquet`             | Per-cell centroid, area, counts           |
| `counts`          | `segmentation_counts.{loom,tsv}`  | `feature_matrix.h5`         | Cell × gene matrix                        |
| `polygons_2d`     | `segmentation_polygons_2d.json`   | `cell_boundaries.parquet`   | Optional — absent with `--polygon-format none` |
| `polygons_3d`     | `segmentation_polygons_3d.json`   | `cell_boundaries_3d.parquet`| Optional — 3D input only                  |
| `plot`            | `segmentation_plot.html`          | `segmentation_plot.html`    | Optional — requires `--plot`              |
| `report`          | `diagnostic_report.html`          | `diagnostic_report.html`    | Optional — requires `--plot`              |
| `params`          | `segmentation_params.dump.toml`   | `run_params.toml`           | Fully resolved parameters, reusable as `config` |
| `log`             | `segmentation_log.log`            | `run.log`                   | Includes per-iteration convergence        |

Note that `--count-matrix-format` only applies to `legacy`; `parquet` always writes
`feature_matrix.h5`.

## Example

```groovy
include { BAYSOR_RUN } from '../modules/sanger-cellgeni/baysor/run/main.nf'

workflow {
    ch_transcripts = Channel.of([[id: 'xenium_sample'], file('transcripts.parquet'), []])
    BAYSOR_RUN(ch_transcripts, [])
}
```

```groovy
process {
    withName: BAYSOR_RUN {
        cpus      = 16
        memory    = 64.GB
        ext.args  = [
            '-x x_location -y y_location -z z_location -g feature_name',
            '--min-qv 20',
            '-m 30',
            '-s 8',
            '--n-clusters 6',
            '--plot',
        ].join(' ')
    }
}
```

## Key `ext.args` flags

Full list via `docker run --rm quay.io/cellgeni/baysor:cpp-0.8.2 baysor run --help`.

| Flag                              | Default        | Description                                                                       |
| --------------------------------- | -------------- | --------------------------------------------------------------------------------- |
| `-x`, `-y`, `-z`, `-g`            | `x`/`y`/`z`/`gene` | Coordinate and gene column names                                              |
| `--qv-column`                     | `qv`           | Quality-value column used by `--min-qv`                                           |
| `-s`, `--scale`                   | estimated      | Approximate cell radius. Setting it disables scale estimation from centers        |
| `--scale-std`                     | `25%`          | Std of cell radius across cells; a number or a percentage of `--scale`            |
| `-m`, `--min-molecules-per-cell`  | —              | Minimum molecules for a cell to be considered real                                |
| `--prior-segmentation-confidence` | `0.2`          | How strongly to trust the prior, in `[0, 1]`                                      |
| `--cluster-method`                | `mrf`          | Molecule clustering prior: `mrf`, `louvain`, `leiden`, `none`                     |
| `--n-clusters`                    | `4` (mrf)      | Target number of molecule clusters / major cell types                             |
| `--min-qv`                        | —              | Drop molecules below this quality value at load time                              |
| `--exclude-genes`                 | —              | Comma-separated genes or patterns to drop, e.g. `'Blank*,MALAT1'`                 |
| `--force-2d`                      | off            | Ignore the z column                                                               |
| `--iters`                         | `500`          | Maximum optimisation iterations                                                   |
| `--tol`                           | `0.005`        | Convergence tolerance; `0` always runs all `--iters`                              |
| `--output-style`                  | `legacy`       | `legacy` or `parquet` output bundle                                               |
| `--polygon-format`                | `FeatureCollection` | `FeatureCollection`, `GeometryCollection`, or `none`                         |
| `--count-matrix-format`           | `loom`         | `loom` or `tsv` (`legacy` style only)                                             |
| `-p`, `--plot`                    | off            | Write `segmentation_plot.html` and `diagnostic_report.html`                       |
| `--skip-ncv-color`                | off            | Skip the NCV colour embedding — a large speedup when plots are not needed         |

## Migrating from the Julia-era module

The earlier local module (`BAYSOR_RUN` on `khersameesh24/baysor:0.7.1`) wrapped Julia Baysor.
Differences to expect when swapping this module in:

- **Prior and config are no longer bare `path` channels.** `prior_segmentation` moved into the
  input tuple so it can vary per sample; `config` stays a separate channel but is now optional.
- **`scale` is no longer a `val` input.** Pass `-s`/`--scale` through `ext.args`.
- **`--plot` is no longer forced on.** Plotting is a substantial fraction of runtime, so request it
  explicitly via `ext.args` when you want `plot`/`report`.
- **Polygons keep their native `.json` extension.** The old module renamed them to `.geojson`; this
  one does not rewrite tool output. Rename downstream if a consumer requires that suffix.
- **Emit names changed** — `segmentation` → `molecules`, `polygons2d`/`polygons3d` →
  `polygons_2d`/`polygons_3d`, `stats` → `cell_stats`, `loom` → `counts`, `htmls` →
  `plot`/`report`. `polygons_3d` is now `optional: true`, since a 2D run does not produce it.
- **Versions come over the `versions` topic channel**, not a `versions.yml` file, matching the rest
  of this repo.
- **No Conda guard needed.** The old module errored out under `-profile conda`; Baysor is not
  packaged for Conda either way, and this repo's modules are container-only.

## Tests

```bash
nf-test test modules/sanger-cellgeni/baysor/run/tests/main.nf.test
```

The suite covers a 3D run, a `--force-2d` run with TSV counts, an `ext.prior` column prior, a TOML
`config` run, and a stub run, against small synthetic molecule tables committed under `run/tests/`.
