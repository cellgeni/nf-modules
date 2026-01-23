#!/usr/bin/env Rscript

library(DOTr)
library(ggplot2)
library(anndata)
library(reticulate)
library(Matrix)
library(tools)


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop(paste(
    "Usage: Rscript run_dot.R <adata_sc_path> <adata_sp_path> <output_path> [options]",
    "",
    "Options:",
    "  --x-col <name>                 Spatial x coordinate column (default: x_centroid)",
    "  --y-col <name>                 Spatial y coordinate column (default: y_centroid)",
    "  --th-gene-low <num>            Low gene expression threshold (default: 0)",
    "  --th-gene-high <num>           High gene expression threshold (default: 1)",
    "  --identities-col <name>        Column in sc obs for annotations (default: identities)",
    "  --ref-subcluster-size <int>    Reference subcluster size (default: 3)",
    "  --max-genes <int>              Max genes for reference (default: adata_sc$n_vars)",
    "  --verbose <true|false>         Verbose output (default: true)",
    "  --ratios-weight <num>          DOT ratios weight (default: 0)",
    "  --weights-dir <name>           Subfolder for weights output (default: DOT_weights)",
    "  --plot-dir <name>              Subfolder for plot output (default: DOT_plotdata)",
    "  --weights-prefix <name>        Prefix for weights CSV (default: DOT_weights_)",
    "  --plot-prefix <name>           Prefix for plot CSV (default: DOT_plotting_)",
    sep = "\n"
  ))
}

parse_optional_args <- function(args) {
  opts <- list()
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (!grepl("^--", arg)) {
      stop(paste("Unexpected argument:", arg))
    }
    key_val <- sub("^--", "", arg)
    if (grepl("=", key_val, fixed = TRUE)) {
      parts <- strsplit(key_val, "=", fixed = TRUE)[[1]]
      key <- parts[1]
      val <- parts[2]
    } else {
      key <- key_val
      if (i == length(args) || grepl("^--", args[i + 1])) {
        val <- "TRUE"
      } else {
        val <- args[i + 1]
        i <- i + 1
      }
    }
    opts[[key]] <- val
    i <- i + 1
  }
  opts
}

get_opt <- function(opts, name, default) {
  if (!is.null(opts[[name]])) {
    opts[[name]]
  } else {
    default
  }
}

as_bool <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

adata_sc_path <- args[1]
adata_sp_path <- args[2]
output_path <- args[3]
opt_args <- if (length(args) > 3) args[4:length(args)] else character(0)
opts <- parse_optional_args(opt_args)

x_col <- get_opt(opts, "x-col", "x_centroid")
y_col <- get_opt(opts, "y-col", "y_centroid")
th_gene_low <- as.numeric(get_opt(opts, "th-gene-low", 0))
th_gene_high <- as.numeric(get_opt(opts, "th-gene-high", 1))
identities_col <- get_opt(opts, "identities-col", "identities")
ref_subcluster_size <- as.integer(get_opt(opts, "ref-subcluster-size", 3))
max_genes_opt <- get_opt(opts, "max-genes", NA)
verbose <- as_bool(get_opt(opts, "verbose", TRUE))
ratios_weight <- as.numeric(get_opt(opts, "ratios-weight", 0))
weights_dir_name <- get_opt(opts, "weights-dir", "DOT_weights")
plot_dir_name <- get_opt(opts, "plot-dir", "DOT_plotdata")
weights_prefix <- get_opt(opts, "weights-prefix", "DOT_weights_")
plot_prefix <- get_opt(opts, "plot-prefix", "DOT_plotting_")

# Extract base name
sp_base <- file_path_sans_ext(basename(adata_sp_path))

# Define subfolders
weights_dir <- file.path(output_path, weights_dir_name)
plot_dir <- file.path(output_path, plot_dir_name)

# Create directories if they don't exist
dir.create(weights_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Build full output file paths
DOT_weights_path <- file.path(weights_dir, paste0(weights_prefix, sp_base, ".csv"))
DOT_plot_path <- file.path(plot_dir, paste0(plot_prefix, sp_base, ".csv"))

# Load expression data
adata_sc <- read_h5ad(adata_sc_path)
adata_sp <- read_h5ad(adata_sp_path)

# Setup DOT with correct expression data and spatial coordinates
expression <- adata_sp$X
if (!all(c(x_col, y_col) %in% colnames(adata_sp$obs))) {
  stop(paste("Spatial columns not found in adata_sp$obs:", x_col, y_col))
}
spatial <- adata_sp$obs[, c(x_col, y_col), drop = FALSE]
dot.srt <- setup.srt(t(expression), spatial, th.gene.low = th_gene_low, th.gene.high = th_gene_high)

# Select reference cell types from adata_sc
if (!(identities_col %in% colnames(adata_sc$obs))) {
  stop(paste("Identities column not found in adata_sc$obs:", identities_col))
}
annotation_list <- as.character(adata_sc$obs[[identities_col]])
matrix_sc <- as(t(adata_sc$X), "CsparseMatrix")
max_genes <- if (is.na(max_genes_opt)) as.integer(adata_sc$n_vars) else as.integer(max_genes_opt)
dot.ref <- setup.ref(
  ref_data = matrix_sc,
  ref_annotations = annotation_list,
  ref_subcluster_size = ref_subcluster_size,
  max_genes = max_genes,
  verbose = verbose
)

# Create DOT object
dot <- create.DOT(dot.srt, dot.ref)

# Run DOT
dot <- run.DOT.highresolution(dot, ratios_weight = ratios_weight)

# Write DOT results (cell type weights per Xenium cell) to file
write.csv(dot@weights, file=DOT_weights_path, row.names = TRUE)

# Write DOT results for plotting purposes
plot_data <- as.data.frame(spatial)
plot_data$celltype <- colnames(dot@weights)[apply(dot@weights, 1, which.max)]
write.csv(plot_data, file = DOT_plot_path, row.names = TRUE)
