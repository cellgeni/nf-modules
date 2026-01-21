#!/usr/bin/env Rscript

library(DOTr)
library(ggplot2)
library(anndata)
library(reticulate)
library(Matrix)
library(tools)


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_dot.R <adata_sc_path> <adata_sp_path> <output_path>")
}

adata_sc_path <- args[1]
adata_sp_path <- args[2]
output_path <- args[3]

# Extract base name
sp_base <- file_path_sans_ext(basename(adata_sp_path))

# Define subfolders
weights_dir <- file.path(output_path, "DOT_weights")
plot_dir <- file.path(output_path, "DOT_plotdata")

# Create directories if they don't exist
dir.create(weights_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# Build full output file paths
DOT_weights_path <- file.path(weights_dir, paste0("DOT_weights_", sp_base, ".csv"))
DOT_plot_path <- file.path(plot_dir, paste0("DOT_plotting_", sp_base, ".csv"))

# Load expression data
adata_sc <- read_h5ad(adata_sc_path)
adata_sp <- read_h5ad(adata_sp_path)

# Setup DOT with correct expression data and spatial coordinates
expression <- adata_sp$X
spatial <- adata_sp$obs[,c('x_centroid', 'y_centroid')]
dot.srt <- setup.srt(t(expression), spatial, th.gene.low=0, th.gene.high=1)

# Select reference cell types from adata_sc
annotation_list <- as.character(adata_sc$obs$identities) ## Adjust this if celltype annotations are stored differently
matrix_sc <- as(t(adata_sc$X), "CsparseMatrix")
dot.ref <- setup.ref(ref_data = matrix_sc, ref_annotations = annotation_list, ref_subcluster_size = 3, max_genes=as.integer(adata_sc$n_vars), verbose=TRUE)

# Create DOT object
dot <- create.DOT(dot.srt, dot.ref)

# Run DOT
dot <- run.DOT.highresolution(dot, ratios_weight = 0)

# Write DOT results (cell type weights per Xenium cell) to file
write.csv(dot@weights, file=DOT_weights_path, row.names = TRUE)

# Write DOT results for plotting purposes
plot_data <- as.data.frame(spatial)
plot_data$celltype <- colnames(dot@weights)[apply(dot@weights, 1, which.max)]
write.csv(plot_data, file = DOT_plot_path, row.names = TRUE)