#!/usr/local/R/4.5.0/bin/Rscript

#/usr/bin/env Rscript

options(bitmapType='cairo')

CESTAAN_ROOT <- Sys.getenv("PWD")

# Load necessary packages
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(argparse)
})

# Setup command line argument parser
parser <- ArgumentParser(description = "Generate UMAP with Seurat FeaturePlot")
parser$add_argument("gene1", help = "First gene name (or only gene if plotting one)")
parser$add_argument("output_file", help = "Path to save output image (PNG)")
parser$add_argument("--gene2", default = NULL, help = "Optional second gene for blended FeaturePlot")
parser$add_argument("--dataset", action = "store_true", default = FALSE, help = "Choose the dataset you'd like to access")
args <- parser$parse_args()


gene1 <- args$gene1
gene2 <- args$gene2
output_file <- args$output_file
dataset <- args$dataset

#args <- commandArgs(trailingOnly = TRUE)
#gene1 <- args[1]
#output_file <- args[2]
#gene2 <- if (length(args) >= 3) args[3] else NULL

if (args$dataset) {
	seurat_obj <- readRDS(paste(CESTAAN_ROOT, "source_data/Male_herm.rds", sep='/'))
} else {
	seurat_obj <- readRDS(paste(CESTAAN_ROOT, "source_data/WT_daf2.rds", sep='/'))
}
print(class(seurat_obj))

print(Seurat::Reductions(seurat_obj))
DefaultAssay(seurat_obj) <- "SCT"


# Check if gene exists
if (!(gene1 %in% rownames(seurat_obj))) {
	stop(paste("Gene", gene1, "not found in dataset"))
}
if (!is.null(gene2) && !(gene2 %in% rownames(seurat_obj))) {
	stop(paste("Gene", gene2, "not found in dataset"))
}


plot <- NULL
# Generate plot
if (!is.null(args$gene2)) {
  plot <- FeaturePlot(
    seurat_obj,
    features = c(args$gene1, args$gene2),
    min.cutoff = "q9",
    blend = TRUE,
    order = TRUE,
    pt.size = 0.5,
    label= TRUE,
    repel = TRUE,
    raster = FALSE
  ) + ggtitle(paste(args$gene1, "+", args$gene2))
} else {
  plot <- FeaturePlot(
    seurat_obj,
    features = args$gene1,
    min.cutoff = "q9",
    pt.size = 0.5,
    label=TRUE,
    repel=TRUE,
    order = TRUE,
    raster = FALSE
  ) + ggtitle(args$gene1)
}

#print(typeof(plot))
#print("UMAP Plotted")
# Save
#ggsave(filename = output_file,plot=plot, width = 6, height = 5, dpi = 300)
if (!is.null(gene2)) {
	  # Rectangular image for blended plot
	  ggsave(filename = output_file, plot = plot, width = 18, height = 7, dpi = 300)
} else {
	  # Square image for single gene
	  ggsave(filename = output_file, plot = plot, width = 9, height = 7, dpi = 300)
}

#ggsave(filename = args$output_file, plot = plot, width = 7, height = 6, dpi = 300)
