#!/usr/bin/env Rscript

options(bitmapType='cairo')

# Load packages
suppressPackageStartupMessages({
	library(Seurat)
	library(ggplot2)
	library(argparse)
})


# ---- Argument Parsing ----
parser <- ArgumentParser(description = "Generate violin plots from Seurat object")

parser$add_argument("--genes", type = "character", required = TRUE,
		                        help = "Comma-separated gene names (e.g., rab-3,twk-18)")
parser$add_argument("--neurons", type = "character", required = TRUE,
		                        help = "Comma-separated neuron names (e.g., AIY,RIG,ASH)")
parser$add_argument("--output", type = "character", required = TRUE,
		                        help = "Path to output PNG file")
parser$add_argument("--differential", action = "store_true", default = FALSE,
		                        help = "Include to split violin plots by genotype")
parser$add_argument("--dataset", action = "store_true", default = FALSE, help = "Choose between the N2/daf-2 dataset or the Male dataset")
args <- parser$parse_args()

# ---- Parse Arguments ----
genes <- strsplit(args$genes, ",")[[1]]
neurons <- strsplit(args$neurons, ",")[[1]]
output_file <- args$output
split_by <- if (args$differential) "genotype" else NULL
dataset <- args$dataset
# ---- Load Seurat Object ----
if (args$dataset) {
	seurat_obj <- readRDS("/var/www/ctvm1/murphy-lab-project/source_data/Male_herm.rds")
} else {
	seurat_obj <- readRDS("/var/www/ctvm1/murphy-lab-project/source_data/WT_daf2.rds")
}

#seurat_obj <- readRDS("/var/www/ctvm1/murphy-lab-project/source_data/WT_daf2.rds")

plot <- VlnPlot(
		seurat_obj,
		features = genes,
		assay = "SCT",
		idents = neurons,
		split.by = split_by,
		cols = c("red","grey")
	)


# ---- Save to PNG ----
ggsave(filename = output_file, plot = plot, width = 10, height = 6, dpi = 300)












