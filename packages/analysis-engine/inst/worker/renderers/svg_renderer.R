#!/usr/bin/env Rscript

# SVG rendering runs in a clean R process. Bioconductor's GO.db namespace,
# loaded during GSEA, conflicts with svglite in the bundled runtime and can
# otherwise make even small vector plots exhaust the macOS R memory limit.

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 4L) {
    stop("Usage: svg_renderer.R <plot.rds> <output.svg> <width-in> <height-in>", call. = FALSE)
}

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggrepel)
    library(ggpp)
    library(patchwork)
    library(svglite)
})

plot_path <- normalizePath(enc2utf8(arguments[[1]]), mustWork = TRUE)
output_path <- enc2utf8(arguments[[2]])
width <- as.numeric(arguments[[3]])
height <- as.numeric(arguments[[4]])
if (!is.finite(width) || !is.finite(height) || width <= 0 || height <= 0) {
    stop("SVG dimensions must be positive finite numbers.", call. = FALSE)
}

plot_object <- readRDS(plot_path)
svglite::svglite(output_path, width = width, height = height)
tryCatch(
    print(plot_object),
    finally = grDevices::dev.off()
)
if (!file.exists(output_path) || file.info(output_path)$size <= 0L) {
    stop("SVG renderer did not create a valid output file.", call. = FALSE)
}
