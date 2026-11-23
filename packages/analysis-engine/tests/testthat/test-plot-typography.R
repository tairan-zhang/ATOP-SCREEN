load_plot_typography <- function() {
  plot_environment <- new.env(parent = globalenv())
  plot_environment$`%||%` <- function(lhs, rhs) {
    if (is.null(lhs) || length(lhs) == 0L) rhs else lhs
  }
  files <- c(
    file.path("common", "typography.R"),
    file.path("gene", "common.R"),
    file.path("gene", "volcano.R"),
    file.path("gene", "ranking.R"),
    file.path("gsea", "lollipop.R"),
    file.path("sgrna", "paired.R")
  )
  for (file in files) {
    sys.source(
      testthat::test_path("..", "..", "R", "plotting", file),
      envir = plot_environment
    )
  }
  plot_environment
}

expect_uniform_plot_theme <- function(plot, font_size_pt) {
  expect_equal(plot$theme$plot.title$size, font_size_pt)
  expect_equal(plot$theme$axis.text$size, font_size_pt)
  expect_equal(plot$theme$axis.title$size, font_size_pt)
  expect_equal(plot$theme$legend.title$size, font_size_pt)
  expect_equal(plot$theme$legend.text$size, font_size_pt)
}

test_that("the shared typography contract uses one point size everywhere", {
  plotter <- load_plot_typography()
  theme <- plotter$atop_plot_text_theme(8, "Helvetica")

  expect_equal(theme$text$size, 8)
  expect_equal(theme$plot.title$size, 8)
  expect_equal(theme$axis.text$size, 8)
  expect_equal(theme$axis.title$size, 8)
  expect_equal(theme$legend.title$size, 8)
  expect_equal(theme$legend.text$size, 8)
  expect_equal(theme$strip.text$size, 8)
  expect_equal(plotter$atop_font_size_mm(8), 8 / (72.27 / 25.4))
})

test_that("gene, pathway, and sgRNA builders honour the same font parameter", {
  plotter <- load_plot_typography()
  font_size <- 9
  gene_data <- data.frame(
    Gene = c("A", "B", "C"),
    GeneNegativeScore = c(-0.6, -0.3, -0.1),
    P_negative = c(0.001, 0.02, 0.5),
    GenePositiveScore = c(0.8, 0.4, 0.1),
    P_positive = c(0.001, 0.03, 0.6),
    GenePositiveRank = 1:3,
    GeneNegativeRank = 1:3
  )
  gsea_data <- data.frame(
    ID = c("PATHWAY_A", "PATHWAY_B"),
    Description = c("PATHWAY_A", "PATHWAY_B"),
    NES = c(1.8, 1.2),
    pvalue = c(0.001, 0.02),
    p.adjust = c(0.01, 0.04),
    setSize = c(20, 30)
  )
  counts <- data.frame(
    guide = c("g1", "g2"),
    gene = c("A", "A"),
    treatment = c(10, 20),
    control = c(8, 12)
  )

  plots <- list(
    plotter$generate_volcano_plot(gene_data, label_genes = "A", base_font_size = font_size),
    plotter$generate_ranking_plot(gene_data, known_genes = "A", base_font_size = font_size),
    plotter$generate_gsea_lollipop_plot(gsea_data, "Positive", 2, "Pathways", base_font_size = font_size),
    plotter$generate_sgrna_paired_plot(counts, "A", "guide", "gene", "treatment", "control", base_font_size = font_size)
  )

  lapply(plots, expect_uniform_plot_theme, font_size_pt = font_size)
  expect_equal(plots[[1]]$layers[[length(plots[[1]]$layers)]]$aes_params$size, plotter$atop_font_size_mm(font_size))
  expect_equal(plots[[2]]$layers[[2]]$aes_params$size, plotter$atop_font_size_mm(font_size))
})

test_that("volcano legend follows the ATOP-SCREEN bottom layout", {
  plotter <- load_plot_typography()
  gene_data <- data.frame(
    Gene = c("A", "B"),
    GeneNegativeScore = c(-0.8, -0.1),
    P_negative = c(0.001, 0.5),
    GenePositiveScore = c(0.8, 0.1),
    P_positive = c(0.001, 0.5),
    GenePositiveRank = 1:2,
    GeneNegativeRank = 1:2
  )

  plot <- plotter$generate_volcano_plot(gene_data, label_genes = "A", base_font_size = 8)
  expect_equal(plot$theme$legend.position, "bottom")
  expect_equal(plot$theme$legend.direction, "horizontal")
  expect_s3_class(plot$layers[[length(plot$layers)]]$geom, "GeomTextRepel")
})
