load_gene_ranking_plotter <- function() {
  plot_environment <- new.env(parent = globalenv())
  plot_environment$`%||%` <- function(lhs, rhs) {
    if (is.null(lhs) || length(lhs) == 0L) rhs else lhs
  }
  sys.source(
    testthat::test_path("..", "..", "R", "plotting", "common", "typography.R"),
    envir = plot_environment
  )
  sys.source(
    testthat::test_path("..", "..", "R", "plotting", "gene", "common.R"),
    envir = plot_environment
  )
  sys.source(
    testthat::test_path("..", "..", "R", "plotting", "gene", "ranking.R"),
    envir = plot_environment
  )
  plot_environment$generate_ranking_plot
}

test_that("ranking plot marks only the known subset without recycling coordinates", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("ggplot2")
  generate_ranking_plot <- load_gene_ranking_plotter()
  ranking_data <- data.frame(
    Gene = paste0("GENE", 1:5),
    GenePositiveRank = 1:5,
    GenePositiveScore = c(1.8, 1.5, 1.2, 0.9, 0.6),
    GeneNegativeRank = 1:5,
    GeneNegativeScore = -c(1.8, 1.5, 1.2, 0.9, 0.6)
  )

  plot <- generate_ranking_plot(
    ranking_data,
    type = "positive",
    top_n = 5,
    known_genes = c("GENE2", "GENE4")
  )
  built_plot <- expect_no_error(ggplot2::ggplot_build(plot))
  labels <- built_plot$data[[2]]$label
  expect_identical(labels[endsWith(labels, "*")], c("#2*", "#4*"))
})

test_that("negative ranking plot uses subset-aligned marker coordinates", {
  skip_if_not_installed("dplyr")
  skip_if_not_installed("ggplot2")
  generate_ranking_plot <- load_gene_ranking_plotter()
  ranking_data <- data.frame(
    Gene = paste0("GENE", 1:4),
    GenePositiveRank = 1:4,
    GenePositiveScore = c(1.8, 1.5, 1.2, 0.9),
    GeneNegativeRank = 1:4,
    GeneNegativeScore = -c(1.8, 1.5, 1.2, 0.9)
  )

  plot <- generate_ranking_plot(
    ranking_data,
    type = "negative",
    top_n = 4,
    known_genes = "GENE3"
  )
  built_plot <- expect_no_error(ggplot2::ggplot_build(plot))
  labels <- built_plot$data[[2]]$label
  expect_identical(labels[endsWith(labels, "*")], "#3*")
  expect_equal(as.numeric(built_plot$data[[3]]$label), ranking_data$GeneNegativeScore)
})
