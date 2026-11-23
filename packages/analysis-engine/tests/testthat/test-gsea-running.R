load_gsea_running_plotter <- function() {
  plot_environment <- new.env(parent = globalenv())
  plot_environment$`%||%` <- function(lhs, rhs) {
    if (is.null(lhs) || length(lhs) == 0L) rhs else lhs
  }
  sys.source(
    testthat::test_path("..", "..", "R", "plotting", "common", "typography.R"),
    envir = plot_environment
  )
  for (file in c("running_helpers.R", "running_gseavis.R", "running_style.R", "running.R")) {
    sys.source(
      testthat::test_path("..", "..", "R", "plotting", "gsea", file),
      envir = plot_environment
    )
  }
  plot_environment
}

test_that("GseaVis running plot accepts an empty highlighted-gene selection", {
  skip_if_not_installed("GseaVis")
  plotter <- load_gsea_running_plotter()
  example_path <- system.file("extdata", "gseaRes.RDS", package = "GseaVis")
  skip_if(!nzchar(example_path))
  gsea_result <- readRDS(example_path)
  pathway_id <- as.character(as.data.frame(gsea_result)$ID[[1]])

  plot <- expect_no_error(plotter$generate_gsea_enrichment_plot_atop(
    gsea_object = gsea_result,
    pathway_ids = pathway_id,
    highlighted_genes = character()
  ))
  expect_false(is.null(plot))
  expect_s3_class(plot, "patchwork")
})

test_that("highlight genes are validated against core enrichment genes", {
  skip_if_not_installed("GseaVis")
  plotter <- load_gsea_running_plotter()
  example_path <- system.file("extdata", "gseaRes.RDS", package = "GseaVis")
  skip_if(!nzchar(example_path))
  gsea_result <- readRDS(example_path)
  pathway_id <- as.character(as.data.frame(gsea_result)$ID[[1]])

  expect_error(
    plotter$resolve_gsea_highlighted_genes_atop(gsea_result, pathway_id, "NOT_IN_THIS_PATHWAY"),
    "None of the specified highlight genes"
  )
})

test_that("pathway labels are normalized without changing the source object", {
  skip_if_not_installed("GseaVis")
  plotter <- load_gsea_running_plotter()
  example_path <- system.file("extdata", "gseaRes.RDS", package = "GseaVis")
  skip_if(!nzchar(example_path))
  gsea_result <- readRDS(example_path)
  original_id <- as.character(gsea_result@result$ID[[1]])

  normalized <- plotter$normalize_gsea_pathway_labels_atop(gsea_result)
  expect_equal(as.character(gsea_result@result$ID[[1]]), original_id)
  expect_equal(as.character(normalized$object@result$ID[[1]]), toupper(gsub("_", " ", original_id, fixed = TRUE)))
})

test_that("running curves retain hits while removing redundant collinear ranks", {
  plotter <- load_gsea_running_plotter()
  curve <- data.frame(
    x = 1:8,
    runningScore = c(-1, -2, 3, 2, 1, 0, 4, 3),
    position = c(0, 0, 1, 0, 0, 0, 1, 0),
    id = "PATHWAY"
  )

  simplified <- plotter$simplify_gsea_running_curve_data_atop(curve)
  expect_true(all(c(2, 3, 4, 6, 7, 8) %in% simplified$x))
  expect_lt(nrow(simplified), nrow(curve))
})
