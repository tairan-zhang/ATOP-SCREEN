load_visualization_parameter_resolver <- function() {
  worker_environment <- new.env(parent = globalenv())
  worker_environment$`%||%` <- function(lhs, rhs) {
    if (is.null(lhs) || length(lhs) == 0L) rhs else lhs
  }
  sys.source(
    testthat::test_path("..", "..", "inst", "worker", "modules", "visualization_commands.R"),
    envir = worker_environment
  )
  worker_environment$resolve_visualization_parameters
}

test_that("visualization colors and selected labels survive parameter resolution", {
  resolve_visualization_parameters <- load_visualization_parameter_resolver()
  params <- resolve_visualization_parameters(list(
    pathwayIds = c("PATHWAY_A", "PATHWAY_B"),
    pathwayColors = c("#112233", "#445566"),
    labelGenes = c("TP53", "MYC"),
    knownGenes = "TP53",
    highlightedGenes = c("MYC", "KRAS"),
    colUp = "#AA0011",
    gradientHigh = "#BB2233",
    positiveColor = "#CC4455",
    condition1Color = "#101010"
  ))

  expect_equal(params$pathwayIds, c("PATHWAY_A", "PATHWAY_B"))
  expect_equal(params$pathwayColors, c("#112233", "#445566"))
  expect_equal(params$labelGenes, c("TP53", "MYC"))
  expect_equal(params$knownGenes, "TP53")
  expect_equal(params$highlightedGenes, c("MYC", "KRAS"))
  expect_equal(params$colUp, "#AA0011")
  expect_equal(params$gradientHigh, "#BB2233")
  expect_equal(params$positiveColor, "#CC4455")
  expect_equal(params$condition1Color, "#101010")
})

test_that("invalid visualization colors fall back to stable defaults", {
  resolve_visualization_parameters <- load_visualization_parameter_resolver()
  params <- resolve_visualization_parameters(list(
    colUp = "not-a-colour",
    pathwayColors = c("#112233", "bad-pathway-colour")
  ))

  expect_equal(params$colUp, "#C02635")
  expect_equal(params$colDown, "#00609D")
  expect_equal(params$pathwayColors, c("#112233", "#C02635"))
  expect_equal(params$fontSizePt, 8)
  expect_equal(params$pngDpi, 300L)
})

test_that("each plot type receives its requested default dimensions", {
  resolve_visualization_parameters <- load_visualization_parameter_resolver()
  expect_equal(resolve_visualization_parameters(list(plotType = "gene_volcano"))[c("widthCm", "heightCm")], list(widthCm = 6.5, heightCm = 6.5))
  expect_equal(resolve_visualization_parameters(list(plotType = "gene_ranking"))[c("widthCm", "heightCm")], list(widthCm = 11, heightCm = 8.5))
  expect_equal(resolve_visualization_parameters(list(plotType = "sgrna_paired"))[c("widthCm", "heightCm")], list(widthCm = 6.5, heightCm = 6.5))
  expect_equal(resolve_visualization_parameters(list(plotType = "gsea_lollipop"))[c("widthCm", "heightCm")], list(widthCm = 13, heightCm = 11.5))
})

test_that("GSEA running-plot controls are normalized for the ATOP renderer", {
  resolve_visualization_parameters <- load_visualization_parameter_resolver()
  params <- resolve_visualization_parameters(list(
    plotType = "gsea_running",
    gseaTermWidth = 4,
    legendX = 2,
    legendY = -1,
    subplotType = 7,
    showNes = FALSE,
    showPvalue = TRUE,
    pvalueType = "pvalue",
    showRankExplanation = FALSE,
    rankPositiveLabel = "Positive end",
    rankNegativeLabel = "Negative end"
  ))

  expect_equal(params$gseaTermWidth, 8L)
  expect_equal(params$legendX, 1)
  expect_equal(params$legendY, 0)
  expect_equal(params$subplotType, 3L)
  expect_false(params$showNes)
  expect_true(params$showPvalue)
  expect_equal(params$pvalueType, "pvalue")
  expect_false(params$showRankExplanation)
  expect_equal(params$rankPositiveLabel, "Positive end")
  expect_equal(params$rankNegativeLabel, "Negative end")
})

test_that("GSEA lollipop titles are centered over the complete plot", {
  test_environment <- new.env(parent = globalenv())
  test_environment$`%||%` <- function(lhs, rhs) {
    if (is.null(lhs) || length(lhs) == 0L) rhs else lhs
  }
  sys.source(
    testthat::test_path("..", "..", "R", "plotting", "common", "typography.R"),
    envir = test_environment
  )
  sys.source(
    testthat::test_path("..", "..", "R", "plotting", "gsea", "lollipop.R"),
    envir = test_environment
  )
  plot <- test_environment$generate_gsea_lollipop_plot(
    data.frame(
      ID = c("PATHWAY_A", "PATHWAY_B"),
      NES = c(2, 1),
      pvalue = c(0.01, 0.05),
      p.adjust = c(0.02, 0.08),
      setSize = c(20, 10)
    ),
    "Positive", 2, "Positive pathway enrichment"
  )

  expect_equal(plot$theme$plot.title.position, "plot")
  expect_equal(plot$theme$plot.title$hjust, 0.5)
})
