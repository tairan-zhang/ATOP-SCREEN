ncs_env <- new.env(parent = baseenv())
ncs_env$median <- stats::median
sys.source(testthat::test_path("..", "..", "R", "atop", "preprocessing.R"), ncs_env)

test_that("NCS respects the neutral boundary, ignores nonfinite LFCs, and uses raw MAD", {
    settings <- list(ncs_delta0 = 0.1, ncs_delta1 = 0.3, ncs_delta2 = 0.2)
    lfc <- rbind(c(-0.1, 0, 0.1), c(1, 2, 3), c(0.15, NA, Inf), c(NaN, NA, Inf))
    actual <- ncs_env$atop_ncs(lfc, settings)
    expect_equal(unname(actual[, 1]), c(1, 1/3, 2/3, 2/9))
    expect_equal(unname(actual[, 2]), c(1, 1, 1/6, 1/6))
    expect_equal(unname(actual[, 3]), c(1, 0.5, 1, 0.5))
    expect_true(all(is.na(actual[, 4])))
})

test_that("NCS cutoff retains equality and reports completely removed genes", {
    guides <- data.frame(Gene = c("A", "A", "B"), NCS = c(0.1, 0.4, 0))
    filter <- function(cutoff) ncs_env$atop_filter_ncs(guides, list(ncs_threshold = cutoff, gene_col = "Gene"))
    expect_identical(filter(0)$ranked, guides)
    expect_equal(filter(0.4)$ranked$NCS, 0.4)
    expect_equal(filter(0.4)$excluded$NCS, c(0.1, 0))
    expect_equal(filter(0.4)$qc$retained_sgrnas[filter(0.4)$qc$Gene == "B"], 0)
})
