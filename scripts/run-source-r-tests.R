#!/usr/bin/env Rscript

if (!requireNamespace("testthat", quietly = TRUE)) {
  stop(
    "The source-test runner requires the development-only testthat package. ",
    "It is intentionally not bundled in the application runtime."
  )
}

testthat::test_dir(
  "packages/analysis-engine/tests/testthat",
  reporter = "summary"
)
