load_atop_validation <- function() {
  validation_environment <- new.env(parent = globalenv())
  validation_environment$`%||%` <- function(lhs, rhs) {
    if (is.null(lhs) || length(lhs) == 0L) rhs else lhs
  }
  sys.source(
    testthat::test_path("..", "..", "R", "atop", "validation.R"),
    envir = validation_environment
  )
  validation_environment
}

test_that("replicate validation preserves the exact user order", {
  validation <- load_atop_validation()
  expect_equal(
    validation$atop_replicate_columns(c("Treatment_3", "Treatment_1", "Treatment_2"), "Condition", c("Treatment_1", "Treatment_2", "Treatment_3")),
    c("Treatment_3", "Treatment_1", "Treatment_2")
  )
})

test_that("replicate validation rejects duplicates and cross-group reuse", {
  validation <- load_atop_validation()
  expect_error(
    validation$atop_replicate_columns(c("Treatment_1", "Treatment_1"), "Condition", "Treatment_1"),
    "repeats columns"
  )

  input <- data.frame(sgRNA = "g1", Gene = "A", T1 = 1, C1 = 2)
  expect_error(
    validation$atop_validate_request(
      input, "sgRNA", "Gene", NULL,
      "T1", "T1", 100L, 1, 1L, "cpp", FALSE,
      NULL, NULL, 0
    ),
    "cannot reuse"
  )
})
