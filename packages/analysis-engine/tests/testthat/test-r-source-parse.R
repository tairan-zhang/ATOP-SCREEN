test_that("all analysis R modules parse", {
  files <- list.files(testthat::test_path("..", "..", "R"), pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
  expect_gt(length(files), 0)
  for (file in files) expect_no_error(parse(file = file))
})
