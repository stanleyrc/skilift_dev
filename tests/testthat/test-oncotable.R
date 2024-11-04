suppressWarnings(devtools::load_all())


library(testthat)

test_that("process_gencode handles NULL input", {
  result <- process_gencode(NULL)
  expect_error(result, "gencode file must be provided")
})

test_that("process_gencode handles .rds input", {
  test_rds_path <- system.file("extdata/test_data/test_gencode_v29lift37.rds", package = "Skilift")
  result <- process_gencode(test_rds_path)
  expect_true(is(result, "GRanges"))
})
