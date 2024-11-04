suppressWarnings(devtools::load_all())


library(testthat)

test_that("process_gencode handles NULL input", {
  result <- process_gencode(NULL)
  expect_true(is(result, "GRanges"))
})

test_that("process_gencode handles .rds input", {
  test_rds_path <- system.file("extdata/test_data/test_gencode_v29lift37.rds", package = "Skilift")
  result <- process_gencode(test_rds_path)
  expect_true(is(result, "GRanges"))
})

test_that("process_gencode handles GTF input", {
  # Assuming 'test_gencode.gtf' is a valid GTF file path for testing
  test_gtf_path <- "path/to/test_gencode.gtf"
  result <- process_gencode(test_gtf_path)
  expect_true(is(result, "GRanges"))
})
