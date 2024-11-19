suppressWarnings(devtools::load_all())

library(testthat)

# test_all <- function() { testthat::test_file("tests/run-all-tests.R") }

testthat::test_file("tests/testthat/test-Cohort.R")
testthat::test_file("tests/testthat/test-filtered-events.R")
testthat::test_file("tests/testthat/test-metadata.R")
