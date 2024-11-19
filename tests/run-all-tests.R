suppressWarnings(devtools::load_all())

library(testthat)

test <- function() { testthat::test_file("tests/run-all-tests.R") }
