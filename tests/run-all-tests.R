suppressWarnings(devtools::load_all())

setwd("~/git/skilift/")

library(testthat)

# test_all <- function() { testthat::test_file("tests/run-all-tests.R") }

testthat::test_dir("tests/testthat", reporter = testthat::ListReporter)

# testthat::test_file("tests/testthat/test-CaseReportsData.R")
testthat::test_file("tests/testthat/test-Cohort.R")
testthat::test_file("tests/testthat/test-copy-number-graph.R")
testthat::test_file("tests/testthat/test-filtered-events.R")
testthat::test_file("tests/testthat/test-metadata.R")
testthat::test_file("tests/testthat/test-scatterplot.R")
testthat::test_file("tests/testthat/test-segment-width-distribution.R")
testthat::test_file("tests/testthat/test-variant-qc.R")
testthat::test_file("tests/testthat/test-multiplicity.R")
testthat::test_file("tests/testthat/test-signatures.R")
testthat::test_file("tests/testthat/test-lift-wrappers.R")


