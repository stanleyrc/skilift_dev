suppressWarnings(devtools::load_all())

setwd("~/Projects/skilift/")

library(testthat)

# test_all <- function() { testthat::test_file("tests/run-all-tests.R") }

# testthat::test_dir("tests/testthat", reporter = testthat::ListReporter)

# testthat::test_file("tests/testthat/test-CaseReportsData.R")
testthat::test_file("tests/testthat/test-Cohort.R")
testthat::test_file("tests/testthat/test-copy-number-graph.R")
testthat::test_file("tests/testthat/test-filtered-events.R")
testthat::test_file("tests/testthat/test-metadata.R")
testthat::test_file("tests/testthat/test-scatterplot.R") #2
testthat::test_file("tests/testthat/test-segment-width-distribution.R")
testthat::test_file("tests/testthat/test-variant-qc.R")
testthat::test_file("tests/testthat/test-multiplicity.R") #9
testthat::test_file("tests/testthat/test-signatures.R")
testthat::test_file("tests/testthat/test-lift-wrappers.R") #4


