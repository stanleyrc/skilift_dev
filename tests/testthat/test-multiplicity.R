suppressWarnings(devtools::load_all())
library(testthat)
library(data.table)

# test <- function() { testthat::test_file("tests/testthat/test-multiplicity.R") }

# Create mock data in setup
setup({
})

teardown({
})

test_that("", {
})

## integration tests (only works on NYU hpc)
will_test_integration = FALSE
if (will_test_integration) {
test_that("lift_multiplicity works on real cohort", {
    clinical_pairs_path = "~/projects/Clinical_NYU/db/pairs.rds"
    clinical_pairs = readRDS(clinical_pairs_path)
    vip_sample = clinical_pairs[patient_id == "397089", ]
    cohort = Cohort$new(vip_sample)

    lift_multiplicity(cohort, "~/public_html/case-reports-data/")
})

}


## debug
DEBUG = FALSE
if (DEBUG) {

}
