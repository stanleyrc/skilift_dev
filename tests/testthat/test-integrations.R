suppressWarnings(devtools::load_all())

library(testthat)

# test <- function() { testthat::test_file("tests/testthat/test-integrations.R") }

setup({
  library(dplyr)
  all_expected_outputs <- c(
    "coverage.arrow",
    "coverage_cn_boxplot.json",
    "germline_altered_hist.json",
    "purple_sunrise_beta_gamma.png",
    "mutations.json",
    "sage.qc.json",
    "hetsnps.arrow",
    "purple_sunrise.json",
    "somatic_altered_hist.json",
    "coverage_cn_boxplot_original.png",
    "coverage_cn_boxplot_denoised.png",
    "allelic.json",
    "sbs_decomposed_prob.json",
    "hetsnps_major_hist.json",
    "hetsnps_minor_hist.json",
    "ppfit.json",
    "purple_sunrise_beta_gamma.png",
    "purple_sunrise_pp.png",
    "complex.json",
    "hetsnps_major_minor.png",
    "id_decomposed_prob.json",
    "germline_total_hist.json",
    "somatic_total_hist.json",
    "filtered.events.json",
    "metadata.json",
    "mutation_catalog.json",
    "id_mutation_catalog.json"
  )
  # Helper function to check if output files exist
  check_output_files <<- function(
    cohort,
    temp_dir,
    tumor_only = FALSE,
    output_names = all_expected_outputs
  ) {
    for (pair in cohort$inputs$pair) {
      for (output_name in output_names) {
        if (tumor_only) {
          if (grepl("germline", output_name)) {
            next
          }

          if (grepl("mutations.json", output_name)) {
            next
          }
        }
        expect_true(
          file.exists(file.path(temp_dir, pair, output_name)),
          info = paste("Output file", output_name, "for pair", pair, "not found.")
        )
      }
    }

    message("All expected output files are present.")
  }

  check_empty_cols <<- function(dt) {
    dt %>%
      select(where(~ all(is.na(.)))) %>%
      names()
  }
})

test_that("can lift_all heme cohort from gosh outputs csv without errors", {
  # Load real clinical pairs
  output_csv_path = system.file(
    "extdata",
    "test_data",
    "heme_outputs.csv",
    package = "Skilift"
  )

  csv = data.table(read.csv(output_csv_path, stringsAsFactors = FALSE))
  csv = csv[patient_id == csv$patient_id[1], ]
  
  cohort <- suppressWarnings(Cohort$new(
    csv,
    cohort_type = "paired"
  ))

  # Create temp directory for output
  temp_dir <- tempdir()

  suppressWarnings(lift_all(cohort, output_data_dir = temp_dir, cores = 2))
  
  expect_message(
    check_output_files(cohort, temp_dir, tumor_only = TRUE),
    "All expected output files are present."
  )
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

test_that("can lift_all ffpe cohort from gosh outputs csv without errors", {
  # Load real clinical pairs
  output_csv_path = system.file(
    "extdata",
    "test_data",
    "ffpe_outputs.csv",
    package = "Skilift"
  )

  csv = data.table(read.csv(output_csv_path, stringsAsFactors = FALSE))
  csv = csv[patient_id == csv$patient_id[11], ]
  
  cohort <- suppressWarnings(Cohort$new(
    csv,
    cohort_type = "paired"
  ))

  # Create temp directory for output
  temp_dir <- tempdir()

  suppressWarnings(lift_all(cohort, output_data_dir = temp_dir, cores = 2))
  
  expect_message(
    check_output_files(cohort, temp_dir),
    "All expected output files are present."
  )
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})


test_that("can lift_all ffpe cohort from pairs table without errors", {
  # Load real clinical pairs
  clinical_pairs_path = "~/projects/Clinical_NYU/db/pairs.rds"
  clinical_pairs = readRDS(clinical_pairs_path)
  pairs_with_all_outputs = c(
    "397089"
  )
 
  cohort <- suppressWarnings(Cohort$new(
    clinical_pairs[clinical_pairs$patient_id %in% pairs_with_all_outputs],
    col_mapping = list(pair = "patient_id", oncotable = "oncotable__wgs"),
    cohort_type = "paired"
  ))

  # Create temp directory for output
  temp_dir <- tempdir()

  suppressWarnings(lift_all(cohort, output_data_dir = temp_dir, cores = 2))
  
  expect_message(
    check_output_files(cohort, temp_dir),
    "All expected output files are present."
  )
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

test_that("can lift_all ALL cohort from gosh outputs csv without errors", {
  # Load real clinical pairs
  output_csv_path = system.file(
    "extdata",
    "test_data",
    "caroll_outputs.csv",
    package = "Skilift"
  )

  csv = data.table(read.csv(output_csv_path, stringsAsFactors = FALSE))
  csv = csv[patient_id == csv$patient_id[1], ]
  
  
  cohort <- suppressWarnings(Cohort$new(
    csv,
    cohort_type = "paired"
  ))

  # Create temp directory for output
  temp_dir <- tempdir()

  lift_all(cohort, output_data_dir = temp_dir, cores = 2)
  lift_hetsnps(cohort, output_data_dir = temp_dir, cores = 2)
  
  expect_message(
    check_output_files(cohort, temp_dir, tumor_only = TRUE),
    "All expected output files are present."
  )
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

# DEBUG
#########################################

# cohort$inputs
# names(cohort$inputs)
# cohort$inputs %>% select(where(~ all(is.na(.)))) %>% names

