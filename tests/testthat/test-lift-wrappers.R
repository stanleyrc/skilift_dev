suppressWarnings(devtools::load_all())

library(testthat)
library(data.table)

test_that("lift_all handles basic case correctly", {
  # Create temp directory structure
  base_dir <- tempdir()
  output_dir <- file.path(base_dir, "output")
  suppressWarnings(dir.create(output_dir, recursive = TRUE))
  on.exit(unlink(base_dir, recursive = TRUE))
  
  # Setup basic test data
  dt <- data.table(
    pair = c("sample1", "sample2"),
    tumor_type = c("BRCA", "LUAD"),
    disease = c("Breast", "Lung"),
    primary_site = c("Breast", "Lung"),
    inferred_sex = c("F", "M")
  )
  
  cohort <- suppressWarnings(Cohort$new(dt))
  
  # Test basic execution with defaults
  expect_message(
    suppressWarnings(lift_all(
      cohort = cohort,
      output_data_dir = output_dir
    )),
    "Uploading in paired mode"
  )
})

test_that("lift_all passes parameters correctly to individual lifters", {
  # Create temp directory structure
  base_dir <- tempdir()
  output_dir <- file.path(base_dir, "output")
  suppressWarnings(dir.create(output_dir, recursive = TRUE))
  on.exit(unlink(base_dir, recursive = TRUE))
  
  # Setup test data with all required fields
  dt <- data.table(
    pair = c("sample1", "sample2"),
    tumor_type = c("BRCA", "LUAD"),
    disease = c("Breast", "Lung"),
    primary_site = c("Breast", "Lung"),
    inferred_sex = c("F", "M"),
    structural_variants = c("sv1.vcf", "sv2.vcf"),
    tumor_coverage = c("cov1.bw", "cov2.bw"),
    somatic_snvs = c("snv1.vcf", "snv2.vcf"),
    germline_snvs = c("germ1.vcf", "germ2.vcf"),
    het_pileups = c("het1.txt", "het2.txt")
  )
  
  # Create mock files
  sapply(c(dt$structural_variants, dt$tumor_coverage, 
          dt$somatic_snvs, dt$germline_snvs, dt$het_pileups), 
         function(x) writeLines("", file.path(base_dir, x)))
  
  cohort <- suppressWarnings(Cohort$new(dt))
  
  # Test with custom parameters
  custom_annotations <- list(genes = c("BRCA1", "BRCA2"))
  
  expect_message(
    suppressWarnings(lift_all(
      cohort = cohort,
      output_data_dir = output_dir,
      cores = 2
    )),
    "Uploading in paired mode"
  )
})

test_that("lift_all handles missing optional parameters gracefully", {
  # Create temp directory structure
  base_dir <- tempdir()
  output_dir <- file.path(base_dir, "output")
  suppressWarnings(dir.create(output_dir, recursive = TRUE))
  on.exit(unlink(base_dir, recursive = TRUE))
  
  # Setup minimal test data
  dt <- data.table(
    pair = c("sample1"),
    tumor_type = c("BRCA")
  )
  
  cohort <- suppressWarnings(Cohort$new(dt))
  
  # Test with minimal parameters
  expect_message(
    suppressWarnings(lift_all(
      cohort = cohort,
      output_data_dir = output_dir
    )),
    "Uploading in paired mode"
  )
})

test_that("lift_all validates required parameters", {
  # Create temp directory
  base_dir <- tempdir()
  on.exit(unlink(base_dir, recursive = TRUE))
  
  # Setup basic test data
  dt <- data.table(
    pair = c("sample1"),
    tumor_type = c("BRCA")
  )
  
  cohort <- suppressWarnings(Cohort$new(dt))
  
  # Test missing required parameters
  expect_error(lift_all())  # Missing both required parameters
  expect_error(lift_all(cohort = cohort))  # Missing output_data_dir
  expect_error(lift_all(output_data_dir = base_dir))  # Missing cohort
})

test_that("lift_all handles different cohort_types correctly", {
  # Create temp directory structure
  base_dir <- tempdir()
  output_dir <- file.path(base_dir, "output")
  suppressWarnings(dir.create(output_dir, recursive = TRUE))
  on.exit(unlink(base_dir, recursive = TRUE))
  
  # Setup basic test data
  dt <- data.table(
    pair = c("sample1", "sample2"),
    tumor_type = c("BRCA", "LUAD"),
    disease = c("Breast", "Lung"),
    primary_site = c("Breast", "Lung"),
    inferred_sex = c("F", "M")
  )

  cohort <- suppressWarnings(Cohort$new(dt))
  tumor_only_cohort <- suppressWarnings(Cohort$new(dt, cohort_type = "tumor_only"))
  heme_cohort <- suppressWarnings(Cohort$new(dt, cohort_type = "heme"))
  
  # Test paired mode (default)
  expect_message(
    suppressWarnings(lift_all(
      cohort = cohort,
      output_data_dir = output_dir
    )),
    "Uploading in paired mode"
  )
  
  # Test tumor-only mode
  expect_message(
    suppressWarnings(lift_all(
      cohort = tumor_only_cohort,
      output_data_dir = output_dir,
    )),
    "Uploading in tumor_only mode"
  )
  
  # Test heme mode
  expect_message(
    suppressWarnings(lift_all(
      cohort = heme_cohort,
      output_data_dir = output_dir,
      cohort_type = "heme"
    )),
    "Uploading in heme mode"
  )
  
  # Test invalid cohort_type
  expect_error(
    invalid_type_cohort <- suppressWarnings(Cohort$new(dt, cohort_type = "invalid"))
  )
})
