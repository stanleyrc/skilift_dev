suppressWarnings(devtools::load_all())

library(testthat)

# test <- function() { testthat::test_file("tests/testthat/test-plots.R") }

test_that("CaseReportsData constructor handles various directory scenarios", {
  # Setup test directories and files
  base_dir <- tempfile("test_cases")
  dir.create(base_dir)
  on.exit(unlink(base_dir, recursive = TRUE))
  
  # Test 1: Empty directory
  expect_error(
    CaseReportsData$new(base_dir),
    "Case reports directory does not exist:"
  )
  
  # Create test case directories
  case_complete <- file.path(base_dir, "case1")
  case_missing <- file.path(base_dir, "case2")
  case_extra <- file.path(base_dir, "case3")
  case_mixed <- file.path(base_dir, "case4")
  
  dir.create(case_complete)
  dir.create(case_missing)
  dir.create(case_extra)
  dir.create(case_mixed)
  
  # Expected files
  expected_files <- c(
    "metadata.json",
    "coverage.arrow",
    "sage_qc.json",
    "filtered.events.json",
    "hetsnps.arrow",
    "allelic.json",
    "complex.json",
    "ppfit.json",
    "mutation_catalog.json",
    "id_mutation_catalog.json",
    "sbs_decomposed_prob.json",
    "id_decomposed_prob.json"
  )
  
  # Create complete case
  for(file in expected_files) {
    writeLines("", file.path(case_complete, file))
  }
  
  # Create case with missing files
  present_files <- expected_files[1:5]  # Only create first 5 files
  for(file in present_files) {
    writeLines("", file.path(case_missing, file))
  }
  
  # Create case with extra files
  for(file in expected_files) {
    writeLines("", file.path(case_extra, file))
  }
  writeLines("", file.path(case_extra, "extra1.json"))
  writeLines("", file.path(case_extra, "extra2.json"))
  
  # Create case with both missing and extra files
  mixed_files <- c(expected_files[1:8], "extra1.json", "extra2.json")
  for(file in mixed_files) {
    writeLines("", file.path(case_mixed, file))
  }
  
  # Test 2: Directory with all files
  complete_data <- suppressWarnings(CaseReportsData$new(base_dir))
  expect_true(all(!is.na(complete_data$case_reports[pair == "case1"])))
  expect_equal(nrow(complete_data$case_reports), 4)
  
  # Test 3: Check warnings for missing files
  expect_warning(
    CaseReportsData$new(base_dir),
    "Missing ppfit file in cases: case2, case4"
  )
  expect_warning(
    CaseReportsData$new(base_dir),
    "Missing mutation_catalog file in cases: case2, case4"
  )
  
  # Test 4: Check warnings for extra files
  expect_warning(
    CaseReportsData$new(base_dir),
    "Unexpected file found in case case3: extra1"
  )
  expect_warning(
    CaseReportsData$new(base_dir),
    "Unexpected file found in case case4: extra1"
  )
  
  # Test 5: Verify data structure
  data <- suppressWarnings(CaseReportsData$new(base_dir))
  expected_cols <- c(tools::file_path_sans_ext(expected_files), "pair")
  expect_equal(sort(names(data$case_reports)), sort(expected_cols))
  
  # Test 6: Verify missing values are NA
  missing_case <- data$case_reports[pair == "case2"]
  expect_true(all(is.na(missing_case[, 6:12, with = FALSE])))
  
  # Test 7: Verify complete values are not NA
  complete_case <- data$case_reports[pair == "case1"]
  expect_true(all(!is.na(complete_case)))
})
