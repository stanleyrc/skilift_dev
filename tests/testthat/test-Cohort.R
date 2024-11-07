library(testthat)
library(data.table)

test_that("Cohort constructor handles various data.table inputs correctly", {
  # Test empty data.table
  empty_dt <- data.table()
  expect_warning(
    cohort <- Cohort$new(empty_dt),
    "No data could be extracted from input data.table"
  )
  expect_equal(nrow(cohort$inputs), 0)
  
  # Test data.table with missing columns
  partial_dt <- data.table(
    pair = c("sample1", "sample2"),
    tumor = c("tumor1", "tumor2")
  )
  expect_warning(
    cohort <- Cohort$new(partial_dt),
    "No matching column found for"
  )
  expect_equal(ncol(cohort$inputs), 2)
  expect_true(all(c("pair", "tumor") %in% names(cohort$inputs)))
  
  # Test data.table with all columns
  complete_dt <- data.table(
    pair = c("sample1", "sample2"),
    tumor = c("tumor1", "tumor2"),
    normal = c("normal1", "normal2"),
    annotated_bcf = c("bcf1", "bcf2"),
    fusions = c("fusion1", "fusion2"),
    jabba_simple = c("jabba1", "jabba2"),
    karyograph = c("kg1", "kg2"),
    events = c("event1", "event2"),
    signature_counts = c("sig1", "sig2"),
    oncokb_maf = c("maf1", "maf2"),
    oncokb_cna = c("cna1", "cna2")
  )
  expect_silent(cohort <- Cohort$new(complete_dt))
  expect_equal(ncol(cohort$inputs), ncol(complete_dt))
  
  # Test data.table with extra columns
  extra_dt <- copy(complete_dt)
  extra_dt[, extra_col := c("extra1", "extra2")]
  expect_silent(cohort <- Cohort$new(extra_dt))
  expect_equal(ncol(cohort$inputs), ncol(complete_dt))
  
  # Test data.table with extra columns and missing columns
  mixed_dt <- data.table(
    pair = c("sample1", "sample2"),
    tumor = c("tumor1", "tumor2"),
    extra_col = c("extra1", "extra2")
  )
  expect_warning(
    cohort <- Cohort$new(mixed_dt),
    "No matching column found for"
  )
  expect_equal(ncol(cohort$inputs), 2)
  
  # Test data.table with two columns mapping to same input column
  duplicate_dt <- data.table(
    pair = c("sample1", "sample2"),
    sample_id = c("sample1", "sample2")  # This should map to 'pair' as well
  )
  expect_silent(cohort <- Cohort$new(duplicate_dt))
  expect_equal(ncol(cohort$inputs), 1)
  expect_true("pair" %in% names(cohort$inputs))
  
  # Test custom col_mapping for new column
  custom_dt <- data.table(
    custom_col = c("custom1", "custom2")
  )
  custom_mapping <- list(
    new_col = c("custom_col")
  )
  expect_silent(
    cohort <- Cohort$new(custom_dt, col_mapping = custom_mapping)
  )
  expect_true("new_col" %in% names(cohort$inputs))
  
  # Test custom col_mapping that overrides default mapping
  override_dt <- data.table(
    alt_pair = c("sample1", "sample2")
  )
  override_mapping <- list(
    pair = c("alt_pair", "pair", "sample")  # Changed order from default
  )
  expect_silent(
    cohort <- Cohort$new(override_dt, col_mapping = override_mapping)
  )
  expect_true("pair" %in% names(cohort$inputs))
  expect_equal(cohort$inputs$pair, override_dt$alt_pair)
})
