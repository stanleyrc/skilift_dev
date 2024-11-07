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
    tumor_type = c("BRCA", "LUAD")
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
    disease = c("Breast", "Lung"),
    primary_site = c("Breast", "Lung"),
    inferred_sex = c("F", "M"),
    structural_variants = c("sv1", "sv2"),
    tumor_coverage = c("cov1", "cov2"),
    somatic_snvs = c("snv1", "snv2"),
    germline_snvs = c("germ1", "germ2"),
    het_pileups = c("het1", "het2"),
    somatic_snv_cn = c("scn1", "scn2"),
    germline_snv_cn = c("gcn1", "gcn2"),
    somatic_variant_annotations = c("sva1", "sva2"),
    germline_variant_annotations = c("gva1", "gva2"),
    oncokb_snv = c("osnv1", "osnv2"),
    oncokb_cna = c("ocna1", "ocna2"),
    jabba_gg = c("jgg1", "jgg2"),
    karyograph = c("kg1", "kg2"),
    balanced_jabba_gg = c("bjg1", "bjg2"),
    events = c("ev1", "ev2"),
    fusions = c("fus1", "fus2"),
    allelic_jabba_gg = c("ajg1", "ajg2"),
    activities_sbs_signatures = c("ass1", "ass2"),
    matrix_sbs_signatures = c("mss1", "mss2"),
    activities_indel_signatures = c("ais1", "ais2"),
    matrix_indel_signatures = c("mis1", "mis2"),
    hrdetect = c("hrd1", "hrd2"),
    oncotable = c("ot1", "ot2")
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
  
  # Test data.table with alternative column names
  alt_names_dt <- data.table(
    patient_id = c("sample1", "sample2"),
    gridss_somatic = c("sv1", "sv2"),
    sage_somatic_vcf = c("snv1", "snv2")
  )
  expect_silent(cohort <- Cohort$new(alt_names_dt))
  expect_equal(ncol(cohort$inputs), 3)
  expect_true(all(c("pair", "structural_variants", "somatic_snvs") %in% names(cohort$inputs)))
  
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
