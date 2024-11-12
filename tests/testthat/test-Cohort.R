suppressWarnings(devtools::load_all())

library(testthat)
library(data.table)

# test <- function() { testthat::test_file("tests/testthat/test-Cohort.R") }

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
  expect_true(all(c("pair", "tumor_type") %in% names(cohort$inputs)))
  
  # Test data.table with all columns
  complete_dt <- data.table(
    pair = c("sample1", "sample2"),
    tumor_type = c("tumor1", "tumor2"),
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
    tumor_type = c("tumor1", "tumor2"),
    extra_col = c("extra1", "extra2")
  )
  expect_warning(
    cohort <- Cohort$new(mixed_dt),
    "No matching column found for"
  )
  expect_equal(ncol(cohort$inputs), 2)
  
  # Test data.table with alternative column names
  alt_names_dt <- copy(complete_dt)
  alt_names_dt[, structural_variants := NULL]
  alt_names_dt[, gridss_somatic := c("sv1", "sv2")]
  expect_silent(cohort <- Cohort$new(alt_names_dt))
  expect_true("structural_variants" %in% names(cohort$inputs))
  
  # Test custom col_mapping for new column
  custom_dt <- copy(complete_dt)
  custom_dt <- custom_dt[, custom_col := c("custom1", "custom2")]
  custom_mapping <- list(
    structural_variant = c("custom_col")
  )
  expect_silent(cohort <- Cohort$new(custom_dt, col_mapping = custom_mapping))
  expect_true("structural_variant" %in% names(cohort$inputs))
  
  # Test custom col_mapping that overrides default mapping
  override_dt <- copy(complete_dt)
  override_dt[, patient_id := c("patient1", "patient2")]
  override_mapping <- list(
    pair = c("patient_id")  # Changed order from default
  )
  expect_silent(
    cohort <- Cohort$new(override_dt, col_mapping = override_mapping)
  )
  expect_true("pair" %in% names(cohort$inputs))
  expect_equal(cohort$inputs$pair, override_dt$patient_id)

test_that("Cohort validate_inputs correctly identifies missing data", {
  # Setup temp directory for file tests
  temp_dir <- tempdir()
  on.exit(unlink(temp_dir, recursive = TRUE))
  
  # Create some test files
  test_file1 <- file.path(temp_dir, "test1.vcf")
  test_file2 <- file.path(temp_dir, "test2.vcf")
  empty_file <- file.path(temp_dir, "empty.vcf")
  writeLines("content", test_file1)
  writeLines("content", test_file2)
  file.create(empty_file)  # Creates empty file
  
  # Test 1: Data table with missing metadata values
  dt_missing_metadata <- data.table(
    pair = c("sample1", "sample2"),
    tumor_type = c("BRCA", NA),
    structural_variants = c(test_file1, test_file2)
  )
  cohort <- Cohort$new(dt_missing_values)
  missing_data <- cohort$validate_inputs()
  
  expect_true(!is.null(missing_data))
  expect_true(any(missing_data$pair == "sample2" & 
                 missing_data$column == "tumor_type" & 
                 missing_data$reason == "NULL or NA value"))
  
  # Test 2: Data table with missing files
  dt_missing_files <- data.table(
    pair = c("sample1", "sample2"),
    tumor_type = c("BRCA", "LUAD"),
    structural_variants = c(test_file1, "nonexistent.vcf")
  )
  cohort <- Cohort$new(dt_missing_files)
  missing_data <- cohort$validate_inputs()
  
  expect_true(!is.null(missing_data))
  expect_true(any(missing_data$pair == "sample2" & 
                 missing_data$column == "structural_variants" & 
                 missing_data$reason == "File does not exist"))
  
  # Test 3: Data table with mixed missing values and files
  dt_mixed_missing <- data.table(
    pair = c("sample1", "sample2", "sample3"),
    tumor_type = c("BRCA", NA, "LUAD"),
    structural_variants = c(test_file1, "nonexistent.vcf", test_file2)
  )
  cohort <- Cohort$new(dt_mixed_missing)
  missing_data <- cohort$validate_inputs()
  
  expect_true(!is.null(missing_data))
  expect_equal(nrow(missing_data), 3)  # Should find both types of missing data
  expect_true(sum(missing_data$reason == "NULL or NA value") == 2)
  expect_true(sum(missing_data$reason == "File does not exist") == 1)
  
  # Test 4: Data table with no missing values or files
  dt_complete <- data.table(
    pair = c("sample1", "sample2"),
    tumor_type = c("BRCA", "LUAD"),
    structural_variants = c(test_file1, test_file2)
  )
  cohort <- Cohort$new(dt_complete)
  expect_message(
    missing_data <- cohort$validate_inputs(),
    "All inputs are valid - no missing values or files found"
  )
  expect_null(missing_data)
  
  # Test 5: Data table with empty file that exists
  dt_empty_file <- data.table(
    pair = c("sample1"),
    tumor_type = c("BRCA"),
    structural_variants = c(empty_file)
  )
  cohort <- Cohort$new(dt_empty_file)
  expect_message(
    missing_data <- cohort$validate_inputs(),
    "All inputs are valid - no missing values or files found"
  )
  expect_null(missing_data)
  
  # Test 6: Empty data table
  dt_empty <- data.table(pair = character(), tumor_type = character())
  cohort <- Cohort$new(dt_empty)
  expect_error(
    cohort$validate_inputs(),
    "No inputs data available to validate"
  )
})

test_that("Cohort constructor handles pipeline directory inputs correctly", {
  # Create temp directory structure
  base_dir <- tempdir()
  pipeline_dir <- file.path(base_dir, "pipeline_output")
  dir.create(pipeline_dir, recursive = TRUE)
  
  # Helper function to create pipeline report and samplesheet
  setup_metadata <- function(dir) {
    # Create pipeline report
    report_dir <- file.path(dir, "pipeline_info")
    dir.create(report_dir, recursive = TRUE)
    report_content <- c(
      paste0("launchDir: ", dir),
      "input: ./samples.csv"
    )
    writeLines(report_content, file.path(report_dir, "pipeline_report.txt"))
    
    # Create samplesheet
    samplesheet <- data.table(
      patient = c("SAMPLE1", "SAMPLE2"),
      tumor_type = c("BRCA", "LUAD"),
      disease = c("Breast", "Lung"),
      primary_site = c("Breast", "Lung"),
      sex = c("F", "M")
    )
    fwrite(samplesheet, file.path(dir, "samples.csv"))
  }
  
  # Test 1: Empty directory
  empty_dir <- file.path(base_dir, "empty")
  dir.create(empty_dir)
  setup_metadata(empty_dir)
  expect_silent(
    cohort <- suppressWarnings(Cohort$new(empty_dir))
  )
  
  # Test 2: Directory with missing files
  partial_dir <- file.path(base_dir, "partial")
  dir.create(partial_dir, recursive = TRUE)
  setup_metadata(partial_dir)
  dir.create(file.path(partial_dir, "gridss_somatic/SAMPLE1"), recursive = TRUE)
  writeLines("", file.path(partial_dir, "gridss_somatic/SAMPLE1/SAMPLE1.high_confidence_somatic.vcf.bgz"))
  expect_silent(cohort <- Cohort$new(partial_dir))
  expect_true("structural_variants" %in% names(cohort$inputs))
  expect_equal(nrow(cohort$inputs), 2)  # Should have both samples even though only one has data
  
  # Test 3: Directory with all expected files, including hrdetect and fusions
  # Ensure the added test data files are used in the test
  complete_dir <- file.path(base_dir, "complete")
  dir.create(complete_dir, recursive = TRUE)
  setup_metadata(complete_dir)
  
  # Create sample file structure
  sample_dirs <- c("SAMPLE1", "SAMPLE2")
  file_structure <- list(
    "gridss_somatic" = "high_confidence_somatic.vcf.bgz",
    "dryclean_tumor" = "drycleaned.cov.rds",
    "hetpileups" = "sites.txt",
    "jabba" = c("jabba.simple.gg.rds", "karyograph.rds"),
    "events" = "complex.rds",
    "fusions" = "fusions.rds",
    "non_integer_balance" = "balanced.gg.rds",
    "lp_phased_balance" = "balanced.gg.rds",
    "sage/somatic" = "sage.somatic.vcf.gz",
    "snpeff/somatic" = "ann.bcf",
    "snv_multiplicity3" = "est_snv_cn_somatic.rds",
    "signatures/sigprofilerassignment/somatic" = list(
      "sbs_results/Assignment_Solution/Activities" = "Assignment_Solution_Activities.txt",
      "indel_results/Assignment_Solution/Activities" = "Assignment_Solution_Activities.txt",
      "sig_inputs/output/SBS" = "sigmat_results.SBS96.all",
      "sig_inputs/output/ID" = "sigmat_results.ID28.all"
    ),
    "hrdetect" = "hrdetect_results.rds"
  )
  
  # Create files for each sample
  for (sample in sample_dirs) {
    for (dir_name in names(file_structure)) {
      if (is.list(file_structure[[dir_name]])) {
        for (subdir in names(file_structure[[dir_name]])) {
          full_path <- file.path(complete_dir, dir_name, sample, subdir)
          dir.create(full_path, recursive = TRUE)
          writeLines("", file.path(full_path, file_structure[[dir_name]][[subdir]]))
        }
      } else {
        full_path <- file.path(complete_dir, dir_name, sample)
        dir.create(full_path, recursive = TRUE)
        if (length(file_structure[[dir_name]]) > 1) {
          for (file in file_structure[[dir_name]]) {
            writeLines("", file.path(full_path, file))
          }
        } else {
          writeLines("", file.path(full_path, file_structure[[dir_name]]))
        }
      }
    }
  }
  
  # Copy the added test data files to the complete directory
  file.copy("inst/extdata/test_data/hrdetect_results.rds", file.path(complete_dir, "hrdetect/SAMPLE1/hrdetect_results.rds"))
  file.copy("inst/extdata/test_data/oncotable_test_data/fusions.rds", file.path(complete_dir, "fusions/SAMPLE1/fusions.rds"))
  
  expect_silent(cohort <- Cohort$new(complete_dir))
  expect_true("hrdetect" %in% names(cohort$inputs))
  expect_true("fusions" %in% names(cohort$inputs))
  expect_equal(nrow(cohort$inputs), 2)
  expect_true(all(c("structural_variants", "tumor_coverage", "jabba_gg") %in% names(cohort$inputs)))
  
  # Test 4: Directory with extra files
  extra_dir <- file.path(base_dir, "extra")
  dir.create(extra_dir, recursive = TRUE)
  setup_metadata(extra_dir)
  file.copy(complete_dir, extra_dir, recursive = TRUE)
  dir.create(file.path(extra_dir, "extra_folder"))
  writeLines("", file.path(extra_dir, "extra_folder/extra_file.txt"))
  expect_silent(cohort <- Cohort$new(extra_dir))
  expect_equal(nrow(cohort$inputs), 2)
  
  # Test 5: Directory with duplicate column mappings
  duplicate_dir <- file.path(base_dir, "duplicate")
  dir.create(duplicate_dir, recursive = TRUE)
  setup_metadata(duplicate_dir)
  
  # Create two different files that map to structural_variants
  dir.create(file.path(duplicate_dir, "gridss_somatic/SAMPLE1"), recursive = TRUE)
  writeLines("", file.path(duplicate_dir, "gridss_somatic/SAMPLE1/SAMPLE1.high_confidence_somatic.vcf.bgz"))
  dir.create(file.path(duplicate_dir, "svaba_sv/SAMPLE1"), recursive = TRUE)
  writeLines("", file.path(duplicate_dir, "svaba_sv/SAMPLE1/SAMPLE1.sv.vcf"))
  
  expect_silent(cohort <- Cohort$new(duplicate_dir))
  expect_equal(nrow(cohort$inputs), 2)
  expect_true("structural_variants" %in% names(cohort$inputs))
  
  # Cleanup
  unlink(base_dir, recursive = TRUE)
})

# integration tests (only works on NYU hpc)
test_that("Cohort constructor handles real pairs table correctly", {
  clinical_pairs_path = "~/projects/Clinical_NYU/db/pairs.rds"
  clinical_pairs = readRDS(clinical_pairs_path)
  cohort <- Cohort$new(clinical_pairs)
  expect_silent(cohort <- Cohort$new(clinical_pairs))
})

test_that("Cohort construct handles real pipeline directory correctly", {
  pipeline_dir = "/gpfs/data/imielinskilab/projects/Clinical_NYU/nf-casereports-other/"
  cohort <- Cohort$new(pipeline_dir)
  expect_silent(cohort <- Cohort$new(pipeline_dir))
})
