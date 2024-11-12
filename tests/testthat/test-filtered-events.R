suppressWarnings(devtools::load_all())


library(testthat)

# test <- function() { testthat::test_file("tests/testthat/test-filtered-events.R") }

setup({
  ot_test_paths <<- list(
    oncotable = system.file('extdata/test_data/oncotable_test_data/new_oncotable/oncotable.rds', package='Skilift'),
    no_oncokb_unit_oncotable = system.file('extdata/test_data/oncotable_test_data/no_oncokb_unit_oncotable.rds', package='Skilift'),
    unit_oncotable = system.file('extdata/test_data/oncotable_test_data/new_oncotable/unit_oncotable.rds', package='Skilift'),
    annotated_bcf = system.file('extdata/test_data/oncotable_test_data/annotated.bcf', package='Skilift'),
    unit_annotated_bcf = system.file('extdata/test_data/oncotable_test_data/unit_annotated.bcf', package='Skilift'),
    jabba_simple_gg = system.file('extdata/test_data/oncotable_test_data/jabba.simple.gg.rds', package='Skilift'),
    complex = system.file('extdata/test_data/oncotable_test_data/complex.rds', package='Skilift'),
    fusions = system.file('extdata/test_data/oncotable_test_data/fusions.rds', package='Skilift'),
    karyograph = system.file('extdata/test_data/oncotable_test_data/karyograph.rds', package='Skilift'),
    oncokb_maf = system.file('extdata/test_data/oncotable_test_data/oncokb.maf', package='Skilift'),
    oncokb_snvcn_maf = system.file('extdata/test_data/oncotable_test_data/oncokb_snvcn.maf', package='Skilift'),
    oncokb_cna = system.file('extdata/test_data/oncotable_test_data/oncokb_cna.txt', package='Skilift')
  )

  # gencode <<- process_gencode('~/DB/GENCODE/gencode.v29lift37.annotation.nochr.rds')
  test_rds_path <- system.file("extdata/test_data/test_gencode_v29lift37.rds", package = "Skilift")
  gencode <<- process_gencode(test_rds_path)
})

test_that("process_gencode handles NULL input", {
  expect_error(process_gencode(NULL), "gencode file must be provided")
})

test_that("process_gencode handles .rds input", {
  test_rds_path <- system.file("extdata/test_data/test_gencode_v29lift37.rds", package = "Skilift")
  result <- process_gencode(test_rds_path)
  expect_true(is(result, "GRanges"))
})

test_that("collect_gene_fusions handles valid input", {
  fusions_path <- ot_test_paths$fusions
  result_fusions <- collect_gene_fusions(fusions_path, gencode, verbose = FALSE)
  expect_true(nrow(result_fusions) > 0)
  expect_true(all(c("gene", "vartype", "fusion_genes", "track", "type", "source", "fusion_gene_coords") %in% colnames(result_fusions)))
})

test_that("collect_complex_events handles valid input", {
  complex_path <- ot_test_paths$complex
  result_complex <- collect_complex_events(complex_path, verbose = FALSE)
  expect_true(nrow(result_complex) > 0)
  expect_true(all(c("value", "type", "track", "source") %in% colnames(result_complex)))
})

test_that("collect_copy_number_jabba handles valid input", {
  result_jabba <- suppressWarnings(collect_copy_number_jabba(
    ot_test_paths$jabba_simple_gg,
    gencode,
    amp.thresh = 4,
    del.thresh = 0.5,
    verbose = FALSE,
    karyograph = ot_test_paths$karyograph
  ))
  expect_true(nrow(result_jabba) > 0)
  expect_true(all(c("value", "type", "track") %in% colnames(result_jabba)))
})

# not testing this function as signatures needs to be updated to use sigprofiler
# test_that("collect_signatures handles valid input", {
#   signature_counts_path <- system.file('extdata/test_data/oncotable_test_data/signature_counts.txt', package='Skilift')
#   result_signatures <- collect_signatures(signature_counts_path, verbose = FALSE)
#   expect_true(nrow(result_signatures) > 0)
#   expect_true(all(c("value", "type", "track", "source") %in% colnames(result_signatures)))
# })

test_that("collect_gene_mutations handles valid input", {
  annotated_bcf_path <- ot_test_paths$unit_annotated_bcf
  result_mutations <- collect_gene_mutations(annotated_bcf_path, ot_test_paths$jabba_simple_gg, filter = 'PASS', verbose = FALSE)
  expect_true(nrow(result_mutations) > 0)
  expect_true(all(c("value", "type", "track", "source") %in% colnames(result_mutations)))
})

test_that("parse_oncokb_tier correctly assigns tiers", {
  # Create test data
  test_oncokb <- data.table(
    LEVEL_1 = c("drug1,drug2", NA, NA, NA),
    LEVEL_2 = c("drug3,drug4", "drug1,drug2", NA, NA),
    LEVEL_R1 = c(NA, NA, "LevelR1", NA),
    LEVEL_Dx1 = c(NA, NA, NA, NA),
    LEVEL_Px1 = c(NA, NA, NA, NA),
    ONCOGENIC = c("Oncogenic", "Likely Oncogenic", "Unknown", "Likely Neutral")
  )
  
  result <- parse_oncokb_tier(test_oncokb)
  
  # Check tier assignments
  expect_equal(result$tier, c(1, 1, 1, 3))
  expect_equal(result$tier_factor, factor(c("Clinically Actionable", "Clinically Actionable", 
                                          "Clinically Actionable", "VUS"),
                                        levels = c("Clinically Actionable", "Clinically Significant", "VUS")))
  
  # Check string concatenation
  expect_equal(result$tx_string, c("drug1,drug2,drug3,drug4", "drug1,drug2", NA, NA))
  expect_equal(result$rx_string, c(NA, NA, "LevelR1", NA))
  expect_true(all(is.na(result$dx_string)))
  expect_true(all(is.na(result$px_string)))
})

test_that("collect_oncokb handles missing file", {
  result <- collect_oncokb(NULL, verbose = FALSE)
  expect_equal(result$type, NA)
  expect_equal(result$source, "oncokb_maf")
})

test_that("collect_oncokb handles valid input", {
  result <- collect_oncokb(ot_test_paths$oncokb_snvcn_maf, verbose = FALSE)
  
  # Check basic structure
  expect_true(is.data.table(result))
  expect_true(nrow(result) > 0)
  
  # Check required columns exist
  expected_cols <- c("gene", "variant.g", "variant.c", "variant.p", "annotation", "type", "tier", "tier_description", "therapeutics", "resistances", "diagnoses", "prognoses", "distance", "major.count", "minor.count", "major_snv_copies", "minor_snv_copies", "total_copies", "VAF", "track", "source")
  expect_true(all(expected_cols %in% names(result)))
  
  # Check values
  expect_equal(result$track[1], "variants")
  expect_true(all(result$tier %in% 1:3))
  expect_true(all(result$tier_description %in% c("Clinically Actionable", "Clinically Significant", "VUS")))
})

test_that("collect_oncokb_cna handles missing file", {
  result <- collect_oncokb_cna(NULL, verbose = FALSE)
  expect_equal(result$type, NA)
  expect_equal(result$source, "oncokb_cna")
})

test_that("collect_oncokb_cna handles valid input", {
  result <- collect_oncokb_cna(ot_test_paths$oncokb_cna, verbose = FALSE)
  
  # Check basic structure
  expect_true(is.data.table(result))
  expect_true(nrow(result) > 0)
  
  # Check required columns exist
  expected_cols <- c("gene", "value", "type", "tier", "tier_description", 
                    "therapeutics", "resistances", "diagnoses", "prognoses",
                    "track", "source")
  expect_true(all(expected_cols %in% names(result)))
  
  # Check values
  expect_equal(result$track[1], "scna")
  expect_equal(result$source[1], "oncokb_cna")
  expect_true(all(result$tier %in% 1:3))
  expect_true(all(result$tier_description %in% c("Clinically Actionable", "Clinically Significant", "VUS")))
  expect_true(all(result$type %in% c("amp", "homdel", NA)))
})

test_that("oncotable produces expected output", {
  expected_oncotable <- readRDS(ot_test_paths$unit_oncotable)
  result_oncotable <- suppressWarnings(oncotable(
    pair = "397089",
    annotated_bcf = ot_test_paths$unit_annotated_bcf,
    fusions = ot_test_paths$fusions,
    jabba_rds = ot_test_paths$jabba_simple_gg,
    complex = ot_test_paths$complex,
    signature_counts = NULL,  # Assuming signature_counts is not available in test paths
    gencode = gencode,
    verbose = TRUE,
    karyograph = ot_test_paths$karyograph,
    oncokb_maf = ot_test_paths$oncokb_snvcn_maf,
    oncokb_cna = ot_test_paths$oncokb_cna
  ))

  # saveRDS(result_oncotable, ot_test_paths$unit_oncotable)
  expect_equal(result_oncotable, expected_oncotable)
})

test_that("create_oncotable handles multiple samples correctly", {
  # Create test cohort data.table
  test_cohort <- data.table(
    pair = c("397089", "397090"),  # Second pair is fake to test error handling
    annotated_bcf = c(
      ot_test_paths$unit_annotated_bcf,
      "non_existent_file.bcf"
    ),
    fusions = c(
      ot_test_paths$fusions,
      ot_test_paths$fusions
    ),
    jabba_simple = c(
      ot_test_paths$jabba_simple_gg,
      ot_test_paths$jabba_simple_gg
    ),
    karyograph = c(
      ot_test_paths$karyograph,
      ot_test_paths$karyograph
    ),
    events = c(
      ot_test_paths$complex,
      ot_test_paths$complex
    ),
    signature_counts = c(NA, NA),  # Optional
    oncokb_maf = c(
      ot_test_paths$oncokb_snvcn_maf,
      ot_test_paths$oncokb_snvcn_maf
    ),
    oncokb_cna = c(
      ot_test_paths$oncokb_cna,
      ot_test_paths$oncokb_cna
    )
  )

  # Create temporary directory for output
  temp_dir <- tempdir()
  
  # Run create_oncotable
  suppressWarnings(create_oncotable(
    cohort = test_cohort,
    amp_thresh_multiplier = 1.5,
    gencode = system.file("extdata/test_data/test_gencode_v29lift37.rds", package = "Skilift"),
    outdir = temp_dir,
    cores = 1
  ))

  # Check that output files were created for successful sample
  expect_true(file.exists(file.path(temp_dir, "397089", "oncotable.rds")))
  expect_true(file.exists(file.path(temp_dir, "397089", "oncotable.txt")))
  
  # Test that successful result matches expected
  result_oncotable <- readRDS(file.path(temp_dir, "397089", "oncotable.rds"))
  expected_oncotable <- readRDS(ot_test_paths$unit_oncotable)
  expect_equal(result_oncotable, expected_oncotable)
  unlink(temp_dir, recursive = TRUE)
})

test_that("create_filtered_events creates correct output", {
  # Create temp directory for output
  temp_dir <- tempdir()
  out_file <- file.path(temp_dir, "filtered_events.json")
  
  # Test with return_table = TRUE to check data structure
  result <- create_filtered_events(
    pair = "test_sample",
    oncotable = ot_test_paths$unit_oncotable,
    jabba_gg = ot_test_paths$jabba_simple_gg,
    out_file = out_file,
    return_table = TRUE
  )
  
  # Test data.table structure
  expect_true(is.data.table(result))
  expect_true(all(c(
    "gene", "fusion_genes", "id", "vartype", "type", "Variant_g", 
    "Variant", "Genome_Location", "fusion_gene_coords", "Tier",
    "therapeutics", "resistances", "diagnoses", "prognoses", "dosage",
    "seqnames", "start", "end", "sample"
  ) %in% names(result)))
  
  # Test JSON file creation
  expect_true(file.exists(out_file))
  json_content <- jsonlite::fromJSON(out_file)
  expect_true(is.data.frame(json_content))
  
  # Test content validation
  expect_true(all(result$type %in% c("trunc", "missense", "synonymous", "amp", "homdel", "fusion")))
  expect_true(all(!is.na(result$gene)))
  expect_true(all(!is.na(result$type)))
  
  # Test with missing oncokb inputs
  result_no_oncokb <- create_filtered_events(
    pair = "test_sample",
    oncotable = ot_test_paths$no_oncokb_unit_oncotable,
    jabba_gg = ot_test_paths$jabba_simple_gg,
    out_file = file.path(temp_dir, "filtered_events_no_oncokb.json"),
    return_table = TRUE
  )

  expect_true(all(result_no_oncokb$tier == ""))
  expect_true(all(result_no_oncokb$dosage == ""))

  # Test without return_table
  result_no_return <- create_filtered_events(
    pair = "test_sample",
    oncotable = ot_test_paths$unit_oncotable,
    jabba_gg = ot_test_paths$jabba_simple_gg,
    out_file = file.path(temp_dir, "filtered_events_no_return.json"),
    return_table = FALSE
  )
  expect_null(result_no_return)
  
  # Test error handling
  expect_error(
    suppressWarnings(create_filtered_events(
      pair = "test_sample",
      oncotable = "nonexistent_file.rds",
      jabba_gg = ot_test_paths$jabba_simple_gg,
      out_file = file.path(temp_dir, "should_not_exist.json")
    )),
    "cannot open the connection"
  )

  unlink(temp_dir, recursive = TRUE)
})

test_that("lift_filtered_events handles various input scenarios", {
  # Create temp directory for output
  temp_dir <- tempdir()
  
  # Create test Cohort objects
  
  # 1. Complete cohort with all required inputs
  complete_inputs <- data.table(
    pair = c("397089"),
    oncotable = c(ot_test_paths$unit_oncotable),
    jabba_gg = c(ot_test_paths$jabba_simple_gg)
  )
  complete_cohort <- suppressWarnings(Cohort$new(complete_inputs))
  
  # 2. Cohort with missing inputs
  missing_inputs <- data.table(
    pair = c("397089", "397090"),
    oncotable = c(ot_test_paths$unit_oncotable, NA),
    jabba_gg = c(ot_test_paths$jabba_simple_gg, ot_test_paths$jabba_simple_gg)
  )
  missing_cohort <- suppressWarnings(Cohort$new(missing_inputs))
  
  # Test: Non-existent output directory
  non_existent_dir <- file.path(temp_dir, "non_existent")
  expect_false(dir.exists(non_existent_dir))
  lift_filtered_events(complete_cohort, non_existent_dir)
  expect_true(dir.exists(non_existent_dir))
  expect_true(dir.exists(file.path(non_existent_dir, "397089")))
  expect_true(file.exists(file.path(non_existent_dir, "397089", "filtered.events.json")))
  
  # Test: Existing output directory
  existing_dir <- file.path(temp_dir, "existing")
  dir.create(existing_dir)
  lift_filtered_events(complete_cohort, existing_dir)
  expect_true(file.exists(file.path(existing_dir, "397089", "filtered.events.json")))
  
  # Test: Invalid input (not a Cohort object)
  expect_error(
    lift_filtered_events(data.table(), temp_dir),
    "Input must be a Cohort object"
  )
  
  # Test: Missing required columns
  invalid_cohort <- suppressWarnings(Cohort$new(data.table(pair = "test")))
  expect_error(
    lift_filtered_events(invalid_cohort, temp_dir),
    "Missing required columns in cohort: oncotable, jabba_gg"
  )
  
  # Test: Cohort with missing inputs for some samples
  warning_msg <- capture_warnings(
    lift_filtered_events(missing_cohort, temp_dir)
  )
  expect_match(
    warning_msg,
    "Error processing 397090: invalid 'description' argument",
    all = FALSE
  )
  
  # Test: Output file content validation
  json_content <- jsonlite::fromJSON(
    file.path(existing_dir, "397089", "filtered.events.json")
  )
  expect_true(is.data.frame(json_content))
  expect_true(all(c(
    "gene", "fusion_genes", "id", "vartype", "type",
    "Variant_g", "Variant", "Genome_Location", "fusion_gene_coords"
  ) %in% names(json_content)))
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})


# test_that("oncotable produces expected output (fail-safe test)", {
#   expected_oncotable <- readRDS(ot_test_paths$oncotable)
#   gencode <- process_gencode('~/DB/GENCODE/gencode.v29lift37.annotation.nochr.rds')
#   result_oncotable <- suppressWarnings(oncotable(
#     pair = "397089",
#     annotated_bcf = ot_test_paths$annotated_bcf,
#     fusions = ot_test_paths$fusions,
#     jabba_rds = ot_test_paths$jabba_simple_gg,
#     complex = ot_test_paths$complex,
#     signature_counts = NULL,  # Assuming signature_counts is not available in test paths
#     gencode = gencode,
#     verbose = TRUE
#   ))
#
#   expect_equal(result_oncotable, expected_oncotable)
# })

