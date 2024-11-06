suppressWarnings(devtools::load_all())


library(testthat)

# test <- function() { testthat::test_file("tests/testthat/test-oncotable.R") }

setup({
  ot_test_paths <<- list(
    oncotable = system.file('extdata/test_data/oncotable_test_data/new_oncotable/oncotable.rds', package='Skilift'),
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
    oncokb_maf = ot_test_paths$oncokb_snvcn_maf
  ))

  expect_equal(result_oncotable, expected_oncotable)
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

