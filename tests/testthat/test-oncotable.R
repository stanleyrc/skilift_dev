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
    karyograph = system.file('extdata/test_data/oncotable_test_data/karyograph.rds', package='Skilift')
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
  jabba_rds_path <- ot_test_paths$jabba_simple_gg
  result_jabba <- suppressWarnings(collect_copy_number_jabba(jabba_rds_path, gencode, amp.thresh = 4, del.thresh = 0.5, verbose = FALSE))
  expect_true(nrow(result_jabba) > 0)
  expect_true(all(c("value", "type", "track") %in% colnames(result_jabba)))
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
    verbose = TRUE
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

