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
    # oncokb_maf = system.file('extdata/test_data/oncotable_test_data/oncokb.maf', package='Skilift'),
    oncokb_snvcn_maf = system.file('extdata/test_data/oncotable_test_data/oncokb_snvcn.maf', package='Skilift'),
    oncokb_cna = system.file('extdata/test_data/oncotable_test_data/oncokb_cna.txt', package='Skilift'),
    oncokb_fusions = system.file('extdata/test_data/oncotable_test_data/oncokb_fusions.tsv', package='Skilift')
  )

  # gencode <<- process_gencode('/gpfs/data/imielinskilab/DB/GENCODE/gencode.v19.annotation.gtf.nochr.rds')
  test_rds_path <- system.file("extdata/test_data/test_gencode_v29lift37.rds", package = "Skilift")
  gencode <<- process_gencode(test_rds_path)

  cytoband_path = system.file("extdata/data/cytoband.rds", package = "Skilift")
  cytoband <<- readRDS(cytoband_path)
})

test_that("process_gencode handles NULL input", {
  expect_error(process_gencode(NULL), "gencode file must be provided")
})

test_that("process_gencode handles .rds input", {
  test_rds_path <- system.file("extdata/test_data/test_gencode_v29lift37.rds", package = "Skilift")
  result <- process_gencode(test_rds_path)
  expect_true(is(result, "GRanges"))
})

test_that("process_cytoband handles NULL input", {
  expect_error(process_cytoband(NULL), "cytoband file must be provided")
})

test_that("process_cytoband handles .rds input", {
  cytoband_path <- system.file("extdata/data/cytoband.rds", package = "Skilift")
  result <- process_cytoband(cytoband_path)
  expect_true(is(result, "GRanges"))
  expect_true(all(c("band", "stain", "chrom_name", "chromband") %in% names(mcols(result))))
})

test_that("process_cytoband handles text file input", {
  # Create a temporary cytoband text file
  temp_file <- tempfile(fileext = ".txt")
  cyto_data <- data.table(
    seqnames = c("chr1", "chr1", "chr1"),
    start = c(0, 2300000, 5400000),
    end = c(2300000, 5400000, 7200000),
    band = c("p36.33", "p36.32", "p36.31"),
    stain = c("gneg", "gpos25", "gneg")
  )
  fwrite(cyto_data, temp_file, col.names = FALSE)
  
  result <- process_cytoband(temp_file)
  expect_true(is(result, "GRanges"))
  expect_true(all(c("band", "stain", "chrom_name", "chromband") %in% names(mcols(result))))
  expect_equal(length(result), 3)
  
  unlink(temp_file)
})

test_that("process_cytoband handles coarse parameter", {
  temp_file <- tempfile(fileext = ".txt")
  cyto_data <- data.table(
    seqnames = c("chr1", "chr1", "chr1"),
    start = c(0, 2300000, 5400000),
    end = c(2300000, 5400000, 7200000),
    band = c("p36.33", "p36.32", "p36.31"),
    stain = c("gneg", "gpos25", "gneg")
  )
  fwrite(cyto_data, temp_file, col.names = FALSE)
  result_fine <- process_cytoband(temp_file, coarse = FALSE)
  result_coarse <- process_cytoband(temp_file, coarse = TRUE)
  
  # Check that coarse bands have no decimal points
  expect_true(any(grepl("\\.", result_fine$chromband)))
  expect_false(any(grepl("\\.", result_coarse$chromband)))
})

test_that("process_cytoband handles zero-based coordinates", {
  # Create a temporary cytoband text file with zero-based coordinates
  temp_file <- tempfile(fileext = ".txt")
  cyto_data <- data.table(
    seqnames = c("chr1", "chr1", "chr1"),
    start = c(0, 2300000, 5400000),
    end = c(2300000, 5400000, 7200000),
    band = c("p36.33", "p36.32", "p36.31"),
    stain = c("gneg", "gpos25", "gneg")
  )
  fwrite(cyto_data, temp_file, col.names = FALSE)
  
  result <- process_cytoband(temp_file)
  expect_true(min(start(result)) > 0)  # Should be 1-based
  expect_equal(start(result)[1], 1)    # First start position should be 1
  
  unlink(temp_file)
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

test_that("collect_oncokb_fusions handles missing file", {
  result <- collect_oncokb_fusions(NULL, gencode, cytoband, verbose = FALSE)  # Add cytoband
  expect_equal(result$vartype, NA)
  expect_equal(result$source, "oncokb_fusions")
  
  result <- collect_oncokb_fusions("nonexistent.tsv", gencode, cytoband, verbose = FALSE)  # Add cytoband
  expect_equal(result$vartype, NA)
  expect_equal(result$source, "oncokb_fusions")
})

test_that("collect_oncokb_fusions handles valid input", {
  result <- collect_oncokb_fusions(
    ot_test_paths$oncokb_fusions,
    gencode,
    cytoband,
    verbose = FALSE
  )
  
  # Check basic structure
  expect_true(is.data.table(result))
  expect_true(nrow(result) > 0)
  
  # Check required columns exist
  expected_cols <- c("gene", "gene_summary", "role", "value", "vartype", "type",
                    "tier", "tier_description", "variant_summary", "therapeutics",
                    "resistances", "diagnoses", "prognoses", "effect",
                    "effect_description", "fusion_genes", "fusion_gene_coords",
                    "track", "source")
  expect_true(all(expected_cols %in% names(result)))
  
  # Check specific values
  expect_equal(result$track[1], "variants")
  expect_equal(result$source[1], "oncokb_fusions")
  expect_true(all(result$tier %in% 1:3))
  expect_true(all(result$tier_description %in% c("Clinically Actionable", "Clinically Significant", "VUS")))
  expect_true(all(result$vartype %in% c("fusion", "outframe_fusion")))
  
  # Check fusion gene coordinates format
  expect_true(all(grepl("^[0-9XY]+:[0-9]+-[0-9]+.,", result$fusion_gene_coords)))
})

test_that("collect_oncokb_fusions handles empty input file", {
  # Create temporary empty TSV file
  temp_file <- tempfile(fileext = ".tsv")
  empty_dt <- data.table(
    FUSION = character(),
    GENE_SUMMARY = character(),
    Role = character(),
    LEVEL_1 = character(),
    LEVEL_2 = character(),
    LEVEL_R1 = character(),
    LEVEL_Dx1 = character(),
    LEVEL_Px1 = character(),
    ONCOGENIC = character(),
    silent = logical(),
    in.frame = logical()
  )
  fwrite(empty_dt, temp_file)
  
  result <- collect_oncokb_fusions(temp_file, gencode, cytoband, verbose = FALSE)  # Add cytoband
  expect_equal(result$vartype, NA)
  expect_equal(result$source, "oncokb_fusions")
  
  unlink(temp_file)
})

test_that("collect_oncokb_fusions handles malformed gene names", {
  # Create temporary TSV file with invalid gene names
  temp_file <- tempfile(fileext = ".tsv")
  malformed_dt <- data.table(
    FUSION = c("INVALID1-INVALID2", "NONEXISTENT1-NONEXISTENT2"),
    Hugo_Symbol = c("INVALID1-INVALID2", "NONEXISTENT1-NONEXISTENT2"),
    GENE_SUMMARY = c("Test summary 1", "Test summary 2"),
    VARIANT_SUMMARY = c("Test summary 1", "Test summary 2"),
    min_cn = c(2, 2),
    Role = c("Oncogene", "TSG"),
    LEVEL_1 = c(NA, NA),
    LEVEL_2 = c(NA, NA),
    LEVEL_R1 = c(NA, NA),
    LEVEL_Dx1 = c(NA, NA),
    LEVEL_Px1 = c(NA, NA),
    ONCOGENIC = c("Oncogenic", "Likely Oncogenic"),
    silent = c(FALSE, FALSE),
    in.frame = c(TRUE, FALSE),
    MUTATION_EFFECT = c("Missense", "Truncating"),
    MUTATION_EFFECT_DESCRIPTION = c("Test effect 1", "Test effect 2")
  )
  fwrite(malformed_dt, temp_file)
  
  result <- collect_oncokb_fusions(
    temp_file,
    gencode,
    cytoband,
    verbose = FALSE
  )
  # Should still return a valid data.table with NA coordinates for invalid genes
  expect_true(is.data.table(result))
  expect_true(all(is.na(result$fusion_gene_coords)))
  
  unlink(temp_file)
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

test_that("create_oncotable handles Cohort objects correctly", {
  # Create test Cohort object
  test_inputs <- data.table(
    pair = c("397089", "397090"),  # Second pair is fake to test error handling
    somatic_variant_annotations = c(
      ot_test_paths$unit_annotated_bcf,
      "non_existent_file.bcf"
    ),
    oncokb_fusions = c(
      ot_test_paths$oncokb_fusions,
      ot_test_paths$oncokb_fusions
    ),
    jabba_gg = c(
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
    oncokb_snv = c(
      ot_test_paths$oncokb_snvcn_maf,
      ot_test_paths$oncokb_snvcn_maf
    ),
    oncokb_cna = c(
      ot_test_paths$oncokb_cna,
      ot_test_paths$oncokb_cna
    )
  )
  test_cohort <- suppressWarnings(Cohort$new(test_inputs))

  # Create temporary directory for output
  temp_dir <- tempdir()
  
  # Run create_oncotable and capture the returned cohort
  updated_cohort <- (create_oncotable(
    cohort = test_cohort,
    amp_thresh_multiplier = 1.5,
    gencode = system.file("extdata/test_data/test_gencode_v29lift37.rds", package = "Skilift"),
    cytoband = system.file("extdata/data/cytoband.rds", package = "Skilift"),  # Add cytoband parameter
    outdir = temp_dir,
    cores = 2
  ))

  # Test that input validation works
  expect_error(
    create_oncotable(data.table(), outdir = temp_dir),
    "Input must be a Cohort object"
  )

  # Test that returned object is a Cohort
  expect_true(inherits(updated_cohort, "Cohort"))
  
  # Test that oncotable column exists in updated cohort
  expect_true("oncotable" %in% names(updated_cohort$inputs))
  
  # Test that output files were created for successful sample
  expect_true(file.exists(file.path(temp_dir, "397089", "oncotable.rds")))
  expect_true(file.exists(file.path(temp_dir, "397089", "oncotable.txt")))
  
  # Test that oncotable path is correctly set in cohort for successful sample
  expect_equal(
    updated_cohort$inputs[pair == "397089", oncotable],
    file.path(temp_dir, "397089", "oncotable.rds")
  )
  
  # Test with empty cohort
  empty_cohort <- suppressWarnings(Cohort$new(data.table(
    pair = character(),
    somatic_variant_annotations = character(),
    jabba_gg = character()
  )))
  expect_error(
    create_oncotable(
      cohort = empty_cohort,
      outdir = temp_dir
    ),
    "No samples found in the cohort"
  )

  # Cleanup
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
  expected_colnames <- c(
    "gene", "gene_summary", "role", "vartype", "type", "Tier",
    "variant_summary", "therapeutics", "resistances", "diagnoses", "prognoses",
    "effect", "effect_description", "fusion_genes", "fusion_gene_coords",
    "Variant_g", "Variant", "estimated_altered_copies", "segment_cn", "ref",
    "alt", "VAF", "Genome_Location", "id", "seqnames", "start", "end", "sample"
  )
  expect_true(all(expected_colnames %in% names(result)))
  
  # Test JSON file creation
  expect_true(file.exists(out_file))
  json_content <- jsonlite::fromJSON(out_file)
  expect_true(is.data.frame(json_content))
  
  # Test content validation
  expect_true(all(result$type %in% c("trunc", "missense", "synonymous", "SCNA", "fusion")))
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
  expect_true(all(result_no_oncokb$total_copies == ""))

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
  
  # Test: Existing output directory
  existing_dir <- file.path(temp_dir, "existing")
  dir.create(existing_dir)
  lift_filtered_events(complete_cohort, existing_dir)
  expect_true(file.exists(file.path(existing_dir, "397089", "filtered.events.json")))
  
  # Test: Output file content validation
  json_content <- jsonlite::fromJSON(
    file.path(existing_dir, "397089", "filtered.events.json")
  )
  expect_true(is.data.frame(json_content))
  expected_colnames <- c(
    "gene", "vartype", "type", "Tier", "effect", "Variant_g", "Variant",
    "estimated_altered_copies", "segment_cn", "ref", "alt", "VAF",
    "Genome_Location", "id", "seqnames", "start", "end", "gene_summary",
    "role", "variant_summary", "effect_description", "therapeutics",
    "resistances", "prognoses", "fusion_genes", "fusion_gene_coords"
  )
  expect_true(all(expected_colnames %in% names(json_content)))
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

# integration test (only works on NYU hpc)
will_run_integrations = FALSE
if (will_run_integrations) {

test_that("create_oncotable works on real cohort", {
  # Load real clinical pairs
  clinical_pairs_path = "~/projects/Clinical_NYU/db/pairs.rds"
  clinical_pairs = readRDS(clinical_pairs_path)
  cohort <- suppressWarnings(Cohort$new(clinical_pairs[13:15], col_mapping = c("oncokb_snv" = "oncokb_maf")))

  # Create temp directory for output
  temp_dir <- tempdir()

  # Run create_oncotable and capture the returned cohort
  updated_cohort <- (create_oncotable(
    cohort = cohort,
    amp_thresh_multiplier = 1.5,
    gencode = system.file("extdata/test_data/test_gencode_v29lift37.rds", package = "Skilift"),
    outdir = temp_dir,
    cores = 15
  ))

  # Test that returned object is a Cohort
  expect_true(inherits(updated_cohort, "Cohort"))

  # Test that oncotable column exists in updated cohort
  expect_true("oncotable" %in% names(updated_cohort$inputs))

  # Test that output files were created for successful sample
  test_pair_id <- "C46E0A12-024F-11EF-9EC7-77C4AC892824"
  expect_true(file.exists(file.path(temp_dir, test_pair_id, "oncotable.rds")))

  # Test that oncotable path is correctly set in cohort for successful sample
  expect_equal(
    updated_cohort$inputs[pair == test_pair_id, oncotable],
    file.path(temp_dir, test_pair_id, "oncotable.rds")
  )

  # Cleanup
  unlink(temp_dir, recursive = TRUE)

})

# test_that("oncotable produces expected output (fail-safe test)", {
#   expected_oncotable <- readRDS(ot_test_paths$oncotable)
#   gencode <- process_gencode('~/DB/GENCODE/gencode.v29lift37.annotation.nochr.rds')
#   result_oncotable <- suppressWarnings(oncotable(
#     pair = "397089",
#     annotated_bcf = ot_test_paths$annotated_bcf,
#     oncokb_fusions = ot_test_paths$fusions,
#     jabba_rds = ot_test_paths$jabba_simple_gg,
#     complex = ot_test_paths$complex,
#     signature_counts = NULL,  # Assuming signature_counts is not available in test paths
#     gencode = gencode,
#     verbose = TRUE
#   ))
#
#   expect_equal(result_oncotable, expected_oncotable)
# })
}

## Debug
DEBUG <- FALSE
if (DEBUG) 
{

# why is create_oncotable failing after merge?
# conclusion: missing min_cn in new oncokb_cna, some typos, wrong gencode file

clinical_pairs_path = "/gpfs/data/imielinskilab/projects/Clinical_NYU/db/pairs.20241119_105207.rds"
clinical_pairs = readRDS(clinical_pairs_path)
cohort <- suppressWarnings(Cohort$new(clinical_pairs[14]))
row <- cohort$inputs[1,]
gc = process_gencode('/gpfs/data/imielinskilab/DB/GENCODE/gencode.v19.annotation.gtf.nochr.rds')
devtools::load_all()
ot <- oncotable(
    pair = row$pair,
    somatic_variant_annotations = row$somatic_variant_annotations,
    oncokb_fusions = clinical_pairs[14]$oncokb_fusions,
    jabba_gg = row$jabba_gg,
    karyograph = row$karyograph,
    events = row$events,
    signature_counts = row$signature_counts,
    oncokb_snv = row$oncokb_snv,
    oncokb_cna = row$oncokb_cna,
    gencode = gc,
    verbose = TRUE,
    amp.thresh = 2,
    filter = "PASS",
    del.thresh = 0.5
)

## update oncotable_test_data with new oncokb_*
file.copy(row$oncokb_snv, ot_test_paths$oncokb_snvcn_maf, overwrite = TRUE)
file.copy(row$oncokb_cna, ot_test_paths$oncokb_cna, overwrite = TRUE)
file.copy(row$oncokb_fusions, ot_test_paths$oncokb_fusions, overwrite = TRUE)

## create a subsetted gencode for testing
# subset gene list (two genes split by "-"): 
# [1] "DCDC2-RGL1"       "DNAH6-NCAM1"      "FOXP2-MDFIC"      "HDGFRP3-ADAMTSL3"
# [5] "MAPK4-DCC"        "MET-MDFIC"        "PIGF-EML4"        "SMAP1-FAM184A"   
# [9] "SPECC1-BCAS3"     "SPTBN1-EML4"      "SPTBN1-PIGF"      "STRN-CLHC1"      
gencode <- process_gencode('/gpfs/data/imielinskilab/DB/GENCODE/gencode.v19.annotation.gtf.nochr.rds')
# subset gencode on gRanges metadata col gene_name
subset_genes = c("DCDC2", "RGL1", "DNAH6", "NCAM1", "FOXP2", "MDFIC", "HDGFRP3", "ADAMTSL3", "MAPK4", "DCC", "MET", "MDFIC", "PIGF", "EML4", "SMAP1", "FAM184A", "SPECC1", "BCAS3", "SPTBN1", "EML4", "SPTBN1", "PIGF", "STRN", "CLHC1")
gencode_subset <- gencode[gencode$gene_name %in% subset_genes]

# save to file
test_rds_path <- system.file("extdata/test_data/test_gencode_v29lift37.rds", package = "Skilift")
saveRDS(gencode_subset, test_rds_path)

## extract genecode gene locations to its own file
gencode <- process_gencode('/gpfs/data/imielinskilab/DB/GENCODE/gencode.v19.annotation.gtf.nochr.rds')
gencode_dt <- as.data.table(gencode)
# extract gene locations by type == "gene"
gene_locations <- gencode_dt[type == "gene", .(gene_name, seqnames, start, end)]
saveRDS(gene_locations, '~/git/skilift/inst/extdata/data/gene_locations.rds')

## extract cytoband to its own file
cytoband_path = "/gpfs/data/imielinskilab/DB/UCSC/hg19.cytoband.txt"
cytoband <- process_cytoband(cytoband_path)
saveRDS(cytoband, '~/git/skilift/inst/extdata/data/cytoband.rds')
}
