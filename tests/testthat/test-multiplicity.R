suppressWarnings(devtools::load_all())
library(testthat)
library(data.table)

# test <- function() { testthat::test_file("tests/testthat/test-multiplicity.R") }

setup({
    # Create mock SNV data.table
    mock_snv_dt <<- data.table(
        seqnames = c("chr1", "chr2", "chr3"),
        start = c(100, 200, 300),
        end = c(100, 200, 300),
        total_copies = c(2, 3, 4),
        VAF = c(0.5, 0.3, 0.4),
        gene = c("GENE1", "GENE2", "GENE3"),
        variant.c = c("c.100A>T", "c.200G>C", "c.300T>A"),
        variant.p = c("p.K100M", "p.R200T", "p.L300P"),
        alt = c(10, 15, 20),
        ref = c(10, 35, 30),
        FILTER = c("PASS", "PASS", "LowQual")
    )
    
    # Mock settings JSON content
    mock_settings_json <<- tempfile(fileext = ".json")
    mock_settings_content <- list(
        coordinates = list(
            sets = list(
                hg19 = data.frame(
                    chromosome = c("chr1", "chr2", "chr3"),
                    startPoint = c(1, 1, 1),
                    endPoint = c(1000, 1000, 1000)
                )
            )
        )
    )
    jsonlite::write_json(mock_settings_content, mock_settings_json)

    # Create mock Cohort class
    MockCohort <<- R6::R6Class(
        "Cohort",
        public = list(
            inputs = NULL,
            reference_name = "hg19",
            initialize = function() {
                self$inputs <- data.table(
                    pair = c("sample1", "sample2"),
                    somatic_snv_cn = c(mock_rds_path, mock_rds_path),
                    germline_snv_cn = c(mock_rds_path, mock_rds_path)
                )
            }
        )
    )
    
    temp_dir <<- tempdir()
    # Create temp output directory for lift_multiplicity tests
    test_output_dir <<- file.path(temp_dir, "test_output")
    dir.create(test_output_dir, recursive = TRUE)

    # Create temp RDS file with mock data
    mock_rds_path <<- file.path(temp_dir, "mock_snv.rds")
    saveRDS(mock_snv_dt, mock_rds_path)
})

teardown({
    if(exists("mock_rds_path")) {
        if(file.exists(mock_rds_path)) {
            file.remove(mock_rds_path)
        }
    }
    if(exists("mock_settings_json")) {
        if(file.exists(mock_settings_json)) {
            file.remove(mock_settings_json)
        }
    }
    if(exists("test_output_dir")) {
        if(dir.exists(test_output_dir)) {
            unlink(test_output_dir, recursive = TRUE)
        }
    }
    rm(list = c("mock_snv_dt", "mock_rds_path", "temp_dir"), envir = .GlobalEnv)
})

test_that("create_multiplicity processes RDS input correctly", {
    result <- create_multiplicity(mock_rds_path)
    
    expect_s3_class(result, "data.table")
    expect_equal(nrow(result), 3)
    expect_true("annotation" %in% names(result))
    expect_error(create_multiplicity(mock_snv_dt), "Input must be a file path to an RDS file")
})

test_that("multiplicity_to_intervals creates correct structure", {
    # First create multiplicity data using existing function
    mult_data <- create_multiplicity(mock_rds_path)
    
    # Test the function
    result <- multiplicity_to_intervals(
        mult_data,
        field = "total_copies",
        settings = mock_settings_json,
        reference_name = "hg19"
    )
    
    # Structure tests
    expect_type(result, "list")
    expect_named(result, c("settings", "intervals", "connections"))
    expect_s3_class(result$intervals, "data.table")
    expect_s3_class(result$connections, "data.table")
    
    # Content tests
    expect_equal(nrow(result$intervals), 3)
    expect_named(result$intervals, c("chromosome", "startPoint", "endPoint", 
                                   "iid", "title", "type", "y"))
    expect_equal(result$intervals$y, c(2, 3, 4))  # matches total_copies from mock data
    expect_true(all(result$intervals$type == "interval"))
})

test_that("multiplicity_to_intervals handles invalid inputs appropriately", {
    mult_data <- create_multiplicity(mock_rds_path)
    
    # Test invalid field name
    expect_error(
        multiplicity_to_intervals(
            mult_data,
            field = "nonexistent_field",
            settings = mock_settings_json
        )
    )
    
    # Test invalid settings file
    expect_error(
        multiplicity_to_intervals(
            mult_data,
            settings = "nonexistent_file.json"
        )
    )
})
test_that("lift_multiplicity handles somatic mutations correctly", {
    cohort <- MockCohort$new()
    
    # Run function
    lift_multiplicity(
        cohort, 
        is_germline = FALSE, 
        output_data_dir = test_output_dir,
        node_metadata = c("gene", "variant.c", "variant.p", "vaf")  # Only use columns that exist in mock data
    )
    
    # Check output files exist
    expect_true(dir.exists(file.path(test_output_dir, "sample1")))
    expect_true(dir.exists(file.path(test_output_dir, "sample2")))
    expect_true(file.exists(file.path(test_output_dir, "sample1", "mutations.json")))
    expect_true(file.exists(file.path(test_output_dir, "sample2", "mutations.json")))
    
    # Check content of output files
    json_content <- jsonlite::read_json(file.path(test_output_dir, "sample1", "mutations.json"))
    expect_type(json_content, "list")
    expect_named(json_content, c("settings", "intervals", "connections"))
})

test_that("lift_multiplicity handles germline mutations correctly", {
    cohort <- MockCohort$new()
    
    # Run function
    lift_multiplicity(
        cohort, 
        is_germline = TRUE, 
        output_data_dir = test_output_dir,
        node_metadata = c("gene", "variant.c", "variant.p", "vaf")  # Only use columns that exist in mock data
    )
    
    # Check output files exist
    expect_true(file.exists(file.path(test_output_dir, "sample1", "germline_mutations.json")))
    expect_true(file.exists(file.path(test_output_dir, "sample2", "germline_mutations.json")))
    
    # Check content of output files
    json_content <- jsonlite::read_json(file.path(test_output_dir, "sample1", "germline_mutations.json"))
    expect_type(json_content, "list")
    expect_named(json_content, c("settings", "intervals", "connections"))
})

test_that("lift_multiplicity validates input correctly", {
    # Test invalid cohort object
    expect_error(lift_multiplicity(list(), output_data_dir = test_output_dir),
                "Input must be a Cohort object")
    
    # Test missing required column
    bad_cohort <- MockCohort$new()
    bad_cohort$inputs$somatic_snv_cn <- NULL
    expect_error(lift_multiplicity(bad_cohort, output_data_dir = test_output_dir),
                "Missing required column in cohort: somatic_snv_cn")
    
    # Test missing reference name
    bad_cohort <- MockCohort$new()
    bad_cohort$reference_name <- NULL
    expect_error(lift_multiplicity(bad_cohort, output_data_dir = test_output_dir),
                "Reference name not found in cohort object")
})

will_test_integration = FALSE
if (will_test_integration) {
test_that("lift_multiplicity works on real cohort", {
    clinical_pairs_path = "~/projects/Clinical_NYU/db/pairs.rds"
    clinical_pairs = readRDS(clinical_pairs_path)
    vip_sample = clinical_pairs[patient_id == "397089", ]
    cohort = Cohort$new(vip_sample, col_mapping = list(pair = "patient_id"))
    cohort$inputs$somatic_snv_cn

    temp_dir = tempdir()
    lift_multiplicity(cohort, output_data_dir = temp_dir, cores = 4)
    expect_true(dir.exists(temp_dir))
    expect_true(file.exists(file.path(temp_dir, "397089", "mutations.json")))
    unlink(temp_dir, recursive = TRUE)
})

}


## debug
DEBUG = FALSE
if (DEBUG) {

}
