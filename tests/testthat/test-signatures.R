suppressWarnings(devtools::load_all())
library(testthat)
library(data.table)

# test <- function() { testthat::test_file("tests/testthat/test-signatures.R") }

# Create mock data in setup
setup({
    # Mock mutation catalog data
    sbs_matrix <- data.table(
        MutationType = c("A[C>A]A", "C[T>G]T"),
        Count = c(10, 20)
    )
    assign("sbs_matrix", sbs_matrix, envir = parent.frame())
    
    indel_matrix <- data.table(
        MutationType = c("DEL:C:1:0", "INS:T:2:1"),
        Count = c(15, 25)
    )
    assign("indel_matrix", indel_matrix, envir = parent.frame())
    
    temp_output_dir <- tempdir()
    assign("temp_output_dir", temp_output_dir, envir = parent.frame())
    sbs_matrix_file <- tempfile(fileext = ".txt")
    fwrite(sbs_matrix, sbs_matrix_file)
    assign("sbs_matrix_file", sbs_matrix_file, envir = parent.frame())
    
    indel_matrix_file <- tempfile(fileext = ".txt")
    fwrite(indel_matrix, indel_matrix_file)
    assign("indel_matrix_file", indel_matrix_file, envir = parent.frame())

    # Mock decomposed mutation probabilities data
    decomposed_sbs_probs <- data.table(
        MutationType = c("A[C>A]A", "C[T>G]T"),
        Type = c("C>A", "T>G"),
        SBS1 = c(0.1, 0.2),
        SBS2 = c(0.3, 0.4)
    )
    assign("decomposed_sbs_probs", decomposed_sbs_probs, envir = parent.frame())
    
    decomposed_indel_probs <- data.table(
        MutationType = c("DEL:C:1:0", "INS:T:2:1"),
        Type = c("DEL", "INS"),
        ID1 = c(0.15, 0.25),
        ID2 = c(0.35, 0.45)
    )
    assign("decomposed_indel_probs", decomposed_indel_probs, envir = parent.frame())
    
    decomposed_sbs_file <- tempfile(fileext = ".txt")
    fwrite(decomposed_sbs_probs, decomposed_sbs_file)
    assign("decomposed_sbs_file", decomposed_sbs_file, envir = parent.frame())
    
    decomposed_indel_file <- tempfile(fileext = ".txt")
    fwrite(decomposed_indel_probs, decomposed_indel_file)
    assign("decomposed_indel_file", decomposed_indel_file, envir = parent.frame())

    # Mock cohort data
    mock_inputs <- data.frame(
        pair = c("sample1", "sample2"),
        matrix_sbs_signatures = c(sbs_matrix_file, sbs_matrix_file),
        matrix_indel_signatures = c(indel_matrix_file, indel_matrix_file),
        decomposed_sbs_signatures = c(decomposed_sbs_file, decomposed_sbs_file),
        decomposed_indel_signatures = c(decomposed_indel_file, decomposed_indel_file)
    )
    mock_cohort <- structure(
        list(inputs = mock_inputs),
        class = "Cohort"
    )
    assign("mock_cohort", mock_cohort, envir = parent.frame())
    
})

teardown({
    # Clean up temporary files
    if(exists("decomposed_sbs_file")) unlink(decomposed_sbs_file)
    if(exists("decomposed_indel_file")) unlink(decomposed_indel_file)
    if(exists("sbs_matrix_file")) unlink(sbs_matrix_file)
    if(exists("indel_matrix_file")) unlink(indel_matrix_file)
    if(exists("temp_output_dir")) {
        unlink(file.path(temp_output_dir, "test_output"), recursive = TRUE)
    }
})

test_that("read_decomposed_mutationtype_probabilities correctly processes SBS data", {
    # Test SBS processing
    result <- read_decomposed_mutationtype_probabilities(decomposed_sbs_file, is_indel = FALSE)
    
    # Check structure
    expect_s3_class(result, "data.frame")
    expect_named(result, c("signature", "tnc", "p"))
    
    # Check content
    expect_equal(nrow(result), 4)  # 2 types × 2 signatures
    expect_true(all(result$signature %in% c("SBS1", "SBS2")))
    expect_true(all(result$tnc %in% c("C>A", "T>G")))
    expect_true(all(result$p >= 0 & result$p <= 1))
})

test_that("read_decomposed_mutationtype_probabilities correctly processes indel data", {
    # Test indel processing
    result <- read_decomposed_mutationtype_probabilities(decomposed_indel_file, is_indel = TRUE)
    
    # Check structure
    expect_s3_class(result, "data.frame")
    expect_named(result, c("signature", "insdel", "p"))
    
    # Check content
    expect_equal(nrow(result), 4)  # 2 types × 2 signatures
    expect_true(all(result$signature %in% c("ID1", "ID2")))
    expect_true(all(result$insdel %in% c("DEL", "INS")))
    expect_true(all(result$p >= 0 & result$p <= 1))
})

test_that("read_decomposed_mutationtype_probabilities handles missing files appropriately", {
    # Test with non-existent file
    expect_error(
        read_decomposed_mutationtype_probabilities("nonexistent_file.txt", is_indel = FALSE),
        "File 'nonexistent_file.txt' does not exist or is non-readable."
    )
})

test_that("read_decomposed_mutationtype_probabilities handles malformed input appropriately", {
    # Create malformed data file
    malformed_decomposed_probs <- data.table(
        BadColumn = c("A", "B"),
        Value = c(1, 2)
    )
    malformed_decomposed_file <- tempfile(fileext = ".txt")
    fwrite(malformed_decomposed_probs, malformed_decomposed_file)
    
    # Test with malformed file
    expect_error(
        read_decomposed_mutationtype_probabilities(malformed_decomposed_file, is_indel = FALSE),
        "One or more values in 'measure.vars' is invalid."
    )
    
    # Clean up
    unlink(malformed_decomposed_file)
})


test_that("create_mutations_catalog correctly processes SBS data", {
    result <- create_mutations_catalog(sbs_matrix_file, is_indel = FALSE)
    
    # Check structure
    expect_s3_class(result, "data.frame")
    expect_named(result, "data")
    
    # Check content of the inner data frame
    inner_data <- result$data[[1]]
    expect_equal(nrow(inner_data), 2)
    expect_equal(inner_data$id, 1:2)
    expect_equal(inner_data$tnc, c("A[C>A]A", "C[T>G]T"))
    expect_equal(inner_data$mutations, c(10, 20))
})

test_that("create_mutations_catalog correctly processes indel data", {
    result <- create_mutations_catalog(indel_matrix_file, is_indel = TRUE)
    
    # Check structure
    expect_s3_class(result, "data.frame")
    expect_named(result, "data")
    
    # Check content of the inner data frame
    inner_data <- result$data[[1]]
    expect_equal(nrow(inner_data), 2)
    expect_equal(inner_data$id, 1:2)
    expect_equal(inner_data$insdel, c("DEL:C:1:0", "INS:T:2:1"))
    expect_equal(inner_data$mutations, c(15, 25))
})

test_that("create_mutations_catalog handles missing files appropriately", {
    expect_error(
        create_mutations_catalog("nonexistent_file.txt", is_indel = FALSE),
        "File 'nonexistent_file.txt' does not exist or is non-readable."
    )
})

test_that("create_mutations_catalog handles malformed input appropriately", {
    # Create malformed data file
    malformed_matrix <- data.table(
        BadColumn = c("A", "B")
        # Missing count column
    )
    malformed_file <- tempfile(fileext = ".txt")
    fwrite(malformed_matrix, malformed_file)
    
    # Test with malformed file
    expect_error(
        create_mutations_catalog(malformed_file, is_indel = FALSE),
        "undefined columns selected"
    )
    
    # Clean up
    unlink(malformed_file)
})

test_that("lift_signatures validates input correctly", {
    # Test invalid cohort object
    invalid_cohort <- list(inputs = data.frame())
    expect_error(
        lift_signatures(invalid_cohort, temp_output_dir),
        "Input must be a Cohort object"
    )
    
    # Test missing required columns
    invalid_inputs <- data.frame(pair = c("sample1"))
    invalid_cohort <- structure(list(inputs = invalid_inputs), class = "Cohort")
    expect_error(
        lift_signatures(invalid_cohort, temp_output_dir),
        "Missing required columns in cohort: "
    )
})

test_that("lift_signatures creates output directories correctly", {
    test_output_dir <- file.path(temp_output_dir, "test_output")
    
    # Clean up any existing directory
    if (dir.exists(test_output_dir)) unlink(test_output_dir, recursive = TRUE)
    
    # Test directory creation
    lift_signatures(mock_cohort, test_output_dir)
    expect_true(dir.exists(test_output_dir))
    expect_true(dir.exists(file.path(test_output_dir, "sample1")))
    expect_true(dir.exists(file.path(test_output_dir, "sample2")))
    
    # Clean up
    unlink(test_output_dir, recursive = TRUE)
})

test_that("lift_signatures processes files correctly", {
    test_output_dir <- file.path(temp_output_dir, "test_output")
    
    # Clean up any existing directory
    if (dir.exists(test_output_dir)) unlink(test_output_dir, recursive = TRUE)
    
    # Run the function
    lift_signatures(mock_cohort, test_output_dir)
    
    # Check for expected output files for first sample
    sample_dir <- file.path(test_output_dir, "sample1")
    expect_true(file.exists(file.path(sample_dir, "mutation_catalog.json")))
    expect_true(file.exists(file.path(sample_dir, "id_mutation_catalog.json")))
    expect_true(file.exists(file.path(sample_dir, "sbs_decomposed_prob.json")))
    expect_true(file.exists(file.path(sample_dir, "id_decomposed_prob.json")))
    
    # Verify content of one output file
    catalog_content <- jsonlite::read_json(
        file.path(sample_dir, "mutation_catalog.json")
    )
    expect_type(catalog_content[[1]], "list")
    expect_true("data" %in% names(catalog_content[[1]]))
    
    # Clean up
    unlink(test_output_dir, recursive = TRUE)
})

test_that("lift_signatures handles missing input files gracefully", {
    # Create cohort with non-existent files
    bad_inputs <- data.frame(
        pair = "sample1",
        matrix_sbs_signatures = "nonexistent.txt",
        matrix_indel_signatures = "nonexistent.txt",
        decomposed_sbs_signatures = "nonexistent.txt",
        decomposed_indel_signatures = "nonexistent.txt"
    )
    bad_cohort <- structure(list(inputs = bad_inputs), class = "Cohort")
    
    # Should run without error but produce warning
    expect_warning(
        lift_signatures(bad_cohort, temp_output_dir),
        "Missing input file for sample1: nonexistent.txt"
    )
})

## integration tests (only works on NYU hpc)
will_test_integration = FALSE
if (will_test_integration) {
test_that("lift_signatures works on real cohort", {
    clinical_pairs_path = "~/projects/Clinical_NYU/db/pairs.rds"
    clinical_pairs = readRDS(clinical_pairs_path)
    vip_sample = clinical_pairs[patient_id == "10732386", ]
    cohort = Cohort$new(vip_sample)
    # /gpfs/data/imielinskilab/projects/Clinical_NYU/nf-casereports-VIP/signatures/sigprofilerassignment/somatic/397089___B24-1267_vs_397089___B23-2915/sbs_results/Assignment_Solution/Activities/Decomposed_MutationType_Probabilities.txt
    activity = cohort$inputs$activities_sbs_signatures
    matrix = cohort$inputs$matrix_sbs_signatures
    decomposed = cohort$inputs$decomposed_sbs_signatures
    fread(activity)
    fread(matrix)
    fread(decomposed)

    devtools::load_all()
    head(read_decomposed_mutationtype_probabilities(decomposed))
    head(create_mutations_catalog(matrix))
})

}


## debug
DEBUG = FALSE
if (DEBUG) {

}
