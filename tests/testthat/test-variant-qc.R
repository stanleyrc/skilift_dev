suppressWarnings(devtools::load_all())


library(testthat)

# test <- function() { testthat::test_file("tests/testthat/test-variant-qc.R") }

setup({
    # Helper function to create mock VCF for testing
    create_mock_vcf <<- function(
        num_variants = 3,
        include_normal = TRUE,
        all_pass = TRUE,
        reference = "hg19"
    ) {
        reference_line = sprintf("##reference=%s", reference)
        header <- c(
            "##fileformat=VCFv4.2",
            reference_line,
            "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic variant\">",
            "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">",
            "##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"Allele frequency\">",
            "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">",
            "##FORMAT=<ID=ABQ,Number=1,Type=Float,Description=\"Average base quality\">",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNORMAL\tTUMOR"
        )
        
        if (!include_normal) {
            header[length(header)] <- "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTUMOR"
        }
        
        variants <- character(num_variants)
        for(i in 1:num_variants) {
            filter <- if(all_pass) "PASS" else sample(c("PASS", "LowQual"), 1)
            qual <- sample(50:100, 1)
            
            # Create mock format fields
            tumor_dp <- sample(50:200, 1)
            tumor_af <- runif(1, 0.1, 0.5)
            tumor_ad <- c(round(tumor_dp * (1-tumor_af)), round(tumor_dp * tumor_af))
            tumor_abq <- sample(20:40, 1)
            
            if (include_normal) {
                normal_dp <- sample(30:100, 1)
                normal_af <- runif(1, 0, 0.1)
                normal_ad <- c(round(normal_dp * (1-normal_af)), round(normal_dp * normal_af))
                format_str <- sprintf(
                    "DP:AF:AD:ABQ\t%d:%.3f:%d,%d:30\t%d:%.3f:%d,%d:%d",
                    normal_dp, normal_af, normal_ad[1], normal_ad[2],
                    tumor_dp, tumor_af, tumor_ad[1], tumor_ad[2], tumor_abq
                )
            } else {
                format_str <- sprintf(
                    "DP:AF:AD:ABQ\t%d:%.3f:%d,%d:%d",
                    tumor_dp, tumor_af, tumor_ad[1], tumor_ad[2], tumor_abq
                )
            }
            
            variants[i] <- sprintf(
                "chr%d\t%d\t.\t%s\t%s\t%d\t%s\tSOMATIC\tDP:AF:AD:ABQ\t%s",
                sample(1:22, 1),
                i * 1000,
                sample(c("A", "C", "G", "T"), 1),
                sample(c("A", "C", "G", "T"), 1),
                qual,
                filter,
                format_str
            )
        }
        
        # Write to temporary VCF file
        temp_vcf <- tempfile(fileext = ".vcf")
        writeLines(c(header, variants), temp_vcf)
        return(temp_vcf)
    }
})

# integration test (only works on NYU)
will_run_integrations = FALSE
if (will_run_integrations) {
    test_that("create_oncotable works on real cohort", {
        # Load real clinical pairs
        clinical_pairs_path = "~/projects/Clinical_NYU/db/pairs.rds"
        clinical_pairs = readRDS(clinical_pairs_path)
        cohort <- suppressWarnings(Cohort$new(clinical_pairs[14:15]))

        # Create temp directory for output
        temp_dir <- tempdir()

        # lift_variant_qc
        lift_variant_qc(cohort, temp_dir, cores = 2)

        expect_true(file.exists(file.path(temp_dir, cohort$inputs$pair[1], "sage.qc.json")))
        expect_true(file.exists(file.path(temp_dir, cohort$inputs$pair[2], "sage.qc.json")))

        # Cleanup
        unlink(temp_dir, recursive = TRUE)
    })
}

# Tests for lift_variant_qc
test_that("lift_variant_qc handles basic case correctly", {
    # Create mock Cohort object
    mock_inputs <- data.table(
        pair = c("sample1", "sample2"),
        somatic_snvs = c(
            create_mock_vcf(num_variants = 3),
            create_mock_vcf(num_variants = 4)
        )
    )
    mock_cohort <- structure(
        list(inputs = mock_inputs),
        class = "Cohort"
    )
    
    # Create temp directory for output
    temp_dir <- tempfile()
    
    # Run function
    lift_variant_qc(mock_cohort, temp_dir, cores = 1)
    
    # Check outputs
    expect_true(dir.exists(temp_dir))
    expect_true(dir.exists(file.path(temp_dir, "sample1")))
    expect_true(dir.exists(file.path(temp_dir, "sample2")))
    expect_true(file.exists(file.path(temp_dir, "sample1", "sage.qc.json")))
    expect_true(file.exists(file.path(temp_dir, "sample2", "sage.qc.json")))
    
    # Check content of JSON files
    json1 <- fromJSON(file.path(temp_dir, "sample1", "sage.qc.json"))
    json2 <- fromJSON(file.path(temp_dir, "sample2", "sage.qc.json"))
    
    expect_equal(nrow(json1), 3)
    expect_equal(nrow(json2), 4)
    
    # Clean up
    unlink(temp_dir, recursive = TRUE)
    unlink(mock_inputs$somatic_snvs)
})

test_that("lift_variant_qc handles invalid inputs", {
    # Test with non-Cohort object
    expect_error(
        lift_variant_qc(list(), tempdir()),
        "Input must be a Cohort object"
    )
    
    # Test with missing required column
    mock_inputs <- data.table(
        pair = c("sample1", "sample2")
        # Missing somatic_snvs column
    )
    mock_cohort <- structure(
        list(inputs = mock_inputs),
        class = "Cohort"
    )
    expect_error(
        lift_variant_qc(mock_cohort, tempdir()),
        "Missing required column 'somatic_snvs' in cohort"
    )
})

test_that("lift_variant_qc handles errors in individual samples", {
    # Create mock Cohort with one valid and one invalid VCF
    mock_inputs <- data.table(
        pair = c("sample1", "sample2"),
        somatic_snvs = c(
            create_mock_vcf(num_variants = 3),
            "nonexistent.vcf"
        )
    )
    mock_cohort <- structure(
        list(inputs = mock_inputs),
        class = "Cohort"
    )
    
    # Create temp directory for output
    temp_dir <- tempfile()
    
    # Should complete with warning but not error
    expect_warning(
        lift_variant_qc(mock_cohort, temp_dir, cores = 1),
        "Error processing sample2"
    )
    
    # Check that valid sample was processed
    expect_true(file.exists(file.path(temp_dir, "sample1", "sage.qc.json")))
    
    # Clean up
    unlink(temp_dir, recursive = TRUE)
    unlink(mock_inputs$somatic_snvs[1])
})

test_that("lift_variant_qc works with parallel processing", {
    # Create mock Cohort with multiple samples
    vcf_files <- replicate(4, create_mock_vcf(num_variants = 3))
    mock_inputs <- data.table(
        pair = paste0("sample", 1:4),
        somatic_snvs = vcf_files
    )
    mock_cohort <- structure(
        list(inputs = mock_inputs),
        class = "Cohort"
    )
    
    # Create temp directory for output
    temp_dir <- tempfile()
    
    # Run with multiple cores
    lift_variant_qc(mock_cohort, temp_dir, cores = 2)
    
    # Check all outputs exist
    for(i in 1:4) {
        expect_true(
            file.exists(file.path(temp_dir, paste0("sample", i), "sage.qc.json"))
        )
    }
    
    # Clean up
    unlink(temp_dir, recursive = TRUE)
    unlink(vcf_files)
})

# Unit tests for create_variant_qc
test_that("create_variant_qc processes paired tumor-normal VCF correctly", {
    # Create mock paired VCF
    vcf_file <- create_mock_vcf(num_variants = 5, include_normal = TRUE)
    
    # Run function
    result <- create_variant_qc(vcf_file, "hg19")
    
    # Check structure
    expect_s3_class(result, "data.table")
    expect_equal(nrow(result), 5)
    
    # Check columns for paired analysis
    expected_cols <- c(
        "chromosome", "position", "reference", "alternate", "filter",
        "mapping_quality", "tumor_depth", "normal_depth", "normal_alt_counts",
        "tumor_alt_counts", "tumor_abq", "tumor_vaf", "normal_vaf", "id"
    )
    expect_true(all(expected_cols %in% names(result)))
    
    # Check data types
    expect_type(result$chromosome, "character")
    expect_type(result$tumor_depth, "double")
    expect_type(result$normal_depth, "double")
    expect_type(result$tumor_vaf, "double")
    
    # Check value ranges
    expect_true(all(result$tumor_vaf >= 0 & result$tumor_vaf <= 1))
    expect_true(all(result$normal_vaf >= 0 & result$normal_vaf <= 1))
    expect_true(all(result$filter == "PASS"))
    
    # Clean up
    unlink(vcf_file)
})

test_that("create_variant_qc processes tumor-only VCF correctly", {
    # Create mock tumor-only VCF
    vcf_file <- create_mock_vcf(num_variants = 5, include_normal = FALSE)
    
    # Run function
    result <- create_variant_qc(vcf_file, "hg19")
    
    # Check structure
    expect_s3_class(result, "data.table")
    expect_equal(nrow(result), 5)
    
    # Check columns for tumor-only analysis
    expected_cols <- c(
        "chromosome", "position", "reference", "alternate", "filter",
        "mapping_quality", "tumor_depth", "tumor_alt_counts", 
        "tumor_abq", "tumor_vaf", "id"
    )
    expect_true(all(expected_cols %in% names(result)))
    
    # Verify normal-specific columns are not present
    expect_false("normal_depth" %in% names(result))
    expect_false("normal_vaf" %in% names(result))
    
    # Clean up
    unlink(vcf_file)
})

test_that("create_variant_qc handles different reference genomes", {
    # Test with hg38
    vcf_file <- create_mock_vcf(num_variants = 3, reference = "hg38")
    result <- create_variant_qc(vcf_file, "hg38")
    expect_s3_class(result, "data.table")
    unlink(vcf_file)
    
    # Test with default hg19
    vcf_file <- create_mock_vcf(num_variants = 3, reference = "hg19")
    result <- create_variant_qc(vcf_file)
    expect_s3_class(result, "data.table")
    unlink(vcf_file)
})

test_that("create_variant_qc handles VCF with non-PASS variants", {
    # Create mock VCF with mix of PASS and non-PASS variants
    vcf_file <- create_mock_vcf(num_variants = 10, all_pass = FALSE)
    
    # Run function
    result <- create_variant_qc(vcf_file, "hg19")
    
    # Check that only PASS variants are included
    expect_true(all(result$filter == "PASS"))
    expect_true(nrow(result) <= 10)  # Should be less if some variants were filtered
    
    # Clean up
    unlink(vcf_file)
})

test_that("create_variant_qc handles invalid inputs gracefully", {
    # Test with non-existent file
    expect_error(create_variant_qc("nonexistent.vcf"))
})


# integration test (only works on NYU)
will_run_integrations = FALSE
if (will_run_integrations) {

test_that("create_oncotable works on real cohort", {
  # Load real clinical pairs
  clinical_pairs_path = "~/projects/Clinical_NYU/db/pairs.rds"
  clinical_pairs = readRDS(clinical_pairs_path)
  cohort <- suppressWarnings(Cohort$new(clinical_pairs[14:15]))

  # Create temp directory for output
  temp_dir <- tempdir()

  # lift_variant_qc
  lift_variant_qc(cohort, temp_dir, cores = 2)

  expect_true(file.exists(file.path(temp_dir, cohort$inputs$pair[1], "sage.qc.json")))
  expect_true(file.exists(file.path(temp_dir, cohort$inputs$pair[2], "sage.qc.json")))

  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})
}
