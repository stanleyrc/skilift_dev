suppressWarnings(devtools::load_all())

library(testthat)

# test <- function() { testthat::test_file("tests/testthat/test-segment-width-distribution.R") }

setup({
    # Helper function to create mock gGraph object
    create_mock_ggraph <<- function(seqnames = c("1", "2", "X")) {
        # Create mock nodes GRanges
        nodes <- GenomicRanges::GRanges(
            seqnames = seqnames,
            ranges = IRanges::IRanges(
                start = seq(1, length(seqnames) * 1000, by = 1000),
                end = seq(1000, length(seqnames) * 1000, by = 1000)
            ),
            strand = "*"
        )
        
        # Add mock color metadata
        mcols(nodes)$col <- rep("#FF0000", length(seqnames))
        mcols(nodes)$type <- rep("REF", length(seqnames))
        
        # Create mock edges
        edges <- data.table(
            n1 = c(1, 2),
            n2 = c(2, 3),
            n1.side = c(1, 1),
            n2.side = c(0, 0)
        )
        
        # Create gGraph using gG
        gg <- gGnome::gG(
            nodes = nodes,
            edges = edges,
            meta = list(
                gr.colorfield = "type",
                purity = 0.8,
                ploidy = 2.1
            )
        )
        
        # Save to temp file and return path
        temp_file <- tempfile(fileext = ".rds")
        saveRDS(gg, temp_file)
        return(temp_file)
    }
})

test_that("get_segstats handles missing inputs correctly", {
    # Test missing balanced_jabba_gg
    expect_error(
        get_segstats(
            balanced_jabba_gg = NULL,
            tumor_coverage = "some_path.rds"
        ),
        "Please provide a valid path to a non-integer balanced gGraph file"
    )
    
    # Test missing tumor_coverage
    mock_gg_path <- suppressWarnings(create_mock_ggraph())
    expect_error(
        get_segstats(
            balanced_jabba_gg = mock_gg_path,
            tumor_coverage = NULL
        ),
        "Please provide a valid path to a coverage file"
    )
    
    # Clean up
    unlink(mock_gg_path)
})

## need better mocks to run these tests:
# test_that("get_segstats handles valid inputs correctly", {
#     # Create mock balanced jabba ggraph
#     mock_gg_path <- create_mock_ggraph()
#     
#     # Create mock coverage data with more bins and valid coverage values
#     mock_cov <- GenomicRanges::GRanges(
#         seqnames = rep(c("1", "2", "X"), each = 100),  # More bins per chromosome
#         ranges = IRanges::IRanges(
#             start = seq(1, 300000, by = 1000),
#             end = seq(1000, 300000, by = 1000)
#         )
#     )
#     
#     # Add realistic coverage values that will pass the variance checks
#     set.seed(42)  # For reproducibility
#     coverage_values <- rnorm(300, mean = 2, sd = 0.2)  # Generate realistic coverage values
#     mcols(mock_cov)$foreground.X <- coverage_values
#     
#     # Ensure some variance but not too much
#     # Add a few NAs but keep them below the default max.na threshold
#     na_indices <- sample(1:300, 10)
#     mcols(mock_cov)$foreground.X[na_indices] <- NA
#     
#     # Save mock coverage to temp file
#     mock_cov_path <- tempfile(fileext = ".rds")
#     saveRDS(mock_cov, mock_cov_path)
#     
#     clinical_pairs_path = "~/projects/Clinical_NYU/db/pairs.rds"
#     clinical_pairs = readRDS(clinical_pairs_path)
#     cohort <- suppressWarnings(Cohort$new(clinical_pairs[14]))
#
#     cov = readRDS(cohort$inputs$tumor_coverage[1])
#     mock_gg_path = cohort$inputs$balanced_jabba_gg[1]
#     saveRDS(cov, mock_cov_path)
#     # Test function
#     result <- get_segstats(
#         balanced_jabba_gg = mock_gg_path,
#         tumor_coverage = mock_cov_path,
#         cores = 1
#     )
#     
#     result <- get_segstats(
#         balanced_jabba_gg = system.file("extdata", "test_data", "test.gg.rds", package = "Skilift"),
#         tumor_coverage = system.file("extdata", "test_data", "test.cov.rds", package = "Skilift"),
#         cores = 1,
#         coverage_field = "foreground"
#     )
#     # Check results
#     expect_true(is.data.table(result))
#     expect_true(all(c("seqnames", "start", "end", "raw_mean", "raw_var") %in% names(result)))
#     
#     # Clean up
#     unlink(c(mock_gg_path, mock_cov_path))
# })

# test_that("get_segstats handles NaN values in coverage correctly", {
#     # Create mock balanced jabba ggraph
#     mock_gg_path <- create_mock_ggraph()
#     
#     # Create mock coverage data with NaN
#     mock_cov <- GenomicRanges::GRanges(
#         seqnames = c("1", "2", "X"),
#         ranges = IRanges::IRanges(
#             start = c(1, 1001, 2001),
#             end = c(1000, 2000, 3000)
#         )
#     )
#     mcols(mock_cov)$foreground.X <- c(2.1, NaN, 2.0)
#     
#     # Save mock coverage to temp file
#     mock_cov_path <- tempfile(fileext = ".rds")
#     saveRDS(mock_cov, mock_cov_path)
#     
#     # Test function
#     result <- get_segstats(
#         balanced_jabba_gg = mock_gg_path,
#         tumor_coverage = mock_cov_path,
#         cores = 1
#     )
#     
#     # Check results
#     expect_true(is.data.table(result))
#     expect_true(!any(is.nan(result$raw_mean)))
#     
#     # Clean up
#     unlink(c(mock_gg_path, mock_cov_path))
# })

# integration test (only works on NYU)
will_run_integrations = FALSE
if (will_run_integrations) {

test_that("lift_segment_width_distribution works on real cohort", {
  # Load real clinical pairs
  clinical_pairs_path = "~/projects/Clinical_NYU/db/pairs.rds"
  clinical_pairs = readRDS(clinical_pairs_path)
  cohort <- suppressWarnings(Cohort$new(clinical_pairs[14:15]))
  
  # ggraph = readRDS(cohort$inputs$balanced_jabba_gg[1])
  

  # Create temp directory for output
  temp_dir <- tempdir()

  # lift_copy_number_graph
  suppressWarnings(lift_segment_width_distribution(cohort, output_data_dir = temp_dir, cores = 2))
  
  expect_true(file.exists(file.path(temp_dir, cohort$inputs$pair[1], "ppfit.json")))
  expect_true(file.exists(file.path(temp_dir, cohort$inputs$pair[2], "ppfit.json")))
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})
}

