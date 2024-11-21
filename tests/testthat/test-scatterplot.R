suppressWarnings(devtools::load_all())


library(testthat)

# test <- function() { testthat::test_file("tests/testthat/test-scatterplot.R") }

setup({
  local_mock_granges <- function(env = parent.frame()) {
    mock_gr <- GRanges(
      seqnames = "1",
      ranges = IRanges(start = 1, end = 100),
      foreground.X = 0.5,
      color = "#FF0000"
    )
    
    mock_gr_path <- tempfile(fileext = ".rds")
    saveRDS(mock_gr, mock_gr_path)
    
    mock_output_file <- tempfile(fileext = ".arrow")
    
    list(gr_path = mock_gr_path, output_file = mock_output_file)
  }

  local_mock_json <- function(env = parent.frame()) {
    mock_json <- tempfile(fileext = ".json")
    jsonlite::write_json(
      list(coordinates = list(default = "hg19", sets = list(hg19 = list(list(chromosome = "1", color = "#FF0000"))))),
      mock_json
    )
    mock_json
  }

  mock_data <<- local_mock_granges()
  mock_json_path <<- local_mock_json()
  mock_gr_path <<- mock_data$gr_path
  mock_output_file <<- mock_data$output_file
  mock_arrow_table <<- data.table(x = 1, y = 0.5, color = 2062260)
  alt_mock_arrow_table <<- data.table(x = 1, y = 0.5, color = 16711680)
  mock_hetsnps_data <<- data.table(
    seqnames = c("1", "1", "2"), 
    start = c(100, 200, 300),
    end = c(101, 201, 301),
    ref.count.t = c(10, 15, 20),
    alt.count.t = c(20, 5, 10),
    alt.frac.n = c(0.4, 0.3, 0.6)
  )
  # Create GRanges with proper seqlengths
  gr <- GRanges(
    seqnames = mock_hetsnps_data$seqnames,
    ranges = IRanges(
      start = mock_hetsnps_data$start,
      end = mock_hetsnps_data$end
    )
  )
  seqlengths(gr) <- c("1" = 249250621, "2" = 243199373)  # hg19 chromosome lengths
  
  # Convert back to data.table and write to file
  mock_hetsnps_file <<- tempfile(fileext = ".txt")
  fwrite(mock_hetsnps_data, mock_hetsnps_file)
  
  # Create mask with matching seqlengths
  mock_mask <- GRanges(
    seqnames = c("1", "2"),
    ranges = IRanges(start = c(1, 1), end = c(1000, 1000))
  )
  seqlengths(mock_mask) <- c("1" = 249250621, "2" = 243199373)  # hg19 chromosome lengths
  mock_mask_file <<- tempfile(fileext = ".rds")
  saveRDS(mock_mask, mock_mask_file)
})

test_that("color2numeric converts hex color to numeric correctly", {
  expect_equal(color2numeric("#FF0000"), 16711680)
  expect_equal(color2numeric("#00FF00"), 65280)
  expect_equal(color2numeric("#0000FF"), 255)
  expect_equal(color2numeric("#FFFFFF"), 16777215)
  expect_equal(color2numeric("#000000"), 0)
})

test_that("get_ref_metadata retrieves metadata correctly", {
  result <- get_ref_metadata(mock_json_path, "hg19")
  expect_equal(nrow(result), 1)
  expect_equal(result$chromosome, list("1"))
  expect_equal(result$color, list("#FF0000"))
})

test_that("make_allelic_hetsnps processes hetsnps correctly", {
  result <- make_allelic_hetsnps(
    het_pileups = mock_hetsnps_file,
    min_normal_freq = 0.2,
    max_normal_freq = 0.8
  )
  
  expect_s4_class(result, "GRanges")
  
  # Convert to data.table for easier testing
  result_dt <- as.data.table(result)
  
  # Should have 6 rows (2 rows per site - major and minor)
  expect_equal(nrow(result_dt), 6)
  
  # Check allele assignments
  major_counts <- result_dt[allele == "major", count]
  minor_counts <- result_dt[allele == "minor", count]
  expect_true(all(major_counts > minor_counts))
  
  # Test filtering
  expect_false(any(mock_hetsnps_data$alt.frac.n <= 0.2))
  expect_false(any(mock_hetsnps_data$alt.frac.n >= 0.8))
})

test_that("make_allelic_hetsnps handles invalid input", {
  expect_error(
    make_allelic_hetsnps(het_pileups = "nonexistent_file.txt"),
    "Please provide a valid path to a hetsnps file"
  )
  
  expect_error(
    make_allelic_hetsnps(het_pileups = NULL),
    "Please provide a valid path to a hetsnps file"
  )
})

test_that("subsample_hetsnps works correctly", {
  result <- suppressWarnings(subsample_hetsnps(
    het_pileups = mock_hetsnps_file,
    mask = mock_mask_file,
    sample_size = 2
  ))
  
  expect_s4_class(result, "GRanges")
  
  result_dt <- as.data.table(result)
  
  # Check color assignments
  expect_true(all(result_dt[allele == "minor"]$col == "blue"))
  
  # Check subsampling
  unique_sites <- unique(result_dt[, .(seqnames, start, end)])
  expect_lte(nrow(unique_sites), 2)
})

test_that("subsample_hetsnps handles missing mask", {
  expect_warning(
    result <- subsample_hetsnps(
      het_pileups = mock_hetsnps_file,
      mask = NULL,
      sample_size = 2
    ),
    "No mask provided, using default mask"
  )
  
  expect_s4_class(result, "GRanges")
})

test_that("subsample_hetsnps handles invalid input", {
  expect_error(
    subsample_hetsnps(het_pileups = NULL),
    "Please provide a valid path to a hetsnps file"
  )
})


test_that("granges_to_arrow_scatterplot creates an arrow table correctly", {
  result <- granges_to_arrow_scatterplot(
    gr_path = mock_gr_path,
    field = "foreground.X",
    ref = 'hg19',
    cov.color.field = "color",
    bin.width = NA
  )
  
  result_dt <- as.data.table(result)
  expect_equal(result_dt, alt_mock_arrow_table)
  
})

test_that("granges_to_arrow_scatterplot handles missing color field gracefully", {
  result <- granges_to_arrow_scatterplot(
    gr_path = mock_gr_path,
    field = "foreground.X",
    ref = 'hg19',
    cov.color.field = NULL,
  )
  
  result_dt <- as.data.table(result)
  color_chr1_numeric <- 2062260
  expect_equal(result_dt$color, color_chr1_numeric)
  
})

test_that("create_scatterplot_arrow handles path input correctly", {
  plot_metadata <- data.table(
    patient.id = "test_patient",
    source = "coverage.arrow",
    x = list(mock_gr_path),
    ref = "hg19",
    field = "foreground.X",
    overwrite = TRUE
  )
  
  datadir <- tempdir()
  
  create_scatterplot_arrow(plot_metadata, datadir)
  
  expected_output_path <- file.path(datadir, plot_metadata$patient.id, plot_metadata$source)
  expect_true(file.exists(expected_output_path))
  
  arrow_table <- arrow::read_feather(expected_output_path)
  expect_equal(arrow_table, mock_arrow_table)
  
})

test_that("create_scatterplot_arrow handles path not in list correctly", {
  plot_metadata <- data.table(
    patient.id = "test_patient",
    source = "coverage.arrow",
    x = mock_gr_path,
    ref = "hg19",
    field = "foreground.X",
    overwrite = TRUE
  )
  
  datadir <- tempdir()
  
  create_scatterplot_arrow(plot_metadata, datadir)
  
  expected_output_path <- file.path(datadir, plot_metadata$patient.id, plot_metadata$source)
  expect_true(file.exists(expected_output_path))
  
  arrow_table <- arrow::read_feather(expected_output_path)
  expect_equal(arrow_table, mock_arrow_table)
})

test_that("create_scatterplot_arrow handles GRanges input correctly", {
  plot_metadata <- data.table(
    patient.id = "test_patient",
    source = "coverage.arrow",
    x = list(GRanges(
      seqnames = "1",
      ranges = IRanges(start = 1, end = 100),
      foreground.X = 0.5,
      color = "#FF0000"
    )),
    ref = "hg19",
    field = "foreground.X",
    overwrite = TRUE
  )
  
  datadir <- tempdir()
  
  create_scatterplot_arrow(plot_metadata, datadir)
  
  expected_output_path <- file.path(datadir, plot_metadata$patient.id, plot_metadata$source)
  expect_true(file.exists(expected_output_path))
  
  arrow_table <- arrow::read_feather(expected_output_path)
  expect_equal(arrow_table, mock_arrow_table)
})

test_that("create_scatterplot_arrow warns when input path does not exist", {
  plot_metadata <- data.table(
    patient.id = "test_patient",
    source = "coverage.arrow",
    x = list("non_existent_path.rds"),
    ref = "hg19",
    field = "foreground.X",
    overwrite = TRUE
  )
  
  datadir <- tempdir()
  
  expect_warning(
    create_scatterplot_arrow(plot_metadata, datadir),
    "Input coverage file does not exist for name:"
  )
})

test_that("create_scatterplot_arrow handles multiple entries in data.table correctly", {
  plot_metadata <- data.table(
    patient.id = c("test_patient1", "test_patient2"),
    source = c("coverage1.arrow", "coverage2.arrow"),
    x = list(mock_gr_path, mock_gr_path),
    ref = c("hg19", "hg19"),
    field = c("foreground.X", "foreground.X"),
    overwrite = TRUE
  )
  
  datadir <- tempdir()
  
  create_scatterplot_arrow(plot_metadata[1], datadir)
  create_scatterplot_arrow(plot_metadata[2], datadir)
  
  expected_output_path1 <- file.path(datadir, plot_metadata$patient.id[1], plot_metadata$source[1])
  expected_output_path2 <- file.path(datadir, plot_metadata$patient.id[2], plot_metadata$source[2])
  
  expect_true(file.exists(expected_output_path1))
  expect_true(file.exists(expected_output_path2))
  
  arrow_table1 <- arrow::read_feather(expected_output_path1)
  arrow_table2 <- arrow::read_feather(expected_output_path2)
  
  expect_equal(arrow_table1, mock_arrow_table)
  expect_equal(arrow_table2, mock_arrow_table)
})

test_that("lift_denoised_coverage handles various input scenarios", {
  # Create temp directory for output and ensure it's clean
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))
  
  # Create a valid mock coverage file
  mock_gr <- GRanges(
    seqnames = "1",
    ranges = IRanges(start = 1, end = 100),
    foreground.X = 0.5
  )
  mock_coverage_file <- tempfile(fileext = ".rds")
  saveRDS(mock_gr, mock_coverage_file)
  on.exit(unlink(mock_coverage_file), add = TRUE)
  
  # 1. Complete cohort with all required inputs
  complete_inputs <- data.table(
    pair = c("test_patient"),
    tumor_coverage = c(mock_coverage_file)
  )
  complete_cohort <- suppressWarnings(Cohort$new(complete_inputs))
  
  # 2. Cohort with missing inputs
  missing_inputs <- data.table(
    pair = c("test_patient", "missing_patient"),
    tumor_coverage = c(mock_coverage_file, NA)
  )
  missing_cohort <- suppressWarnings(Cohort$new(missing_inputs))
  
  # Test: Non-existent output directory
  non_existent_dir <- file.path(temp_dir, "non_existent")
  unlink(non_existent_dir, recursive = TRUE) # Ensure directory doesn't exist
  expect_false(dir.exists(non_existent_dir))
  lift_denoised_coverage(complete_cohort, non_existent_dir)
  expect_true(dir.exists(non_existent_dir))
  expect_true(dir.exists(file.path(non_existent_dir, "test_patient")))
  expect_true(file.exists(file.path(non_existent_dir, "test_patient", "coverage.arrow")))
  
  # Test: Invalid input (not a Cohort object)
  expect_error(
    lift_denoised_coverage(data.table(), temp_dir),
    "Input must be a Cohort object"
  )
  
  # Test: Missing required columns
  invalid_cohort <- suppressWarnings(Cohort$new(data.table(pair = "test")))
  expect_error(
    lift_denoised_coverage(invalid_cohort, temp_dir),
    "Missing required columns in cohort: tumor_coverage"
  )
  
  # Test: Cohort with missing inputs for some samples
  warning_msg <- capture_warnings(
    lift_denoised_coverage(missing_cohort, temp_dir)
  )
  expect_match(
    warning_msg,
    "Tumor coverage file missing for missing_patient",
    all = FALSE
  )
  
  # Test: Existing output directory
  existing_dir <- file.path(temp_dir, "existing")
  dir.create(existing_dir)
  lift_denoised_coverage(complete_cohort, existing_dir)
  expect_true(file.exists(file.path(existing_dir, "test_patient", "coverage.arrow")))
})

test_that("lift_hetsnps handles various input scenarios", {
  # Create temp directory for output and ensure it's clean
  temp_dir <- tempfile()
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))
  
  # Create a valid mock hetsnps file
  mock_hetsnps_data <- data.table(
    seqnames = c("1", "1", "2"),
    start = c(100, 200, 300),
    end = c(101, 201, 301),
    ref.count.t = c(10, 15, 20),
    alt.count.t = c(20, 5, 10),
    alt.frac.n = c(0.4, 0.3, 0.6)
  )
  mock_hetsnps_file <- tempfile(fileext = ".txt")
  fwrite(mock_hetsnps_data, mock_hetsnps_file)
  on.exit(unlink(mock_hetsnps_file), add = TRUE)
  
  # Create mock mask file
  mock_mask <- GRanges(
    seqnames = c("1", "2"),
    ranges = IRanges(start = c(1, 1), end = c(1000, 1000))
  )
  seqlengths(mock_mask) <- c("1" = 249250621, "2" = 243199373)
  mock_mask_file <- tempfile(fileext = ".rds")
  saveRDS(mock_mask, mock_mask_file)
  on.exit(unlink(mock_mask_file), add = TRUE)
  
  # 1. Complete cohort with all required inputs
  complete_inputs <- data.table(
    pair = c("test_patient"),
    het_pileups = c(mock_hetsnps_file)
  )
  complete_cohort <- suppressWarnings(Cohort$new(complete_inputs))
  
  # 2. Cohort with missing inputs
  missing_inputs <- data.table(
    pair = c("test_patient", "missing_patient"),
    het_pileups = c(mock_hetsnps_file, NA)
  )
  missing_cohort <- suppressWarnings(Cohort$new(missing_inputs))
  
  # Test: Non-existent output directory
  non_existent_dir <- file.path(temp_dir, "non_existent")
  unlink(non_existent_dir, recursive = TRUE) # Ensure directory doesn't exist
  expect_false(dir.exists(non_existent_dir))
  
  # Set environment variable for mask file
  withr::with_envvar(
    new = c("SKILIFT_MASK" = mock_mask_file),
    suppressWarnings(lift_hetsnps(complete_cohort, non_existent_dir))
  )
  
  expect_true(dir.exists(non_existent_dir))
  expect_true(dir.exists(file.path(non_existent_dir, "test_patient")))
  expect_true(file.exists(file.path(non_existent_dir, "test_patient", "hetsnps.arrow")))
  
  # Test: Invalid input (not a Cohort object)
  expect_error(
    lift_hetsnps(data.table(), temp_dir),
    "Input must be a Cohort object"
  )
  
  # Test: Missing required columns
  invalid_cohort <- suppressWarnings(Cohort$new(data.table(pair = "test")))
  expect_error(
    lift_hetsnps(invalid_cohort, temp_dir),
    "Missing required columns in cohort: het_pileups"
  )
  
  # Test: Cohort with missing inputs for some samples
  warning_msg <- capture_warnings(
    lift_hetsnps(missing_cohort, temp_dir)
  )
  expect_match(
    warning_msg,
    "Het pileups file missing for missing_patient",
    all = FALSE
  )
  
  # Test: Existing output directory
  existing_dir <- file.path(temp_dir, "existing")
  dir.create(existing_dir)
  
  withr::with_envvar(
    new = c("SKILIFT_MASK" = mock_mask_file),
    suppressWarnings(lift_hetsnps(complete_cohort, existing_dir))
  )
  
  output_file <- file.path(existing_dir, "test_patient", "hetsnps.arrow")
  expect_true(file.exists(output_file))
  
  # Test: Output file content validation
  if (file.exists(output_file)) {
    arrow_table <- arrow::read_feather(output_file)
    expect_true(all(c("x", "y", "color") %in% names(arrow_table)))
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

  # lift_denoised_coverage
  lift_denoised_coverage(cohort, temp_dir, cores = 2)

  expect_true(file.exists(file.path(temp_dir, cohort$inputs$pair[1], "coverage.arrow")))
  expect_true(file.exists(file.path(temp_dir, cohort$inputs$pair[2], "coverage.arrow")))

  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})

test_that("lift_hetsnps works on real cohort", {
  # Load real clinical pairs
  clinical_pairs_path = "~/projects/Clinical_NYU/db/pairs.rds"
  clinical_pairs = readRDS(clinical_pairs_path)
  cohort <- suppressWarnings(Cohort$new(clinical_pairs[14:15]))
  
  # Create temp directory for output
  temp_dir <- tempdir()
  
  # lift_hetsnps
  suppressWarnings(lift_hetsnps(cohort, temp_dir, cores = 2))
  
  expect_true(file.exists(file.path(temp_dir, cohort$inputs$pair[1], "hetsnps.arrow")))
  expect_true(file.exists(file.path(temp_dir, cohort$inputs$pair[2], "hetsnps.arrow")))
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})
}

