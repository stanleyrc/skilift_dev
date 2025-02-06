suppressWarnings(devtools::load_all())


library(testthat)

# test <- function() { testthat::test_file("tests/testthat/test-copy-number-graph.R") }

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

test_that("lift_copy_number_graph processes total copy number correctly", {
    # Create mock cohort
    mock_events_path <- create_mock_ggraph()
    mock_inputs <- data.table(
        pair = "TEST001",
        events = mock_events_path
    )
    mock_cohort <- list(
        inputs = mock_inputs,
        reference_name = "hg19"
    )
    class(mock_cohort) <- "Cohort"
    
    # Create temp output directory
    temp_dir <- tempdir()
    
    # Run function
    suppressWarnings(
        lift_copy_number_graph(
            mock_cohort,
            output_data_dir = temp_dir,
            is_allelic = FALSE
        )
    )
    
    # Check output
    expect_true(dir.exists(file.path(temp_dir, "TEST001")))
    expect_true(file.exists(file.path(temp_dir, "TEST001", "complex.json")))
    
    # Clean up
    unlink(c(mock_events_path, temp_dir), recursive = TRUE)
})

test_that("lift_copy_number_graph processes allelic copy number correctly", {
    # Create mock cohort
    mock_allelic_path <- create_mock_ggraph()
    mock_inputs <- data.table(
        pair = "TEST001",
        allelic_jabba_gg = mock_allelic_path
    )
    mock_cohort <- list(
        inputs = mock_inputs,
        reference_name = "hg19"
    )
    class(mock_cohort) <- "Cohort"
    
    # Create temp output directory
    temp_dir <- tempdir()
    
    # Run function
    suppressWarnings(
        lift_copy_number_graph(
            mock_cohort,
            output_data_dir = temp_dir,
            is_allelic = TRUE
        )
    )
    
    # Check output
    expect_true(dir.exists(file.path(temp_dir, "TEST001")))
    expect_true(file.exists(file.path(temp_dir, "TEST001", "allelic.json")))
    
    # Clean up
    unlink(c(mock_allelic_path, temp_dir), recursive = TRUE)
})

test_that("lift_copy_number_graph handles missing required columns", {
    # Create mock cohort without required column
    mock_inputs <- data.table(pair = "TEST001")
    mock_cohort <- list(inputs = mock_inputs)
    class(mock_cohort) <- "Cohort"
    
    # Should error due to missing column
    expect_error(
        lift_copy_number_graph(mock_cohort, output_data_dir = tempdir(), is_allelic = FALSE),
        "Missing required columns in cohort: events"
    )
})

test_that("lift_copy_number_graph handles invalid input type", {
    # Test with invalid input type
    invalid_cohort <- list(inputs = data.table(pair = "TEST001"))
    
    expect_error(
        lift_copy_number_graph(invalid_cohort, tempdir()),
        "Input must be a Cohort object"
    )
})

test_that("lift_copy_number_graph handles missing files", {
    # Create mock cohort with non-existent file
    mock_inputs <- data.table(
        pair = "TEST001",
        events = "nonexistent.rds"
    )
    mock_cohort <- list(inputs = mock_inputs)
    class(mock_cohort) <- "Cohort"
    
    # Should warn about missing file but not error
    expect_warning(
        lift_copy_number_graph(mock_cohort, output_data_dir = tempdir(), is_allelic = FALSE),
        "Copy number graph file missing for TEST001"
    )
})

test_that("lift_copy_number_graph processes multiple samples correctly", {
    # Create mock cohort with multiple samples
    mock_events_path1 <- create_mock_ggraph()
    mock_events_path2 <- create_mock_ggraph()
    mock_inputs <- data.table(
        pair = c("TEST001", "TEST002"),
        events = c(mock_events_path1, mock_events_path2)
    )
    mock_cohort <- list(
        inputs = mock_inputs,
        reference_name = "hg19"
    )
    class(mock_cohort) <- "Cohort"
    
    # Create temp output directory
    temp_dir <- tempdir()
    
    # Run function with multiple cores
    suppressWarnings(
        lift_copy_number_graph(
            mock_cohort,
            output_data_dir = temp_dir,
            is_allelic = FALSE,
            cores = 2
        )
    )
    
    # Check output for both samples
    expect_true(file.exists(file.path(temp_dir, "TEST001", "complex.json")))
    expect_true(file.exists(file.path(temp_dir, "TEST002", "complex.json")))
    
    # Clean up
    unlink(c(mock_events_path1, mock_events_path2, temp_dir), recursive = TRUE)
})


# integration test (only works on NYU)
will_run_integrations = FALSE
if (will_run_integrations) {

test_that("lift_hetsnps works on real cohort", {
  # Load real clinical pairs
  clinical_pairs_path = "~/projects/Clinical_NYU/db/pairs.rds"
  clinical_pairs = readRDS(clinical_pairs_path)
  cohort <- suppressWarnings(Cohort$new(clinical_pairs[14:15]))
  
  ggraph = readRDS(cohort$inputs$events[1])

  # Create temp directory for output
  temp_dir <- tempdir()

  # lift_copy_number_graph
  suppressWarnings(lift_copy_number_graph(cohort, output_data_dir = temp_dir, cores = 2))
  
  expect_true(file.exists(file.path(temp_dir, cohort$inputs$pair[1], "complex.json")))
  expect_true(file.exists(file.path(temp_dir, cohort$inputs$pair[2], "complex.json")))
  
  # Cleanup
  unlink(temp_dir, recursive = TRUE)
})
}

# debug
DEBUG = FALSE
if (DEBUG) 
{
# why is $json not outputting a file
# hypothesis: an issue with the json call
# conclusion: the annotations list had to be unlisted when passed to $json

clinical_pairs_path = "~/projects/Clinical_NYU/db/pairs.rds"
clinical_pairs = readRDS(clinical_pairs_path)
cohort <- suppressWarnings(Cohort$new(clinical_pairs[14]))

ggraph = readRDS(cohort$inputs$events[1])

output_data_dir = tempdir()
out_filename = "complex.json"
cn_column = "events"
row <- cohort$inputs[1]
pair_dir <- file.path(output_data_dir, row$pair)

if (!dir.exists(pair_dir)) {
    dir.create(pair_dir, recursive = TRUE)
}

out_file <- file.path(pair_dir, out_filename)

ggraph_path <- row[[cn_column]]

    message(sprintf("Reading gGraph for %s", row$pair))
    ggraph <- readRDS(ggraph_path)
    
    if (!any(class(ggraph) == "gGraph")) {
        warning(sprintf("Input for %s is not a gGraph object", row$pair))
        return(NULL)
    }
    
    # Check sequence names overlap with reference
    seq_lengths <- gGnome::parse.js.seqlengths(
        Skilift:::default_settings_path,
        js.type = "PGV",
        ref = cohort$reference_name
    )
    
    # Reduce gGraph to only sequences that overlap with reference
    ggraph.reduced <- ggraph[seqnames %in% names(seq_lengths)]
    if (length(ggraph.reduced) == 0) {
        warning(sprintf(
            "No overlap between reference sequences and gGraph sequences for %s",
            row$pair
        ))
        return(NULL)
    }
    
    annotations = list(c("bfb", "chromoplexy", "chromothripsis", "del", "dm", "cpxdm", "dup", "pyrgo", "rigma", "simple", "tic", "tyfonas"))
    # Set parameters for json export
    params <- list(
        filename = out_file,
        verbose = TRUE,
        maxcn = 100,
        nfields = if("col" %in% names(mcols(ggraph$nodes$gr))) "col" else NULL,
        annotations = unlist(annotations)
    )
    
    # Generate and write JSON
    message(sprintf("Writing copy number graph JSON for %s", row$pair))

    gGnome::refresh(ggraph[seqnames %in% names(seq_lengths)])$json(
        filename = out_file,
        verbose = TRUE,
        annotations = unlist(annotations),
        maxcn = 100,
        nfields = if("col" %in% names(mcols(ggraph$nodes$gr))) "col" else NULL,
        # cid.field = field
    )
    do.call(gGnome::refresh(ggraph.reduced)$json, params)
}
