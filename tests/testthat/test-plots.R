suppressWarnings(devtools::load_all())


library(testthat)

# test <- function() { testthat::test_file("tests/testthat/test-plots.R") }

setup({
  ot_test_paths <<- list(
    oncotable = system.file('extdata/test_data/oncotable_test_data/new_oncotable/oncotable.rds', package='Skilift'),
    unit_oncotable = system.file('extdata/test_data/oncotable_test_data/new_oncotable/unit_oncotable.rds', package='Skilift'),
    jabba_simple_gg = system.file('extdata/test_data/oncotable_test_data/jabba.simple.gg.rds', package='Skilift')
  )
})

test_that("filtered_events_json creates correct output", {
  # Create temp directory for output
  temp_dir <- tempdir()
  out_file <- file.path(temp_dir, "filtered_events.json")
  
  # Test with return_table = TRUE to check data structure
  result <- filtered_events_json(
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
  
  # Test with different input combinations
  # Test without return_table
  result_no_return <- filtered_events_json(
    pair = "test_sample",
    oncotable = ot_test_paths$unit_oncotable,
    jabba_gg = ot_test_paths$jabba_simple_gg,
    out_file = file.path(temp_dir, "filtered_events_no_return.json"),
    return_table = FALSE
  )
  expect_null(result_no_return)
  
  # Test error handling
  expect_error(
    suppressWarnings(filtered_events_json(
      pair = "test_sample",
      oncotable = "nonexistent_file.rds",
      jabba_gg = ot_test_paths$jabba_simple_gg,
      out_file = file.path(temp_dir, "should_not_exist.json")
    )),
    "cannot open the connection"
  )

  unlink(temp_dir, recursive = TRUE)
})

