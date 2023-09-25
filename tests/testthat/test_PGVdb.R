library(testthat)

setup({
  library(parallel)
  library(R6)
  library(data.table)
  library(jsonlite)
  library(httr)
  devtools::load_all("../../gGnome/gGnome")
})

devtools::load_all(".")
context("PGVdb")

load_paths <- function() {
  datafiles.json <- system.file("extdata", "pgv", "public", "datafiles.json", package = "PGVdb")
  datadir <- system.file("extdata", "pgv", "public", "data", package = "PGVdb")
  settings <- system.file("extdata", "pgv", "public", "settings.json", package = "PGVdb")

  list(
    datafiles = datafiles.json,
    datadir = datadir,
    settings = settings
  )
}

reset_pgvdb  <- function() {
  devtools::load_all(".")
  paths <- load_paths()
  default_datafiles_json_path <- system.file("extdata", "pgv", "public", "datafiles0.json", package = "PGVdb")
  file.copy(default_datafiles_json_path, paths$datafiles, overwrite=TRUE)
  endpoint <- "http://10.1.29.225:8000/"
  pgvdb <- PGVdb$new(paths$datafiles, paths$datadir, paths$settings, higlass_metadata=list(endpoint=endpoint))
  return(pgvdb)
}


test_that("PGVdb initializes correctly", {
  pgvdb <- reset_pgvdb()
  expect_equal(nrow(pgvdb$metadata), 1)
  expect_equal(nrow(pgvdb$plots), 10)
})

test_that("load_json works correctly", {
  pgvdb <- reset_pgvdb()
  expect_error(pgvdb$load_json("bad_path.json"))

  paths <- load_paths()
  pgvdb$load_json(paths$datafiles)

  expect_equal(nrow(pgvdb$metadata), 1)
  expect_equal(nrow(pgvdb$plots), 10)
})

test_that("to_datatable returns correct output", {
  pgvdb <- reset_pgvdb()
  dt <- pgvdb$to_datatable()

  expect_s3_class(dt, "data.table")
  expect_gte(nrow(dt), 10)

  dt_filtered <- pgvdb$to_datatable(list("patient.id", "DEMO"))

  expect_equal(nrow(dt_filtered), 10)
})

test_that("add_plots loads from filepath correctly", {
  pgvdb <- reset_pgvdb()

  new_cov <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg19",
    tags = c("tags1", "tags2", "tags3"),
    x = system.file("extdata", "test_data", "test.cov.rds", package = "PGVdb"),
    field = "cn",
    visible = TRUE,
    type = "scatterplot",
    overwrite = TRUE
  )
  pgvdb$add_plots(new_cov)
  expect_equal(nrow(pgvdb$plots), 11)

  new_genome <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg38",
    x = system.file("extdata", "test_data", "test.gg.rds", package = "PGVdb"),
    visible = TRUE
  )
  pgvdb$add_plots(new_genome)
  expect_equal(nrow(pgvdb$plots), 12)

  new_walk <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg19",
    x = system.file("extdata", "test_data", "test.gw.rds", package = "PGVdb"),
    visible = TRUE
  )
  pgvdb$add_plots(new_walk)
  expect_equal(nrow(pgvdb$plots), 13)

  new_bw <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg19",
    x = list(list(server = "http://higlass.io", uuid = "AOR9BgKaS4WPX7esuBM4sQ")),
    visible = TRUE
  )
  pgvdb$add_plots(new_bw)
  expect_equal(nrow(pgvdb$plots), 14)

  new_json <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg19",
    type = "walk",
    x = system.file("extdata", "test_data", "walks.json", package = "PGVdb"),
    visible = TRUE
  )
  pgvdb$add_plots(new_json)
  expect_equal(nrow(pgvdb$plots), 15)

  new_bigwig_granges  <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg38",
    type = "bigwig",
    field = "foreground",
    x = system.file("extdata", "test_data", "test_bigwig_granges.rds", package = "PGVdb"),
    visible = TRUE
  )
  pgvdb$add_plots(new_bigwig_granges)
  expect_equal(nrow(pgvdb$plots), 16)
})

test_that("add_plots loads from object correctly", {
  pgvdb <- reset_pgvdb()

  cov = readRDS(system.file("extdata", "test_data", "test.cov.rds", package = "PGVdb"))
  new_cov <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg19",
    x = list(cov),
    field = "cn",
    type = "scatterplot",
    visible = TRUE,
    overwrite = TRUE
  )
  pgvdb$add_plots(new_cov)
  expect_equal(nrow(pgvdb$plots), 11)


  gg = readRDS(system.file("extdata", "test_data", "test.gg.rds", package = "PGVdb"))
  new_genome <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg19",
    x = list(gg),
    visible = TRUE
  )
  pgvdb$add_plots(new_genome)
  expect_equal(nrow(pgvdb$plots), 12)

  gw = readRDS(system.file("extdata", "test_data", "test.gw.rds", package = "PGVdb"))
  new_walk <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg19",
    x = list(gw),
    visible = TRUE
  )
  pgvdb$add_plots(new_walk)
  expect_equal(nrow(pgvdb$plots), 13)

  gr_bw = readRDS(system.file("extdata", "test_data", "test_bigwig_granges.rds", package = "PGVdb"))
  new_bigwig_granges  <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg38",
    type = "bigwig",
    field = "foreground",
    x = list(gr_bw),
    visible = TRUE,
    overwrite = TRUE
  )
  pgvdb$add_plots(new_bigwig_granges)
  expect_equal(nrow(pgvdb$plots), 14)
})

test_that("Bug with rds file", {
  pgvdb <- reset_pgvdb()
  new_plots <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg19",
    x = system.file("extdata", "test_data", "complex_not_added.rds", package = "PGVdb"),
    visible = TRUE
  )
  pgvdb$add_plots(new_plots)
})

test_that("add_plots works correctly with multiple plot filepaths", {
  pgvdb <- reset_pgvdb()

  paths  <- c(
    system.file("extdata", "test_data", "test.cov.rds", package = "PGVdb"),
    system.file("extdata", "test_data", "test.gg.rds", package = "PGVdb"),
    system.file("extdata", "test_data", "test.gw.rds", package = "PGVdb"),
    system.file("extdata", "test_data", "test_bigwig_granges.rds", package = "PGVdb")
  )

  new_plots <- data.table(
    patient.id = "TEST_ADD",
    ref = c("hg19", "hg19", "hg19", "hg38"),
    x = paths,
    field= c("cn", NA, NA, "foreground"),
    type=c("scatterplot", NA, NA, "bigwig"),
    visible = TRUE,
    overwrite = c(TRUE, TRUE, TRUE, TRUE)
  )
  pgvdb$add_plots(new_plots)
  expect_equal(nrow(pgvdb$plots), 14)
})

test_that("add_plots works correctly with multiple plot objects", {
  pgvdb <- reset_pgvdb()
  objects  <- c(
    list(readRDS(system.file("extdata", "test_data", "test.cov.rds", package = "PGVdb"))),
    list(readRDS(system.file("extdata", "test_data", "test.gg.rds", package = "PGVdb"))),
    list(readRDS(system.file("extdata", "test_data", "test.gw.rds", package = "PGVdb"))),
    list(readRDS(system.file("extdata", "test_data", "test_bigwig_granges.rds", package = "PGVdb")))
  )

  new_plots <- data.table(
    patient.id = "TEST_ADD",
    ref = c("hg19", "hg19", "hg19", "hg38"),
    x = objects,
    field = c("cn", NA, NA, "foreground"),
    type = c("scatterplot", NA, NA, "bigwig"),
    visible = TRUE,
    overwrite = c(TRUE, TRUE, TRUE, TRUE)
  )
  pgvdb$add_plots(new_plots)
  expect_equal(nrow(pgvdb$plots), 14)
})

test_that("add_plots works correctly with multiple bigwigs", {
  pgvdb <- reset_pgvdb()
  bigwigs  <- c(
    list(readRDS(system.file("extdata", "test_data", "test_bigwig_granges.rds", package = "PGVdb"))),
    list(readRDS(system.file("extdata", "test_data", "test_bigwig_granges.rds", package = "PGVdb")))
  )

  new_plots <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg38",
    x = bigwigs,
    field = "foreground",
    type = "bigwig",
    visible = TRUE
  )
  pgvdb$add_plots(new_plots)
  expect_equal(nrow(pgvdb$plots), 14)
})

test_that("add_plots works correctly with multiple patients", {
  pgvdb <- reset_pgvdb()
  paths  <- c(
    system.file("extdata", "test_data", "test.cov.rds", package = "PGVdb"),
    system.file("extdata", "test_data", "test.gg.rds", package = "PGVdb"),
    system.file("extdata", "test_data", "test.gw.rds", package = "PGVdb")
  )
  new_plots <- data.table(
    patient.id = c("TEST_ADD1", "TEST_ADD2", "TEST_ADD3"),
    ref = "hg19",
    x = paths,
    field= c("cn", NA, NA),
    type=c("scatterplot", NA, NA),
    visible = TRUE,
    overwrite = TRUE
  )
  pgvdb$add_plots(new_plots)
  expect_equal(nrow(pgvdb$plots), 13)
})

test_that("remove_plots works correctly", {
  pgvdb <- reset_pgvdb()
  paths  <- c(
    system.file("extdata", "test_data", "test.cov.rds", package = "PGVdb"),
    system.file("extdata", "test_data", "test.gg.rds", package = "PGVdb"),
    system.file("extdata", "test_data", "test.gw.rds", package = "PGVdb")
  )
  new_plots <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg19",
    x = paths,
    field = c("cn", NA, NA),
    type=c("scatterplot", NA, NA),
    visible = TRUE,
    overwrite = TRUE
  )
  pgvdb$add_plots(new_plots)

  remove_plot <- data.table(
    patient.id = "TEST_ADD",
    source = "coverage.arrow"
  )

  pgvdb$remove_plots(remove_plot)

  expect_equal(nrow(pgvdb$plots), 12)
})

test_that("remove_plots works correctly when removing patients", {
  pgvdb <- reset_pgvdb()
  paths  <- c(
    system.file("extdata", "test_data", "test.cov.rds", package = "PGVdb"),
    system.file("extdata", "test_data", "test.gg.rds", package = "PGVdb"),
    system.file("extdata", "test_data", "test.gw.rds", package = "PGVdb")
  )
  new_plots <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg19",
    x = paths,
    field = c("cn", NA, NA),
    type=c("scatterplot", NA, NA),
    visible = TRUE,
    overwrite = TRUE
  )
  pgvdb$add_plots(new_plots)

  remove_plot <- data.table(
    patient.id = "TEST_ADD"
  )

  pgvdb$remove_plots(remove_plot)

  expect_equal(nrow(pgvdb$plots), 10)
})

test_that("validate works correctly", {
  # duplicate plots
  pgvdb  <- reset_pgvdb()
  non_dup_pgvdb  <- pgvdb$plots
  pgvdb$plots[, patient.id2 := list(patient.id)]
  setnames(pgvdb$plots, "patient.id2", "patient.id")
  expect_warning(pgvdb$validate())
  expect_equal(non_dup_pgvdb, pgvdb$plots)
})

test_that("listing higlass tilesets works correctly", {
  pgvdb  <- reset_pgvdb()
  pgvdb$higlass_metadata$endpoint <- "http://10.1.29.225:8000/"
  tilesets <- pgvdb$list_higlass_tilesets()
  print(tilesets)
})


test_that("adding to higlass server works correctly", {
  pgvdb <- reset_pgvdb()
  pgvdb$higlass_metadata$endpoint <- "http://10.1.29.225:8000/"
  pgvdb$upload_to_higlass(
    datafile = system.file("extdata", "test_data", "chromSizes.tsv", package = "PGVdb"),
    filetype = "chromsizes-tsv",
    datatype = "chromsizes",
    coordSystem = "hg38",
    name = "hg38"
  )
  pgvdb$upload_to_higlass(
    datafile = system.file("extdata", "test_data", "higlass_test_bigwig.bw", package = "PGVdb"),
    name = "test_bigwig",
    filetype = "bigwig",
    datatype = "vector",
    coordSystem = "hg38",
  )
  expect_equal(nrow(pgvdb$plots), 11)
})

test_that("deleting higlass tileset works correctly", {
  pgvdb <- reset_pgvdb()
  pgvdb$higlass_metadata$endpoint <- "http://10.1.29.225:8000/"

  # flush higlass
  tilesets <- pgvdb$list_higlass_tilesets()
  uuids <- tilesets$uuid
  pgvdb$delete_from_higlass(pgvdb$higlass_metadata$endpoint, uuids = uuids)

  uuid  <- pgvdb$plots[11, "uuid"]
  pgvdb$delete_from_higlass(pgvdb$higlass_metadata$endpoint, uuid = uuid[[1]])
  expect_equal(nrow(pgvdb$plots), 10)
})
    
test_that("init_pgv works correctly", {
  pgvdb <- reset_pgvdb()
  paths  <- c(
    system.file("extdata", "test_data", "test.cov.rds", package = "PGVdb"),
    system.file("extdata", "test_data", "test.gg.rds", package = "PGVdb"),
    system.file("extdata", "test_data", "test.gw.rds", package = "PGVdb")
  )
  new_plots <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg19",
    type = c("scatterplot", "genome", "walk"),
    path = paths,
    source = c("coverage.arrow", "genome.json", "walk.json"),
    field = c("cn", NA, NA),
    visible = TRUE
  )
  pgvdb$add_plots(new_plots)

  pgv_dir  <- "/Users/diders01/projects/pgv_init_test"
  pgvdb$init_pgv(pgv_dir)
})
