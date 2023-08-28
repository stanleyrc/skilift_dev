library(testthat)

setup({
  library(parallel)
  library(R6)
  library(data.table)
  library(jsonlite)
  library(httr)
  devtools::load_all("../../../gGnome/gGnome")
})

devtools::load_all(".")
context("PGVdb")

load_paths <- function() {
  datafiles.json <- system.file("extdata", "pgv", "public", "datafiles.json", package = "PGVdb")
  datadir <- system.file("extdata", "pgv", "public", "data", package = "PGVdb")
  publicdir <- system.file("extdata", "pgv", "public", package = "PGVdb")
  settings <- system.file("extdata", "pgv", "public", "settings.json", package = "PGVdb")

  list(
    datafiles = datafiles.json,
    datadir = datadir,
    publicdir = publicdir,
    settings = settings
  )
}

reset_pgvdb  <- function() {
  devtools::load_all(".")
  paths <- load_paths()
  default_datafiles_json_path <- system.file("extdata", "pgv", "public", "datafiles0.json", package = "PGVdb")
  file.copy(default_datafiles_json_path, paths$datafiles, overwrite = TRUE)
  pgvdb <- PGVdb$new(paths$datafiles, paths$datadir, paths$publicdir, paths$settings)
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
    x = system.file("extdata", "test_data", "test.cov.rds", package = "PGVdb"),
    field = "cn",
    visible = TRUE
  )
  pgvdb$add_plots(new_cov, overwrite = TRUE)

  new_genome <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg19",
    x = system.file("extdata", "test_data", "test.gg.rds", package = "PGVdb"),
    visible = TRUE
  )
  pgvdb$add_plots(new_genome)

  new_walk <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg19",
    x = system.file("extdata", "test_data", "test.gw.rds", package = "PGVdb"),
    visible = TRUE
  )
  pgvdb$add_plots(new_walk)

  new_bw <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg19",
    x = list(list(server = "http://higlass.io", uuid = "AOR9BgKaS4WPX7esuBM4sQ")),
    visible = TRUE
  )
  pgvdb$add_plots(new_bw)
  expect_warning(pgvdb$add_plots(new_bw)) # Try adding duplicate

  new_json <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg19",
    type = "walk",
    x = system.file("extdata", "test_data", "walks.json", package = "PGVdb"),
    visible = TRUE
  )
  pgvdb  <- reset_pgvdb()
  pgvdb$add_plots(new_json)
})

test_that("add_plots loads from object correctly", {
  pgvdb <- reset_pgvdb()

  cov = readRDS(system.file("extdata", "test_data", "test.cov.rds", package = "PGVdb"))
  new_cov <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg19",
    x = list(cov),
    field = "cn",
    visible = TRUE
  )
  pgvdb$add_plots(new_cov)

  gg = readRDS(system.file("extdata", "test_data", "test.gg.rds", package = "PGVdb"))
  new_genome <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg19",
    x = list(gg),
    visible = TRUE
  )
  pgvdb$add_plots(new_genome)

  gw = readRDS(system.file("extdata", "test_data", "test.gw.rds", package = "PGVdb"))
  new_walk <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg19",
    x = list(gw),
    visible = TRUE
  )
  pgvdb  <- reset_pgvdb()
  pgvdb$add_plots(new_walk)

})

test_that("add_plots works correctly with multiple plot filepaths", {
  pgvdb <- reset_pgvdb()
  new_plots <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg19",
    x = system.file("extdata", "test_data", "complex_not_added.rds", package = "PGVdb"),
    visible = TRUE
  )
  pgvdb$add_plots(new_plots)

  paths  <- c(
    system.file("extdata", "test_data", "test.cov.rds", package = "PGVdb"),
    system.file("extdata", "test_data", "test.gg.rds", package = "PGVdb"),
    system.file("extdata", "test_data", "test.gw.rds", package = "PGVdb")
  )
  new_plots <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg19",
    x = paths,
    field= c("cn", NA, NA),
    visible = TRUE
  )
  pgvdb$add_plots(new_plots, overwrite=TRUE)
})

test_that("add_plots works correctly with multiple plot objects", {
  pgvdb <- reset_pgvdb()
  objects  <- c(
    list(readRDS(system.file("extdata", "test_data", "test.cov.rds", package = "PGVdb"))),
    list(readRDS(system.file("extdata", "test_data", "test.gg.rds", package = "PGVdb"))),
    list(readRDS(system.file("extdata", "test_data", "test.gw.rds", package = "PGVdb")))
  )
  new_plots <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg19",
    x = objects,
    field= c("cn", NA, NA),
    visible = TRUE
  )
  pgvdb$add_plots(new_plots, overwrite=TRUE)
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
    visible = TRUE
  )
  pgvdb$add_plots(new_plots, overwrite=TRUE)
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
    visible = TRUE
  )
  pgvdb$add_plots(new_plots, overwrite=TRUE)

  remove_plot <- data.table(
    patient.id = "TEST_ADD",
    source = "coverage.arrow"
  )

  pgvdb$remove_plots(remove_plot)

  expect_equal(nrow(pgvdb$plots), 10)
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
    visible = TRUE
  )
  pgvdb$add_plots(new_plots, overwrite=TRUE)

  remove_plot <- data.table(
    patient.id = "TEST_ADD"
  )

  pgvdb$remove_plots(remove_plot)

  expect_equal(nrow(pgvdb$plots), 10)
})

test_that("validate works correctly", {
  expect_silent(pgvdb$validate())

})

test_that("adding to higlass server works correctly", {
  pgvdb <- reset_pgvdb()
  # upload_to_higlass = function(datafile, filetype, datatype, coordSystem, name) {
  endpoint <- "http://10.1.29.225:8000/api/v1/tilesets/"
  pgvdb$upload_to_higlass(
    endpoint,
    datafile = system.file("extdata", "test_data", "chromSizes.tsv", package = "PGVdb"),
    filetype = "chromsizes-tsv",
    datatype = "chromsizes",
    coordSystem = "hg38",
    name = "hg38"
  )
  pgvdb$upload_to_higlass(
    endpoint,
    datafile = system.file("extdata", "test_data", "higlass_test_bigwig.bw", package = "PGVdb"),
    name = "test_bigwig",
    filetype = "bigwig",
    datatype = "vector",
    coordSystem = "hg38",
    uuid="test"
  )
})

test_that("deleting higlass tileset works correctly", {
  pgvdb <- reset_pgvdb()
  # upload_to_higlass = function(datafile, filetype, datatype, coordSystem, name) {
  endpoint <- "http://10.1.29.225:8000/api/v1/tilesets/"
  pgvdb$delete_from_higlass(endpoint, uuid = "VgApTJAsRHipjwtHJNX1IA")
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
