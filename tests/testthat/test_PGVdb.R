library(testthat)

setup({
  library(parallel)
  library(R6)
  library(data.table)
  library(jsonlite)
  library(httr)
                                        #  devtools::load_all("../../../gGnome")
  devtools::load_all("/gpfs/commons/groups/imielinski_lab/home/sclarke/git/gGnome_dev")
  setDTthreads(1)
})

##devtools::load_all(".")
devtools::load_all("/gpfs/commons/groups/imielinski_lab/home/sclarke/git/pgvdb_dev")
context("PGVdb")

load_paths <- function() {
  ## datafiles.json <- system.file("extdata", "pgv", "public", "datafiles.json", package = "PGVdb")
  ## datadir <- system.file("extdata", "pgv", "public", "data", package = "PGVdb")
  ## settings <- system.file("extdata", "pgv", "public", "settings.json", package = "PGVdb")

    datafiles.json <- "/gpfs/commons/groups/imielinski_lab/home/sclarke/git/pgv_testing/public/datafiles.json"
    default_datafiles_json_path <- "/gpfs/commons/groups/imielinski_lab/home/sclarke/git/pgv_testing/public/datafiles0.json"
    datadir <- "/gpfs/commons/groups/imielinski_lab/home/sclarke/git/pgv_testing/public/data"
    settings <- "/gpfs/commons/groups/imielinski_lab/home/sclarke/git/pgv_testing/public/settings.json"
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
  expect_equal(nrow(pgvdb$plots), 13)
})

## I am not sure how this was used before so I am just reading them in for all tests
##default_datafiles_json_path <- "/gpfs/commons/groups/imielinski_lab/home/sclarke/git/pgv_testing/public/datafiles.json"
## datafiles.json <- "/gpfs/commons/groups/imielinski_lab/home/sclarke/git/pgv_testing/public/datafiles.json"
## empty_datafiles_json_path <- "/gpfs/commons/groups/imielinski_lab/home/sclarke/git/pgv_testing/public/datafiles_empty.json"
## datadir <- "/gpfs/commons/groups/imielinski_lab/home/sclarke/git/pgv_testing/public/data"
## settings <- "/gpfs/commons/groups/imielinski_lab/home/sclarke/git/pgv_testing/public/settings.json"


test_that("PGVdb initializes from empty datafiles.json", {
  devtools::load_all(".")
  ## datafiles.json <- "/Users/diders01/projects/pgvdb/inst/extdata/pgv/public/datafiles_empty.json"
  datadir <- system.file("extdata", "pgv", "public", "data", package = "PGVdb")
  settings <- system.file("extdata", "pgv", "public", "settings.json", package = "PGVdb")
  ## datafiles.json <- "/gpfs/commons/groups/imielinski_lab/home/sclarke/git/pgv_for_testingpgvdb/public/datafiles.json"
  ## datadir <- "/gpfs/commons/groups/imielinski_lab/home/sclarke/git/pgv_for_testingpgvdb/public/data"
  ## settings <- "/gpfs/commons/groups/imielinski_lab/home/sclarke/git/pgv_for_testingpgvdb/public/settings.json"
  ## default_datafiles_json_path <- system.file("extdata", "pgv", "public", "datafiles0.json", package = "PGVdb")
  empty_datafiles_json_path <- "/gpfs/commons/groups/imielinski_lab/home/sclarke/git/pgv_testing/public/datafiles_empty.json"
  endpoint <- "http://10.1.29.225:8000/"
  ## pgvdb <- PGVdb$new(datafiles.json, datadir, settings, higlass_metadata=list(endpoint=endpoint))
  pgvdb <- PGVdb$new(empty_datafiles_json_path, datadir, settings, higlass_metadata=list(endpoint=endpoint))
  expect_equal(nrow(pgvdb$metadata), 0)
  expect_equal(nrow(pgvdb$plots), 0)
})


######THIS FAILS - not sure why the bad_path.json works
test_that("load_json works correctly", {
  pgvdb <- reset_pgvdb()
  expect_error(pgvdb$load_json("bad_path.json"))

  paths <- load_paths()
  pgvdb$load_json(paths$datafiles)
  pgvdb$plots
  ##expect_equal(nrow(pgvdb$metadata), 1)
  expect_equal(nrow(pgvdb$plots), 13)
})


##passes
test_that("to_datatable returns correct output", {
  pgvdb <- reset_pgvdb()
  dt <- pgvdb$to_datatable()

  expect_s3_class(dt, "data.table")
  ## expect_gte(nrow(dt), 10)
  nrow1 = nrow(dt) ## demo data has been updated so this test whether it is equal is not relevant

  dt_filtered <- pgvdb$to_datatable(list("patient.id", "DEMO"))

  expect_equal(nrow(dt_filtered), nrow1)
})


#######only warnings - passes
test_that("You can add plots to empty datafiles.json", {
  devtools::load_all(".")
  pgvdb <- reset_pgvdb()
  ## datafiles.json <- "/Users/diders01/projects/pgvdb/inst/extdata/pgv/public/datafiles_empty.json"
  empty_datafiles_json_path <- "/gpfs/commons/groups/imielinski_lab/home/sclarke/git/pgv_testing/public/datafiles_empty.json"
  datafiles.json <- empty_datafiles_json_path
  datadir <- system.file("extdata", "pgv", "public", "data", package = "PGVdb")
  settings <- system.file("extdata", "pgv", "public", "settings.json", package = "PGVdb")
  default_datafiles_json_path <- system.file("extdata", "pgv", "public", "datafiles0.json", package = "PGVdb")
  endpoint <- "http://10.1.29.225:8000/"
  ##  pgvdb <- PGVdb$new(datafiles.json, datadir, settings, higlass_metadata=list(endpoint=endpoint))
  nrow1 = nrow(pgvdb$plots)
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
  ## expect_equal(nrow(pgvdb$metadata), 1)
  ## expect_equal(nrow(pgvdb$plots), 1)
  expect_equal(nrow(pgvdb$metadata), 2)
  expect_equal(nrow(pgvdb$plots), 14)
})

test_that("add_plots loads from filepath correctly", {
    pgvdb <- reset_pgvdb()
    nrow1 = nrow(pgvdb$plots)
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
    ## expect_equal(nrow(pgvdb$plots), 11)
    expect_equal(nrow(pgvdb$plots), nrow1) ##14 now

    new_genome <- data.table(
        patient.id = "TEST_ADD",
        ref = "hg38",
        x = system.file("extdata", "test_data", "test.gg.rds", package = "PGVdb"),
        visible = TRUE
    )
    pgvdb$add_plots(new_genome)
    expect_equal(nrow(pgvdb$plots), nrow1+1) ## 15 now

    new_walk <- data.table(
        patient.id = "TEST_ADD",
        ref = "hg19",
        x = system.file("extdata", "test_data", "test.gw.rds", package = "PGVdb"),
        visible = TRUE
    )
    pgvdb$add_plots(new_walk)
    expect_equal(nrow(pgvdb$plots), nrow1+2) ##16 now

    new_bw <- data.table(
        patient.id = "TEST_ADD",
        ref = "hg19",
        x = list(list(server = "http://higlass.io", uuid = "AOR9BgKaS4WPX7esuBM4sQ")),
        visible = TRUE
    )
    pgvdb$add_plots(new_bw)
    expect_equal(nrow(pgvdb$plots), nrow1+3) ##17 now

    new_json <- data.table(
        patient.id = "TEST_ADD",
        ref = "hg19",
        type = "walk",
        x = system.file("extdata", "test_data", "walks.json", package = "PGVdb"),
        visible = TRUE
    )
    pgvdb$add_plots(new_json)
    expect_equal(nrow(pgvdb$plots), nrow1+4) ##18 now

    ##not sure why this added to plots because test_bigwig_granges.rds does not exist anymore
    new_bigwig_granges  <- data.table(
        patient.id = "TEST_ADD",
        ref = "hg38",
        type = "bigwig",
        field = "foreground",
        x = system.file("extdata", "test_data", "test_bigwig_granges.rds", package = "PGVdb"), ## this file does not exists
        visible = TRUE
    )
    pgvdb$add_plots(new_bigwig_granges)
    expect_equal(nrow(pgvdb$plots), nrow1+5)
})

test_that("add_plots loads from object correctly", {
    pgvdb <- reset_pgvdb()
    nrow1 = nrow(pgvdb$plots)
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
    expect_equal(nrow(pgvdb$plots), nrow1)


    gg = readRDS(system.file("extdata", "test_data", "test.gg.rds", package = "PGVdb"))
    new_genome <- data.table(
        patient.id = "TEST_ADD",
        ref = "hg19",
        x = list(gg),
        visible = TRUE
    )
    pgvdb$add_plots(new_genome)
    expect_equal(nrow(pgvdb$plots), nrow1 + 1)

    gw = readRDS(system.file("extdata", "test_data", "test.gw.rds", package = "PGVdb"))
    new_walk <- data.table(
        patient.id = "TEST_ADD",
        ref = "hg19",
        x = list(gw),
        visible = TRUE
    )
    pgvdb$add_plots(new_walk)
    expect_equal(nrow(pgvdb$plots), nrow1 + 2)
    ##this one fails because the file does not exist
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
    expect_equal(nrow(pgvdb$plots), nrow1 + 3)
})


## this is failing- not sure what the expectation is
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

## this is failing but going to be testing this with new tests
test_that("add_plots works correctly with multiple plot filepaths", {
    pgvdb <- reset_pgvdb()
    nrow1 = nrow(pgvdb$plots)
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
    expect_equal(nrow(pgvdb$plots), nrow1 + 4)
})


##fails not sure where this data is
test_that("add_plots works correctly with multiple plot objects", {
    pgvdb <- reset_pgvdb()
    nrow = nrow(pgvdb$plots)
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
    expect_equal(nrow(pgvdb$plots), nrow1 + 4)
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

test_that("mixing higlass upload with adding plots works correctly", {
  pgvdb  <- reset_pgvdb()
  pgvdb$higlass_metadata$endpoint <- "http://10.1.29.225:8000"
  # Mix
  paths  <- c(
    system.file("extdata", "test_data", "test_higlass_mix_granges1.rds", package = "PGVdb"),
    system.file("extdata", "test_data", "test_higlass_mix_granges2.rds", package = "PGVdb"),
    system.file("extdata", "test_data", "test_higlass_mix_complex1.rds", package = "PGVdb"),
    system.file("extdata", "test_data", "test_higlass_mix_complex2.rds", package = "PGVdb")
  )

  gg <- system.file("extdata", "test_data", "test_higlass_mix_complex2.rds", package = "PGVdb")

  new_plots <- data.table(
    patient.id = c("TEST_ADD1", "TEST_ADD2", "TEST_ADD1", "TEST_ADD2"),
    ref = "hg38_chr",
    x = paths,
    field= c("score", "score", NA, NA),
    type=c("bigwig", "bigwig", NA, NA),
    visible = TRUE,
    overwrite = TRUE
  )

  # gGraphs only
  paths  <- c(
    system.file("extdata", "test_data", "test_higlass_mix_complex1.rds", package = "PGVdb"),
    system.file("extdata", "test_data", "test_higlass_mix_complex2.rds", package = "PGVdb")
  )

  new_plots <- data.table(
    patient.id = c("TEST_ADD1", "TEST_ADD2"),
    ref = "hg38_chr",
    x = paths,
    visible = TRUE,
    overwrite = TRUE
  )

  # Bigwigs Only
  paths  <- c(
    system.file("extdata", "test_data", "test_higlass_mix_granges1.rds", package = "PGVdb"),
    system.file("extdata", "test_data", "test_higlass_mix_granges2.rds", package = "PGVdb")
  )

  new_plots <- data.table(
    patient.id = c("TEST_ADD1", "TEST_ADD2"),
    ref = "hg38_chr",
    field= c("score", "score"),
    type=c("bigwig", "bigwig"),
    x = paths,
    visible = TRUE,
    overwrite = TRUE
  )
  pgvdb$add_plots(new_plots)
}

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
  pgvdb$higlass_metadata$endpoint <- "http://10.1.29.225:8000"
  tilesets <- pgvdb$list_higlass_tilesets()
  print(tilesets)
})


test_that("adding to higlass server works correctly", {
  pgvdb <- reset_pgvdb()
  pgvdb$higlass_metadata$endpoint <- "http://10.1.29.225:8000"
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
  pgvdb$higlass_metadata$endpoint <- "http://10.1.29.225:8000"

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

### Debugging

test_that("adding arrows in parallel works correctly", {
  pgv <- reset_pgvdb()
  maska_path = system.file("extdata", "test_data", "maskA_re.rds", package = "PGVdb")
  maska = readRDS(maska_path)
  maska$mask = "mask"

  pairs_path = system.file("extdata", "test_data", "casereport_pairs.rds", package = "PGVdb")
  pairs = readRDS(pairs_path)
  covs.lst = mclapply(pairs$pair, function(pair) {
    cov_gr = readRDS(pairs[pair,decomposed_cov])
    jab = readRDS(pairs[pair,complex])
    cov_gr$foregroundabs = rel2abs(gr=cov_gr,
      purity = jab$meta$purity,
      ploidy = jab$meta$ploidy,
      field="foreground"
    )
    cov_gr2 = rebin(cov_gr, 1e4, field = "foregroundabs")
    cov_gr3 = gr.val(cov_gr2, maska, "mask")
    cov_gr3 = cov_gr3 %Q% (mask != "mask")
    cov_gr3$mask = NULL
    plot_to_add = data.table(patient.id = pair,
      visible = TRUE,
      x = list(cov_gr3),
      type = "scatterplot",
      field = "foregroundabs",
      ref="hg38_chr",
      title = "Coverage rel2abs"
    )
    return(plot_to_add)
  }, mc.cores = 40)

  covs.dt = rbindlist(covs.lst)
  covs.dt[,ref := "hg19"]
  covs.dt[,title := "Masked Coverage rel2abs"]

  pgv <- reset_pgvdb()
  pgv$add_plots(covs.dt[1:1000], cores = 40)
})






#########Stanley new tests
## reset to demo here
library(JaBbA) ## was not working to load later - said gGnome could not be found- never had this happen so loaded here
library(testthat)

setup({
    library(parallel)
    library(R6)
    library(data.table)
    library(jsonlite)
    library(httr)
                                        #  devtools::load_all("../../../gGnome")
    devtools::load_all("/gpfs/commons/groups/imielinski_lab/home/sclarke/git/gGnome_dev")
    setDTthreads(1)
})

##devtools::load_all(".")
devtools::load_all("/gpfs/commons/groups/imielinski_lab/home/sclarke/git/pgvdb_dev")
devtools::load_all("/gpfs/commons/groups/imielinski_lab/home/sclarke/git/gGnome_dev")
context("PGVdb")

getPGV = function() {
    public_dir = "~/git/pgv_testing/public/"
    json_file = paste0(public_dir,"datafiles.json")
    datadir = paste0(public_dir,"data/")
    settings = paste0(public_dir,"settings.json")
    higlass.list = list(endpoint = "https://higlass01.nygenome.org/")
    pgv = PGVdb$new(datafiles_json_path = json_file,datadir = datadir,settings = settings,higlass_metadata = higlass.list)
    return(pgv)
}

test.pairs = readRDS("/gpfs/commons/home/sclarke/git/pgvdb_test_data/test_pairs.rds")
setkey(test.pairs,pair)
hg19.seq = readRDS("/gpfs/commons/home/sclarke/git/pgvdb_test_data/hg19.seq") #for bigwigs

test.meta.pairs = readRDS("/gpfs/commons/home/sclarke/git/pgvdb_test_data/test.meta_pairs.rds")
setkey(test.meta.pairs,pair)


##################################################################################################################################################################################################################

test_that("genome graphs add correctly using template ", {
    pgvdb = getPGV()
    nrow1 = nrow(pgvdb$plots)
    genomes.add = genome_temp(patient_id = test.pairs$pair, x = test.pairs$jabba_gg, annotation = NULL, ref = "hg19", order = NULL)
    pgvdb$add_plots(genomes.add, cores = 5)
    expect_equal(nrow(pgvdb$plots), nrow1 + 5)
})

test_that("genome graphs add correctly using template with type = NULL ", {
    pgvdb = getPGV()
    nrow1 = nrow(pgvdb$plots)
    ## add without type to make sure genome is added as type
    genomes.add = genome_temp(patient_id = test.pairs$pair, x = test.pairs$jabba_gg, annotation = NULL, ref = "hg19", order = NULL, type = NULL)
    pgvdb$add_plots(genomes.add, cores = 5)
    expect_equal(nrow(pgvdb$plots), nrow1 + 5)
})

test_that("genome graphs with annotations (events output) get added correctly", {
    pgvdb = getPGV()
    nrow1 = nrow(pgvdb$plots)
    genomes.add = genome_temp(patient_id = test.pairs$pair, x = test.pairs$complex, ref = "hg19", order = NULL)
    pgvdb$add_plots(genomes.add, cores = 5)
    expect_equal(nrow(pgvdb$plots), nrow1 + 5)
})

## this also tests arrow_temp because it is within the rebinning function (cov2arrow_pgv)
test_that("coverage plots as arrows get uploaded correctly", {
    pgvdb = getPGV()
    nrow1 = nrow(pgvdb$plots)
    add.lst = mclapply(test.pairs$pair, function(pair) {
        add.dt = cov2arrow_pgv(patient.id = pair, dryclean_cov = test.pairs[pair,tumor_dryclean_cov], ref = "hg19")
        return(add.dt)
    }, mc.cores = 5)
    rebin.cov.dt = rbindlist(add.lst)
    ## genomes.add = arrow_temp(patient_id = test.pairs$pair, x = test.pairs$tumor_dryclean_cov, ref = "hg19", order = NULL)
    pgvdb$add_plots(rebin.cov.dt, cores = 5)
    expect_equal(nrow(pgvdb$plots), nrow1 + 5)
})

## just doing one for bigwigs, some finicky seqlengths things- in the future implement a way to force the seqlengths to match- would be easier once we do not have _chr references
test_that("coverage plots as bigwigs upload correctly", {
    pgvdb = getPGV()
    nrow1 = nrow(pgvdb$plots)
    pair = test.pairs$pair[1]
    cov.gr = readRDS(test.pairs[pair,tumor_dryclean_cov])
    cov.gr2 = as.data.table(cov.gr) %>% GRanges(.,seqlengths = hg19.seq) %>% trim ## only works when trimmming - should implement into upload itself but may make it slower?
    bw.add = bw_temp(patient_id = pair, x = list(cov.gr2), ref = "hg19", order = NULL)
    pgvdb$add_plots(bw.add, cores = 1)
    expect_equal(nrow(pgvdb$plots), nrow1 + 1)
})


## need Flow and JaBbA for this at the moment- had to keep this as the original directory and not use the test location to get coverage from job output
test_that("ppfit plots upload correctly", {
    library(Flow); library(skitools) ## skitools needed for one function (dirr)
    pgvdb = getPGV()
    nrow1 = nrow(pgvdb$plots)
    ppfit.add = ppfit_temp(patient_id = test.pairs$pair, x = test.pairs$balanced_gg, ref = "hg19")
    pgvdb$add_plots(ppfit.add, cores = 5)
    expect_equal(nrow(pgvdb$plots), nrow1 + 5)
})


## need Flow and JaBbA for this at the moment- had to keep this as the original directory and not use the test location to get coverage from job output
test_that("allelic plots upload correctly", {
    pgvdb = getPGV()
    nrow1 = nrow(pgvdb$plots)
    allelic.add = genome_temp(patient_id = test.pairs$pair, x = test.pairs$balanced_gg_rds, ref = "hg19", type = "allelic", annotation = NULL)
    allelic.add = allelic.add[c(1,3:5),] ## I double checked this but for some reason the second test sample balanced_gg_rds is not actually an allelic graph so it fails (it's good that it failed but not for tests)
    pgvdb$add_plots(allelic.add, cores = 4)
    expect_equal(nrow(pgvdb$plots), nrow1 + 4)
})


test_that("mutation plots upload correctly", {
    pgvdb = getPGV()
    nrow1 = nrow(pgvdb$plots)
    mutations.add = mutations_temp(patient_id = test.pairs$pair,field = "total_copies", x = test.pairs$somatic_snv_cn, ref = "hg19")
    pgvdb$add_plots(mutations.add, cores = 5)
    expect_equal(nrow(pgvdb$plots), nrow1 + 5)
})


## one sample tests for meta_data_json and filtered_events_json

test_that("filtered events jsons created", {
    pgvdb = getPGV() ## just used to get path for output
    pair = test.meta.pairs$pair
    out.file = gsub("settings.json","test_filtered_events.json",pgvdb$settings)
    filtered_events_json(pair = pair,
                     oncotable = test.meta.pairs[pair,oncotable],
                     jabba_gg = test.meta.pairs[pair,complex],
                     out_file = out.file,
                     cgc_file = "/gpfs/commons/groups/imielinski_lab/DB/COSMIC/v99_GRCh37/cancer_gene_census_fixed.csv",
                     temp_fix = TRUE)
    expect_true(file.exists(out.file))
})



test_that("metadata json created", {
    pgvdb = getPGV() ## just used to get path for output
    pair = test.meta.pairs$pair
    out.file = gsub("settings.json","test_metadata.json",pgvdb$settings)
    meta.dt = meta_data_json(pair = pair,
                             out_file = out.file,
                             coverage = test.meta.pairs[pair,decomposed_cov],
                             jabba_gg = test.meta.pairs[pair,complex],
                             svaba_somatic_vcf = test.meta.pairs[pair,svaba_somatic_vcf],
                             seqnames_loh = c(1:22),
                             karyograph = test.meta.pairs[pair,karyograph_rds],
                             vcf = test.meta.pairs[pair,strelka2_somatic_filtered_variants],
                             tumor_type = test.meta.pairs[pair,tumor_type_final],
                             disease = test.meta.pairs[pair,disease],
                             primary_site = test.meta.pairs[pair,primary_site_simple],
                             inferred_sex = test.meta.pairs[pair,inferred_sex],
                             seqnames_genome_width = c(1:22,"X","Y"),
                             write_json = TRUE,
                             overwrite = FALSE)
    expect_true(file.exists(out.file))
})

##create distributions
test_that("metadata json created", {
    pgvdb = getPGV() ## just used to get path for output
    pair = test.meta.pairs$pair
    input.folder = "/gpfs/commons/groups/imielinski_lab/home/sclarke/git/pgvdb_test_data/case_report_data/"
    out.folder = gsub("settings.json","test_distributions_output",pgvdb$settings)
    cmd = paste0("mkdir ",out.folder)
    system(command = cmd)
    create_distributions(input.folder, out.folder, filter_patients = NULL)
    ##make sure all files exist
    files = paste0(out.folder,"/",c("coverageVariance.json", "ploidy.json", "snvCount.json", "tmb.json", "lohFraction.json", "purity.json", "svCount.json"))
    expect_true(all(file.exists(files)))
})


##create strelka_qc() test
test_that("strelka.qc.json created", {
    pgvdb = getPGV() ## just used to get path for output
    pair = test.meta.pairs$pair
    out.file = gsub("settings.json","test_strelka.qc.json",pgvdb$settings)
    strelka.qc = strelka_qc(vcf = test.meta.pairs[pair,svaba_somatic_vcf],
			    seqnames_genome_width = c(1:22,"X","Y"),
                            outfile = out.file,
                            write_json = TRUE,
                            return_table = TRUE)
    expect_true(file.exists(out.file))
})


