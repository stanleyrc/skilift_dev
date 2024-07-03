library(testthat)

setup({
  library(parallel)
  library(R6)
  library(data.table)
  library(jsonlite)
  library(httr)
  devtools::load_all("~/git/gGnome")
  setDTthreads(1)
})

devtools::load_all(".")
context("Skilift")

load_paths <- function() {
  package_name <- "Skilift"
  relative_path <- "extdata/pgv/public/datafiles.json"

  publicdir <- system.file("extdata", "pgv", "public", package = "Skilift")
  datafiles.json <- file.path(publicdir, "datafiles.json")
  empty.datafiles.json <- file.path(publicdir, "empty_datafiles.json")
  datadir <- system.file("extdata", "pgv", "public", "data", package = "Skilift")
  settings <- system.file("extdata", "pgv", "public", "settings.json", package = "Skilift")

  list(
    publicdir = publicdir,
    datafiles = datafiles.json,
    empty_datafiles = empty.datafiles.json,
    datadir = datadir,
    settings = settings
  )
}

reset_skilift  <- function() {
    devtools::load_all(".")
    paths <- load_paths()
    default_datafiles_json_path <- system.file("extdata", "pgv", "public", "datafiles0.json", package = "Skilift")
    file.copy(default_datafiles_json_path, paths$datafiles, overwrite=TRUE)
    endpoint <- "http://10.1.29.225:8000/"
    skilift <- Skilift$new(public_dir = paths$publicdir, higlass_metadata=list(endpoint=endpoint))
    return(skilift)
}

test_that("Skilift initializes correctly", {
  skilift <- reset_skilift()
  expect_equal(nrow(skilift$metadata), 1)
  expect_equal(nrow(skilift$plots), 13)
})

test_that("skilift initializes from public dir", {
  devtools::load_all(".")
  endpoint <- "http://10.1.29.225:8000/"
  paths <- load_paths()
  skilift <- Skilift$new(paths$publicdir, higlass_metadata=list(endpoint=endpoint))
  expect_equal(nrow(skilift$metadata), 1)
  expect_equal(nrow(skilift$plots), 13)
})

test_that("skilift initializes from individual paths", {
  devtools::load_all(".")
  endpoint <- "http://10.1.29.225:8000/"
  paths <- load_paths()
  skilift <- Skilift$new(
    datafiles_json_path = paths$datafiles,
    datadir = paths$datadir,
    settings = paths$settings,
    higlass_metadata=list(endpoint=endpoint)
  )
  expect_equal(nrow(skilift$metadata), 1)
  expect_equal(nrow(skilift$plots), 13)
})

test_that("skilift initializes from empty datafiles.json", {
  devtools::load_all(".")
  endpoint <- "http://10.1.29.225:8000/"
  paths <- load_paths()
  suppressWarnings({
    skilift <- Skilift$new(
      datafiles_json_path = paths$empty_datafiles,
      datadir = paths$datadir,
      settings = paths$settings,
      higlass_metadata=list(endpoint=endpoint)
    )
  })
  expect_equal(nrow(skilift$metadata), 0)
  expect_equal(nrow(skilift$plots), 0)
})


test_that("load_json works correctly", {
  skilift <- reset_skilift()
  expect_warning(skilift$load_json("bad_path.json"))

  paths <- load_paths()
  skilift$load_json(paths$datafiles)
  skilift$plots
  expect_equal(nrow(skilift$metadata), 1)
  expect_equal(nrow(skilift$plots), 13)
})


test_that("to_datatable returns correct output", {
  skilift <- reset_skilift()
  dt <- skilift$to_datatable()

  expect_s3_class(dt, "data.table")
  nrow1 = nrow(dt) ## demo data has been updated so this test whether it is equal is not relevant

  dt_filtered <- skilift$to_datatable(list("patient.id", "DEMO"))

  expect_equal(nrow(dt_filtered), nrow1)
})


test_that("You can add plots to empty datafiles.json", {
  skilift <- reset_skilift()
  endpoint <- "http://10.1.29.225:8000/"
  paths <- load_paths()
  skilift <- Skilift$new(public_dir=paths$publicdir, higlass_metadata=list(endpoint=endpoint))
  nrow1 = nrow(skilift$plots)
  new_cov <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg19",
    tags = c("tags1", "tags2", "tags3"),
    x = system.file("extdata", "test_data", "test.cov.rds", package = "Skilift"),
    field = "cn",
    visible = TRUE,
    type = "scatterplot",
    overwrite = TRUE
  )
  suppressWarnings({skilift$add_plots(new_cov)})
  expect_equal(nrow(skilift$metadata), 2)
  expect_equal(nrow(skilift$plots), 14)
})

test_that("add_plots loads from filepath correctly", {
    skilift <- reset_skilift()
    nrow1 = nrow(skilift$plots)
    new_cov <- data.table(
        patient.id = "TEST_ADD",
        ref = "hg19",
        tags = c("tags1", "tags2", "tags3"),
        x = system.file("extdata", "test_data", "test.cov.rds", package = "Skilift"),
        field = "cn",
        visible = TRUE,
        type = "scatterplot",
        overwrite = TRUE
    )
    skilift$add_plots(new_cov)
    expect_equal(nrow(skilift$plots), nrow1+1) ##14 now

    new_genome <- data.table(
        patient.id = "TEST_ADD",
        ref = "hg38",
        x = system.file("extdata", "test_data", "test.gg.rds", package = "Skilift"),
        visible = TRUE
    )
    skilift$add_plots(new_genome)
    expect_equal(nrow(skilift$plots), nrow1+2) ## 15 now

    new_walk <- data.table(
        patient.id = "TEST_ADD",
        ref = "hg19",
        x = system.file("extdata", "test_data", "test.gw.rds", package = "Skilift"),
        visible = TRUE
    )
    skilift$add_plots(new_walk)
    expect_equal(nrow(skilift$plots), nrow1+3) ##16 now

    new_bw <- data.table(
        patient.id = "TEST_ADD",
        ref = "hg19",
        x = list(list(server = "http://higlass.io", uuid = "AOR9BgKaS4WPX7esuBM4sQ")),
        visible = TRUE
    )
    skilift$add_plots(new_bw)
    expect_equal(nrow(skilift$plots), nrow1+4) ##17 now

    new_json <- data.table(
        patient.id = "TEST_ADD",
        ref = "hg19",
        type = "walk",
        x = system.file("extdata", "test_data", "walks.json", package = "Skilift"),
        visible = TRUE
    )
    skilift$add_plots(new_json)
    expect_equal(nrow(skilift$plots), nrow1+5) ##18 now

    # this test will fail outside of the nygc hpc because the bigwig file is gitignored (due to it being large)
    skilift <- reset_skilift()
    new_bigwig_granges  <- data.table(
        patient.id = "TEST_ADD",
        ref = "hg38",
        type = "bigwig",
        field = "foreground",
        x = file.path("/gpfs/commons/groups/imielinski_lab/home/sdider/Projects/skilift/skilift/inst/extdata/test_data/test_bigwig_granges.rds"), 
        visible = TRUE
    )
    skilift$add_plots(new_bigwig_granges)
    expect_equal(nrow(skilift$plots), nrow1+5)
})

test_that("add_plots loads from object correctly", {
    skilift <- reset_skilift()
    nrow1 = nrow(skilift$plots)
    cov = readRDS(system.file("extdata", "test_data", "test.cov.rds", package = "Skilift"))
    new_cov <- data.table(
        patient.id = "TEST_ADD",
        ref = "hg19",
        x = list(cov),
        field = "cn",
        type = "scatterplot",
        visible = TRUE,
        overwrite = TRUE
    )
    skilift$add_plots(new_cov)
    expect_equal(nrow(skilift$plots), nrow1)


    gg = readRDS(system.file("extdata", "test_data", "test.gg.rds", package = "Skilift"))
    new_genome <- data.table(
        patient.id = "TEST_ADD",
        ref = "hg19",
        x = list(gg),
        visible = TRUE
    )
    skilift$add_plots(new_genome)
    expect_equal(nrow(skilift$plots), nrow1 + 1)

    gw = readRDS(system.file("extdata", "test_data", "test.gw.rds", package = "Skilift"))
    new_walk <- data.table(
        patient.id = "TEST_ADD",
        ref = "hg19",
        x = list(gw),
        visible = TRUE
    )
    skilift$add_plots(new_walk)
    expect_equal(nrow(skilift$plots), nrow1 + 2)

    # this test fails outside of the nygc hpc because the bigwig file is gitignored (due to it being large)

    bw_gr_rds_path = file.path("/gpfs/commons/groups/imielinski_lab/home/sdider/Projects/skilift/skilift/inst/extdata/test_data/test_bigwig_granges.rds")
    gr_bw = readRDS(bw_gr_rds_path)
    new_bigwig_granges  <- data.table(
        patient.id = "TEST_ADD",
        ref = "hg38",
        type = "bigwig",
        field = "foreground",
        x = list(gr_bw),
        visible = TRUE,
        overwrite = TRUE
    )
    skilift$add_plots(new_bigwig_granges)
    expect_equal(nrow(skilift$plots), nrow1 + 3)
})

test_that("add_plots loads from gtrack object correctly", {
    skilift <- reset_skilift()
    nrow1 = nrow(skilift$plots)
    cov = readRDS(system.file("extdata", "test_data", "cov_gtrack.rds", package = "Skilift"))
    new_cov <- data.table(
        patient.id = "TEST_ADD",
        ref = "hg19",
        x = list(cov),
        field = "seg.mean",
        type = "scatterplot",
        visible = TRUE,
        overwrite = TRUE
    )
    skilift$add_plots(new_cov)
    expect_equal(nrow(skilift$plots), nrow1 + 1)

    # gGenome gtrack not yet implemented
    # gg_gt = readRDS(system.file("extdata", "test_data", "ggraph_gtrack.rds", package = "Skilift"))
    # seqlengths = seqlengths(gg_gt@data[[1]])
    # nodes = gg_gt@data[[1]]
    # edges = gg_gt@edges[[1]]
    # new_genome <- data.table(
    #     patient.id = "TEST_ADD",
    #     ref = "hg19",
    #     x = list(gg),
    #     visible = TRUE
    # )
    # skilift$add_plots(new_genome)
    # expect_equal(nrow(skilift$plots), nrow1 + 1)

    skilift <- reset_skilift()
    gt_bw = readRDS(system.file("extdata", "test_data", "bigwig_gtrack.rds", package = "Skilift"))
    new_bigwig_granges  <- data.table(
        patient.id = "TEST_ADD",
        ref = "hg38",
        type = "bigwig",
        field = "seg.mean",
        x = list(gt_bw),
        visible = TRUE,
        overwrite = TRUE
    )
    skilift$add_plots(new_bigwig_granges)
    expect_equal(nrow(skilift$plots), nrow1 + 2)
})

test_that("add_plots works correctly with multiple plot filepaths", {
    skilift <- reset_skilift()
    nrow1 = nrow(skilift$plots)
    paths  <- c(
        system.file("extdata", "test_data", "test.cov.rds", package = "Skilift"),
        system.file("extdata", "test_data", "test.gg.rds", package = "Skilift"),
        system.file("extdata", "test_data", "test.gw.rds", package = "Skilift")
    )

    new_plots <- data.table(
        patient.id = "TEST_ADD",
        ref = c("hg19", "hg19", "hg19"),
        x = paths,
        field= c("cn", NA, NA),
        type=c("scatterplot", "genome", NA),
        visible = TRUE,
        overwrite = c(TRUE, TRUE, TRUE)
    )
    skilift$add_plots(new_plots)
    expect_equal(nrow(skilift$plots), nrow1 + 3)
})


test_that("add_plots works correctly with multiple plot objects", {
    skilift <- reset_skilift()
    nrow = nrow(skilift$plots)
    objects  <- c(
        list(readRDS(system.file("extdata", "test_data", "test.cov.rds", package = "Skilift"))),
        list(readRDS(system.file("extdata", "test_data", "test.gg.rds", package = "Skilift"))),
        list(readRDS(system.file("extdata", "test_data", "test.gw.rds", package = "Skilift")))
    )

    new_plots <- data.table(
        patient.id = "TEST_ADD",
        ref = c("hg19", "hg19", "hg19"),
        x = objects,
        field = c("cn", NA, NA),
        type = c("scatterplot", "genome", NA),
        visible = TRUE,
        overwrite = c(TRUE, TRUE, TRUE)
    )
    skilift$add_plots(new_plots)
    expect_equal(nrow(skilift$plots), nrow1 + 3)
})

# this test will fail outside of the nygc hpc because the bigwig file is gitignored (due to it being large)
test_that("add_plots works correctly with multiple bigwigs", {
  skilift <- reset_skilift()

  bw_gr_rds_path = file.path("/gpfs/commons/groups/imielinski_lab/home/sdider/Projects/skilift/skilift/inst/extdata/test_data/test_bigwig_granges.rds")
  bigwigs  <- c(
    list(readRDS(bw_gr_rds_path)),
    list(readRDS(bw_gr_rds_path))
  )

  new_plots <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg38",
    x = bigwigs,
    field = "foreground",
    type = "bigwig",
    visible = TRUE
  )
  skilift$add_plots(new_plots)
  expect_equal(nrow(skilift$plots), 14)
})

test_that("add_plots works correctly with multiple patients", {
  skilift <- reset_skilift()
  nrow = nrow(skilift$plots)
  paths  <- c(
    system.file("extdata", "test_data", "test.cov.rds", package = "Skilift"),
    system.file("extdata", "test_data", "test.gg.rds", package = "Skilift"),
    system.file("extdata", "test_data", "test.gw.rds", package = "Skilift")
  )
  new_plots <- data.table(
    patient.id = c("TEST_ADD1", "TEST_ADD2", "TEST_ADD3"),
    ref = "hg19",
    x = paths,
    field= c("cn", NA, NA),
    type=c("scatterplot", "genome", NA),
    visible = TRUE,
    overwrite = TRUE
  )
  skilift$add_plots(new_plots)
  expect_equal(nrow(skilift$plots), nrow+3)
})

test_that("mixing higlass upload with adding plots works correctly", {
  skilift  <- reset_skilift()
  skilift$higlass_metadata$endpoint <- "http://10.1.29.225:8000"
  # Mix
  paths  <- c(
    system.file("extdata", "test_data", "test_higlass_mix_granges1.rds", package = "Skilift"),
    system.file("extdata", "test_data", "test_higlass_mix_granges2.rds", package = "Skilift"),
    system.file("extdata", "test_data", "test_higlass_mix_complex1.rds", package = "Skilift"),
    system.file("extdata", "test_data", "test_higlass_mix_complex2.rds", package = "Skilift")
  )

  gg <- system.file("extdata", "test_data", "test_higlass_mix_complex2.rds", package = "Skilift")

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
    system.file("extdata", "test_data", "test_higlass_mix_complex1.rds", package = "Skilift"),
    system.file("extdata", "test_data", "test_higlass_mix_complex2.rds", package = "Skilift")
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
    system.file("extdata", "test_data", "test_higlass_mix_granges1.rds", package = "Skilift"),
    system.file("extdata", "test_data", "test_higlass_mix_granges2.rds", package = "Skilift")
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
  skilift$add_plots(new_plots)

})

test_that("remove_plots works correctly", {
  skilift <- reset_skilift()
  nrow = nrow(skilift$plots)
  paths  <- c(
    system.file("extdata", "test_data", "test.cov.rds", package = "Skilift"),
    system.file("extdata", "test_data", "test.gg.rds", package = "Skilift"),
    system.file("extdata", "test_data", "test.gw.rds", package = "Skilift")
  )
  new_plots <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg19",
    x = paths,
    field = c("cn", NA, NA),
    type=c("scatterplot", "genome", NA),
    visible = TRUE,
    overwrite = TRUE
  )
  skilift$add_plots(new_plots)

  remove_plot <- data.table(
    patient.id = "TEST_ADD",
    source = "coverage.arrow"
  )

  skilift$remove_plots(remove_plot)

  expect_equal(nrow(skilift$plots), nrow+3-1)
})

test_that("remove_plots works correctly when removing patients", {
  skilift <- reset_skilift()
  nrow = nrow(skilift$plots)
  paths  <- c(
    system.file("extdata", "test_data", "test.cov.rds", package = "Skilift"),
    system.file("extdata", "test_data", "test.gg.rds", package = "Skilift"),
    system.file("extdata", "test_data", "test.gw.rds", package = "Skilift")
  )
  new_plots <- data.table(
    patient.id = "TEST_ADD",
    ref = "hg19",
    x = paths,
    field = c("cn", NA, NA),
    type=c("scatterplot", "genome", NA),
    visible = TRUE,
    overwrite = TRUE
  )
  skilift$add_plots(new_plots)

  remove_plot <- data.table(
    patient.id = "TEST_ADD"
  )

  skilift$remove_plots(remove_plot)

  expect_equal(nrow(skilift$plots), nrow)
})

test_that("validate works correctly", {
  # duplicate plots
  skilift  <- reset_skilift()
  non_dup_skilift  <- skilift$plots
  skilift$plots[, patient.id2 := list(patient.id)]
  setnames(skilift$plots, "patient.id2", "patient.id")
  expect_warning(skilift$validate())
  expect_equal(non_dup_skilift, skilift$plots)
})

test_that("listing higlass tilesets works correctly", {
  skilift  <- reset_skilift()
  skilift$higlass_metadata$endpoint <- "http://10.1.29.225:8000"
  tilesets <- skilift$list_higlass_tilesets()
  print(tilesets)
})


test_that("adding to higlass server works correctly", {
  skilift <- reset_skilift()
  skilift$higlass_metadata$endpoint <- "http://10.1.29.225:8000"
  skilift$upload_to_higlass(
    datafile = system.file("extdata", "test_data", "chromSizes.tsv", package = "Skilift"),
    filetype = "chromsizes-tsv",
    datatype = "chromsizes",
    coordSystem = "hg38",
    name = "hg38"
  )
  skilift$upload_to_higlass(
    datafile = system.file("extdata", "test_data", "higlass_test_bigwig.bw", package = "Skilift"),
    name = "test_bigwig",
    filetype = "bigwig",
    datatype = "vector",
    coordSystem = "hg38",
  )
  expect_equal(nrow(skilift$plots), 11)
})

test_that("deleting higlass tileset works correctly", {
  skilift <- reset_skilift()
  skilift$higlass_metadata$endpoint <- "http://10.1.29.225:8000"

  # flush higlass
  tilesets <- skilift$list_higlass_tilesets()
  uuids <- tilesets$uuid
  skilift$delete_from_higlass(skilift$higlass_metadata$endpoint, uuids = uuids)

  uuid  <- skilift$plots[11, "uuid"]
  skilift$delete_from_higlass(skilift$higlass_metadata$endpoint, uuid = uuid[[1]])
  expect_equal(nrow(skilift$plots), 10)
})

    
test_that("init_pgv works correctly", {
  skilift <- reset_skilift()
  paths  <- c(
    system.file("extdata", "test_data", "test.cov.rds", package = "Skilift"),
    system.file("extdata", "test_data", "test.gg.rds", package = "Skilift"),
    system.file("extdata", "test_data", "test.gw.rds", package = "Skilift")
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
  skilift$add_plots(new_plots)

  pgv_dir  <- "/Users/diders01/projects/pgv_init_test"
  skilift$init_pgv(pgv_dir)
})

test_that("you can mix pgv and case-report datafiles/datadir", {
  skilift <- reset_skilift()
  skilift <- Skilift$new(public_dir = pgv_public_dir, datadir=casereport_datadir, higlass_metadata=list(endpoint=endpoint))

})

### json generation methods

# Sage
test_that("sage methods work correctly", {
  pgv <- reset_skilift()
  test_meta_pairs = readRDS("/gpfs/commons/home/sclarke/git/pgvdb_test_data/test.meta_pairs.rds")
  sage_vcf_gz = test_meta_pairs[1]$sage_somatic_variants
  sage_dt = sage_qc(sage_vcf_gz, "hg19", write_json = FALSE, return_table = TRUE)
  sage_count = sage_count(sage_vcf_gz, "hg19")

  expect_true(!is.null(sage_dt) && length(sage_dt) > 0)
  expect_true(!is.null(sage_count) && length(sage_count) > 0)
})


### Debugging

test_that("adding arrows in parallel works correctly", {
  pgv <- reset_skilift()
  maska_path = system.file("extdata", "test_data", "maskA_re.rds", package = "Skilift")
  maska = readRDS(maska_path)
  maska$mask = "mask"

  pairs_path = system.file("extdata", "test_data", "casereport_pairs.rds", package = "Skilift")
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

  pgv <- reset_skilift()
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
devtools::load_all("/gpfs/commons/groups/imielinski_lab/home/sclarke/git/skilift_dev")
devtools::load_all("/gpfs/commons/groups/imielinski_lab/home/sclarke/git/gGnome_dev")
context("skilift")

getPGV = function() {
    public_dir = "~/../sclarke/git/pgv_testing/public/"
    json_file = paste0(public_dir,"datafiles.json")
    datadir = paste0(public_dir,"data/")
    settings = paste0(public_dir,"settings.json")
    higlass.list = list(endpoint = "https://higlass01.nygenome.org/")
    pgv = skilift$new(datafiles_json_path = json_file,datadir = datadir,settings = settings,higlass_metadata = higlass.list)
    return(pgv)
}

# Prepare pairs table
test_pairs = readRDS("/gpfs/commons/home/sclarke/git/pgvdb_test_data/test_pairs.rds")
setkey(test_pairs, pair)
hg19_seq = readRDS("/gpfs/commons/home/sclarke/git/pgvdb_test_data/hg19.seq") #for bigwigs

test_meta_pairs = readRDS("/gpfs/commons/home/sclarke/git/pgvdb_test_data/test.meta_pairs.rds")
setkey(test_meta_pairs, pair)


##################################################################################################################################################################################################################

test_that("genome graphs add correctly using template ", {
    skilift = getPGV()
    nrow1 = nrow(skilift$plots)
    genomes.add = genome_temp(patient_id = test_pairs$pair, x = test_pairs$jabba_gg, annotation = NULL, ref = "hg19", order = NULL)
    skilift$add_plots(genomes.add, cores = 5)
    expect_equal(nrow(skilift$plots), nrow1 + 5)
})

test_that("genome graphs add correctly using template with type = NULL ", {
    skilift = getPGV()
    nrow1 = nrow(skilift$plots)
    ## add without type to make sure genome is added as type
    genomes.add = genome_temp(patient_id = test_pairs$pair, x = test_pairs$jabba_gg, annotation = NULL, ref = "hg19", order = NULL, type = NULL)
    skilift$add_plots(genomes.add, cores = 5)
    expect_equal(nrow(skilift$plots), nrow1 + 5)
})

test_that("genome graphs with annotations (events output) get added correctly", {
    skilift = getPGV()
    nrow1 = nrow(skilift$plots)
    genomes.add = genome_temp(patient_id = test_pairs$pair, x = test_pairs$complex, ref = "hg19", order = NULL)
    skilift$add_plots(genomes.add, cores = 5)
    expect_equal(nrow(skilift$plots), nrow1 + 5)
})

## this also tests arrow_temp because it is within the rebinning function (cov2arrow_pgv)
test_that("coverage plots as arrows get uploaded correctly", {
    skilift = getPGV()
    nrow1 = nrow(skilift$plots)
    add.lst = mclapply(test_pairs$pair, function(pair) {
        add.dt = cov2arrow_pgv(patient.id = pair, dryclean_cov = test_pairs[pair,tumor_dryclean_cov], ref = "hg19")
        return(add.dt)
    }, mc.cores = 5)
    rebin.cov.dt = rbindlist(add.lst)
    ## genomes.add = arrow_temp(patient_id = test_pairs$pair, x = test_pairs$tumor_dryclean_cov, ref = "hg19", order = NULL)
    skilift$add_plots(rebin.cov.dt, cores = 5)
    expect_equal(nrow(skilift$plots), nrow1 + 5)
})

## just doing one for bigwigs, some finicky seqlengths things- in the future implement a way to force the seqlengths to match- would be easier once we do not have _chr references
test_that("coverage plots as bigwigs upload correctly", {
    skilift = getPGV()
    nrow1 = nrow(skilift$plots)
    pair = test_pairs$pair[1]
    cov.gr = readRDS(test_pairs[pair,tumor_dryclean_cov])
    cov.gr2 = as.data.table(cov.gr) %>% GRanges(.,seqlengths = hg19.seq) %>% trim ## only works when trimmming - should implement into upload itself but may make it slower?
    bw.add = bw_temp(patient_id = pair, x = list(cov.gr2), ref = "hg19", order = NULL)
    skilift$add_plots(bw.add, cores = 1)
    expect_equal(nrow(skilift$plots), nrow1 + 1)
})


## need Flow and JaBbA for this at the moment- had to keep this as the original directory and not use the test location to get coverage from job output
test_that("ppfit plots upload correctly", {
    library(Flow); library(skitools) ## skitools needed for one function (dirr)
    skilift = getPGV()
    nrow1 = nrow(skilift$plots)
    ppfit.add = ppfit_temp(patient_id = test_pairs$pair, x = test_pairs$balanced_gg, ref = "hg19")
    skilift$add_plots(ppfit.add, cores = 5)
    expect_equal(nrow(skilift$plots), nrow1 + 5)
})


## need Flow and JaBbA for this at the moment- had to keep this as the original directory and not use the test location to get coverage from job output
test_that("allelic plots upload correctly", {
    skilift = getPGV()
    nrow1 = nrow(skilift$plots)
    allelic.add = genome_temp(patient_id = test_pairs$pair, x = test_pairs$balanced_gg_rds, ref = "hg19", type = "allelic", annotation = NULL)
    allelic.add = allelic.add[c(1,3:5),] ## I double checked this but for some reason the second test sample balanced_gg_rds is not actually an allelic graph so it fails (it's good that it failed but not for tests)
    skilift$add_plots(allelic.add, cores = 4)
    expect_equal(nrow(skilift$plots), nrow1 + 4)
})



test_that("mutation plots upload correctly", {
    skilift = getPGV()
    nrow1 = nrow(skilift$plots)
    mutations.add = mutations_temp(patient_id = test_pairs$pair,field = "total_copies", x = test_pairs$somatic_snv_cn, ref = "hg19")
    skilift$add_plots(mutations.add, cores = 5)
    expect_equal(nrow(skilift$plots), nrow1 + 5)
})


## one sample tests for meta_data_json and filtered_events_json

test_that("filtered events jsons created", {
    skilift = getPGV() ## just used to get path for output
    pair = test_meta_pairs$pair
    out.file = gsub("settings.json","test_filtered_events.json",skilift$settings)
    filtered_events_json(pair = pair,
                     oncotable = test_meta_pairs[pair,oncotable],
                     jabba_gg = test_meta_pairs[pair,complex],
                     out_file = out.file,
                     cgc_file = "/gpfs/commons/groups/imielinski_lab/DB/COSMIC/v99_GRCh37/cancer_gene_census_fixed.csv",
                     temp_fix = TRUE)
    expect_true(file.exists(out.file))
})


somatic.filtered.vcf = read.delim(test_meta_pairs$sage_somatic_variants, header=F,comment.char='#') %>% as.data.table
sub.vcf = somatic.filtered.vcf[, c("T_GT", "T_ABQ", "T_AD", "VAF_T", "T_DP", "T_RABQ", "T_RAD", "T_RC_CNT","T_RC_IPC","T_RC_JIT", "T_RC_QUAL", "T_RDP","T_SB") := tstrsplit(V11, ":", fixed = TRUE)]

test_that("metadata json created", {
    skilift = getPGV() ## just used to get path for output
    pair = test_meta_pairs$pair
    out.file = gsub("settings.json","test_metadata.json",skilift$settings)
    meta.dt = meta_data_json(pair = pair,
                             out_file = out.file,
                             coverage = test_meta_pairs[pair,decomposed_cov],
                             jabba_gg = test_meta_pairs[pair,complex],
                             svaba_somatic_vcf = test_meta_pairs[pair,svaba_somatic_vcf],
                             seqnames_loh = c(1:22),
                             karyograph = test.meta.pairs[pair,karyograph_rds],
                             strelka2_vcf = test.meta.pairs[pair,strelka2_somatic_filtered_variants],
                             sage_vcf = test.meta.pairs[pair,sage_somatic_variants],
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
    skilift = getPGV() ## just used to get path for output
    pair = test_meta_pairs$pair
    input.folder = "/gpfs/commons/groups/imielinski_lab/home/sclarke/git/skilift_test_data/case_report_data/"
    out.folder = gsub("settings.json","test_distributions_output",skilift$settings)
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


