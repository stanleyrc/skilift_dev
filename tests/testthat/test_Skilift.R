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
    paths <<- load_paths()
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

### plot generation methods

# coverage
test_that("coverage plots are generated correctly", {
    pgv <- reset_skilift()
    pairs_path = system.file("extdata", "test_data", "casereport_pairs.rds", package = "Skilift")
    pairs = readRDS(pairs_path)
    cov_path = pairs[1]$decomposed_cov
    cov.add = arrow_temp(
        patient_id = pairs[1,pair],
        x = list(readRDS(cov_path)),
        ref = "hg19",
        field = "foreground",
        title = "Coverage",
        order = 50,
        overwrite = TRUE
    )

    create_cov_arrow(plot_metadata = cov.add, datadir = paths$datadir)
})

# sage
test_that("sage methods work correctly", {
  pgv <- reset_skilift()
  test_meta_pairs = readRDS("/gpfs/commons/home/sclarke/git/pgvdb_test_data/test.meta_pairs.rds")
  sage_vcf_gz = test_meta_pairs[1]$sage_somatic_variants
  sage_dt = sage_qc(sage_vcf_gz, "hg19", write_json = FALSE, return_table = TRUE)
  sage_count = sage_count(sage_vcf_gz, "hg19")

  expect_true(!is.null(sage_dt) && length(sage_dt) > 0)
  expect_true(!is.null(sage_count) && length(sage_count) > 0)
})

# hetsnps
test_that("hetsnps.arrow generation works correctly", {
    pgv <- reset_skilift()
    pairs_path = system.file("extdata", "test_data", "casereport_pairs.rds", package = "Skilift")
    pairs = readRDS(pairs_path)
    het_pileups_wgs_path = pairs[1]$het_pileups_wgs
    subset.hets.gr = subsample_hetsnps(het_pileups_wgs=het_pileups_wgs_path)
    genome.add = arrow_temp(
        patient_id = pairs[1,pair],
        x = list(subset.hets.gr),
        ref = "hg19",
        field = "count",
        title = "HET SNPS SUBSETTED",
        order = 50
    )

    genome.add$source = "hetsnps.arrow"

    create_cov_arrow(plot_metadata = genome, datadir = paths$datadir)
})

# ggraph 
test_that("ggraph plots are generated correctly", {
    pgv <- reset_skilift()
    pairs_path = system.file("extdata", "test_data", "casereport_pairs.rds", package = "Skilift")
    pairs = readRDS(pairs_path)
    ggraph_path = pairs[1]$jabba_gg
    ggraph.add = genome_temp(
        patient_id = pairs[1,pair],
        x = ggraph_path,
        ref = "hg19",
        source = "genome.json", # for complex/events output, use source=complex.json
        title = "ggraph",
        order = 50,
        overwrite = TRUE
    )

    create_ggraph_json(plot_metadata = ggraph.add, datadir = paths$datadir)
})

# snv multiplicity
test_that("multiplicity plots are generated correctly", {
    pgv <- reset_skilift()
    pairs_path = system.file("extdata", "test_data", "test_pairs.rds", package = "Skilift")
    pairs = readRDS(pairs_path)
    multiplicity_path = pairs[1]$somatic_snv_cn
    mutations.add = mutations_temp(
        patient_id = pairs[1]$pair,
        field = "total_copies",
        x = multiplicity_path,
        ref = "hg19",
        overwrite = TRUE
    )

    create_somatic_json(plot_metadata = mutations.add, datadir = paths$datadir)
})

# filtered events
# WIP since oncotable will be ported to skilift
# test_that("filtered events table is generated correctly", {
#     oncotable_paths = mclapply(filtered_outputs$patient_id, function(pair) {
#       output_row = outputs[patient_id == pair, ]
#       # construct outdir using sub on jabba_simple path
#       outdir = sub('jabba', 'oncotable', output_row$jabba_simple)
#       outdir = sub('/jabba.simple.rds', '', output_row$jabba_simple)
#       paths_row = data.table(patient_id = pair, oncotable = paste0(outdir, "/oncotable.rds"))
#       create_oncotable(
#           pair = pair,
#           annotated_bcf = output_row$snpeff_bcf,
#           signature_counts = NULL,
#           jabba_simple = output_row$jabba_simple,
#           events = output_row$events,
#           fusions = output_row$fusions,
#           gencode = NULL,
#           amp_thresh_multiplier = 2,
#           outdir = outdir
#       )
#
#       return(paths_row)
#     }, mc.cores = 5)
#     oncotable_paths = rbindlist(oncotable_paths)
#
#     mclapply(filtered_outputs$patient_id, function(pair) {
#       cgc_file = "~/DB/COSMIC/v99_GRCh37/cancer_gene_census_fixed.csv"
#       oncotable_path = oncotable_paths[patient_id == pair, oncotable]
#       outfile = paste0(datadir, "/", pair, "/filtered.events.json")
#       filtered_events_json(
#         pair = pair,
#         oncotable = oncotable_path,
#         jabba_gg = outputs[patient_id == pair, events],
#         out_file = outfile,
#         cgc_file = cgc_file,
#         temp_fix = TRUE
#       )
#     }, mc.cores = 5)
# })

# mutation catalog
test_that("mutation catalog generation works correctly", {
    pgv <- reset_skilift()
    indel_matrix = "~/../spanja/Projects/casereports/clinical_ip/output/ID/clinical_ip.ID28.all"
    create_mutations_catalog_json(
        sig_matrix_path=indel_matrix,
        is_indel=TRUE,
        output_dir="./test_mut_catalog/"
    )
    sbs_matrix = "~/../spanja/Projects/casereports/clinical_ip/output/SBS/clinical_ip.SBS96.all"
    create_mutations_catalog_json(
        sig_matrix_path=sbs_matrix,
        is_indel=TRUE,
        output_dir="./test_mut_catalog/"
    )

    indel_matrix_test = "~/Projects/nf-casereports/tests/test_runs/chr21_test_results/signatures/sigprofilerassignment/somatic/JTS-1501_T_vs_JTS-1501_N/sig_inputs/output/ID/Input_vcffiles.ID28.all"
    create_mutations_catalog_json(
        sig_matrix_path=indel_matrix_test,
        is_indel=TRUE,
        output_dir="./test_mut_catalog_2/"
    )
    sbs_matrix_test = "~/Projects/nf-casereports/tests/test_runs/chr21_test_results/signatures/sigprofilerassignment/somatic/JTS-1501_T_vs_JTS-1501_N/sig_inputs/output/SBS/Input_vcffiles.SBS96.all"
    create_mutations_catalog_json(
        sig_matrix_path=sbs_matrix_test,
        is_indel=FALSE,
        output_dir="./test_mut_catalog_2/"
    )
})

# metadata.json
test_that("metadata.json is generated correctly", {
    skilift <- reset_skilift()
    pairs_path = system.file("extdata", "test_data", "test_meta_pairs.rds", package = "Skilift")
    pairs = readRDS(pairs_path)
    pair = pairs[1, pair]

    output = system.file("extdata", "test_data", package = "Skilift")
    meta_data_json(
        pair = pair,
        outdir = output,
        coverage = pairs[pair, decomposed_cov],
        jabba_gg = pairs[pair, complex],
        seqnames_loh = c(1:22),
        karyograph = pairs[pair, karyograph_rds],
        strelka2_vcf = pairs[pair, strelka2_somatic_filtered_variants], #optional
        sage_vcf = pairs[pair, sage_somatic_variants],
        tumor_type = pairs[pair, tumor_type_final],
        disease = pairs[pair, disease],
        primary_site = pairs[pair, primary_site_simple],
        seqnames_genome_width = c(1:22, "X", "Y"),
        write_json = TRUE,
        make_dir = TRUE,
        overwrite = TRUE
    )

    # read json file
    md_json = fromJSON(paste0(output, "/", pair, "/metadata.json"))
})

# common distributions
test_that("kpi distributions are generated correctly", {

    input_folder = system.file("extdata", "test_data", "casereports_test_datadir", package = "Skilift")
    output_folder = system.file("extdata", "test_data", "casereports_test_common", package = "Skilift")

    skilift <- reset_skilift()
    create_distributions(
        case_reports_datadir = input_folder,
        common_dir = output_folder,
        filter_patients = NULL,
        write_jsons = TRUE,
        cores = 1
    )

    files = paste0(
        output_folder, "/",
        c(
            "coverageVariance.json",
            "ploidy.json",
            "snvCount.json",
            "tmb.json",
            "lohFraction.json",
            "purity.json",
            "svCount.json"
        )
    )
    expect_true(all(file.exists(files)))
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

### Development
test_that("adding new patients to distributions work correctly", {

    casereport_datadir = system.file("extdata", "test_data", "casereports_test_datadir", package = "Skilift")
    common_dir = system.file("extdata", "test_data", "casereports_test_common", package = "Skilift")

    skilift <- reset_skilift()
    distributions = create_distributions(
        case_reports_datadir = casereport_datadir,
        filter_patients = list("0125"),
        write_to_json = TRUE,
        common_dir = common_dir,
        overwrite = TRUE,
        cores = 1
    )

    metadata = system.file("extdata", "test_data", "casereports_test_datadir", "0124", "metadata.json", package = "Skilift")

    add_patient_to_distributions(
        metadata_json_path = metdata,
        common_dir = common_dir,
        overwrite = TRUE
    )
})
