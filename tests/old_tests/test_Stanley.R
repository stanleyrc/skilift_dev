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
    skilift <- reset_skilift()
    pair = test_meta_pairs$pair
    out.file = gsub("settings.json","test_metadata.json",skilift$settings)
    meta.dt = meta_data_json(pair = pair,
                             out_file = out.file,
                             coverage = test_meta_pairs[pair,decomposed_cov],
                             jabba_gg = test_meta_pairs[pair,complex],
                             svaba_somatic_vcf = test_meta_pairs[pair,svaba_somatic_vcf],
                             seqnames_loh = c(1:22),
                             karyograph = test_meta_pairs[pair,karyograph_rds],
                             strelka2_vcf = test_meta_pairs[pair,strelka2_somatic_filtered_variants],
                             sage_vcf = test_meta_pairs[pair,sage_somatic_variants],
                             tumor_type = test_meta_pairs[pair,tumor_type_final],
                             disease = test_meta_pairs[pair,disease],
                             primary_site = test_meta_pairs[pair,primary_site_simple],
                             inferred_sex = test_meta_pairs[pair,inferred_sex],
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


