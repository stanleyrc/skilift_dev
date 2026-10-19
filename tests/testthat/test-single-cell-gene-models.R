suppressWarnings(devtools::load_all())
library(testthat)

## synthetic GENCODE-like GRanges: a coding gene AAA (+, two transcripts, one NMD) and a lncRNA BBB (-), no patient data
synthetic_gencode <- function() {
    rows <- data.table::rbindlist(list(
        list("chr7", 100L, 1000L, "+", "gene", "ENSG1.1", "protein_coding", "AAA", NA, NA, NA, NA),
        list("chr7", 100L, 1000L, "+", "transcript", "ENSG1.1", "protein_coding", "AAA", "ENST1.2", "protein_coding", "AAA-201", "basic,CCDS"),
        list("chr7", 100L, 200L, "+", "exon", "ENSG1.1", "protein_coding", "AAA", "ENST1.2", "protein_coding", "AAA-201", "basic,CCDS"),
        list("chr7", 500L, 600L, "+", "exon", "ENSG1.1", "protein_coding", "AAA", "ENST1.2", "protein_coding", "AAA-201", "basic,CCDS"),
        list("chr7", 900L, 1000L, "+", "exon", "ENSG1.1", "protein_coding", "AAA", "ENST1.2", "protein_coding", "AAA-201", "basic,CCDS"),
        list("chr7", 150L, 200L, "+", "CDS", "ENSG1.1", "protein_coding", "AAA", "ENST1.2", "protein_coding", "AAA-201", "basic,CCDS"),
        list("chr7", 500L, 600L, "+", "CDS", "ENSG1.1", "protein_coding", "AAA", "ENST1.2", "protein_coding", "AAA-201", "basic,CCDS"),
        list("chr7", 900L, 947L, "+", "CDS", "ENSG1.1", "protein_coding", "AAA", "ENST1.2", "protein_coding", "AAA-201", "basic,CCDS"),
        list("chr7", 948L, 950L, "+", "stop_codon", "ENSG1.1", "protein_coding", "AAA", "ENST1.2", "protein_coding", "AAA-201", "basic,CCDS"),
        list("chr7", 100L, 1000L, "+", "transcript", "ENSG1.1", "protein_coding", "AAA", "ENST2.1", "nonsense_mediated_decay", "AAA-202", NA),
        list("chr7", 100L, 200L, "+", "exon", "ENSG1.1", "protein_coding", "AAA", "ENST2.1", "nonsense_mediated_decay", "AAA-202", NA),
        list("chr7", 900L, 1000L, "+", "exon", "ENSG1.1", "protein_coding", "AAA", "ENST2.1", "nonsense_mediated_decay", "AAA-202", NA),
        list("chr7", 5000L, 6000L, "-", "gene", "ENSG3.1", "lncRNA", "BBB", NA, NA, NA, NA),
        list("chr7", 5000L, 6000L, "-", "transcript", "ENSG3.1", "lncRNA", "BBB", "ENST3.1", "lncRNA", "BBB-201", "basic"),
        list("chr7", 5000L, 5100L, "-", "exon", "ENSG3.1", "lncRNA", "BBB", "ENST3.1", "lncRNA", "BBB-201", "basic"),
        list("chr7", 5800L, 6000L, "-", "exon", "ENSG3.1", "lncRNA", "BBB", "ENST3.1", "lncRNA", "BBB-201", "basic")))
    data.table::setnames(rows, c("seqnames", "start", "end", "strand", "type", "gene_id", "gene_type", "gene_name",
                                 "transcript_id", "transcript_type", "transcript_name", "tag"))
    GenomicRanges::makeGRangesFromDataFrame(as.data.frame(rows), keep.extra.columns = TRUE)
}

test_that("gene model tables: CDS includes the stop codon, one canonical per gene", {
    m <- sc_gene_models_table(synthetic_gencode())
    expect_equal(sort(m$genes$name), c("AAA", "BBB"))
    expect_equal(m$genes[name == "AAA"]$chromosome, "7")
    t1 <- m$transcripts[transcript_id == "ENST1"]
    expect_equal(c(t1$cds_start, t1$cds_end), c(150L, 950L))
    expect_true(t1$canonical)
    expect_equal(t1$canonical_basis, "inferred")
    expect_false(m$transcripts[transcript_id == "ENST2"]$canonical)
    expect_true(is.na(m$transcripts[transcript_id == "ENST3"]$cds_start))
    expect_equal(m$transcripts[, sum(canonical), by = gene_id]$V1, c(1L, 1L))
    expect_equal(nrow(m$exons[transcript_id == "ENST1"]), 3L)
})

test_that("MANE ids override the inferred canonical; drop_types keeps canonical transcripts", {
    m <- sc_gene_models_table(synthetic_gencode(), mane = "ENST2")
    expect_true(m$transcripts[transcript_id == "ENST2"]$canonical)
    expect_equal(m$transcripts[transcript_id == "ENST2"]$canonical_basis, "MANE_Select")
    d <- sc_gene_models_table(synthetic_gencode(), drop_types = "nonsense_mediated_decay")
    expect_false("ENST2" %in% d$transcripts$transcript_id)
    d2 <- sc_gene_models_table(synthetic_gencode(), mane = "ENST2", drop_types = "nonsense_mediated_decay")
    expect_true("ENST2" %in% d2$transcripts$transcript_id)
})

test_that("sc_export_gene_models writes per-chromosome JSON and an index", {
    out <- tempfile("models")
    idx <- sc_export_gene_models(out, cfg = list(gencode = "synthetic.rds", canonical_gtf = ""), gr = synthetic_gencode())
    expect_equal(idx$chromosome, "7")
    expect_true(file.exists(file.path(out, "chr7.json")))
    j <- jsonlite::fromJSON(file.path(out, "chr7.json"), simplifyVector = FALSE)
    expect_equal(j$chromosome, "7")
    expect_equal(j$n_genes, 2L)
    expect_equal(vapply(j$genes, `[[`, "", "name"), c("AAA", "BBB"))
    aaa <- j$genes[[1]]
    expect_equal(aaa$strand, "+")
    expect_equal(aaa$transcripts[[1]]$id, "ENST1")
    expect_equal(aaa$transcripts[[1]]$canonical, 1L)
    expect_equal(unlist(aaa$transcripts[[1]]$cds), c(150L, 950L))
    expect_equal(length(aaa$transcripts[[1]]$exons), 3L)
    expect_equal(unlist(aaa$transcripts[[1]]$exons[[2]]), c(500L, 600L))
    expect_null(j$genes[[2]]$transcripts[[1]]$cds)
    expect_equal(j$genes[[2]]$type, "lncRNA")
    ix <- jsonlite::fromJSON(file.path(out, "index.json"))
    expect_equal(ix$chromosomes$n_transcripts, 3L)
    ## the size cap drops non-canonical NMD transcripts and records it
    idx2 <- sc_export_gene_models(out, cfg = list(gencode = "synthetic.rds", canonical_gtf = ""), gr = synthetic_gencode(), max_mb = 1e-4)
    expect_match(idx2$dropped_types, "nonsense_mediated_decay")
    expect_equal(idx2$n_transcripts, 2L)
})

test_that("each chromosome file holds only its own genes", {
    gr <- synthetic_gencode()
    gr8 <- gr[gr$gene_name == "BBB"]
    GenomeInfoDb::seqlevels(gr8) <- "chr8"
    gr8$gene_id <- "ENSG8.1"; gr8$gene_name <- "CCC"; gr8$transcript_id[!is.na(gr8$transcript_id)] <- "ENST8.1"
    out <- tempfile("models")
    idx <- sc_export_gene_models(out, cfg = list(gencode = "synthetic.rds", canonical_gtf = ""), gr = suppressWarnings(c(gr, gr8)))
    expect_equal(idx$chromosome, c("7", "8"))
    expect_equal(idx$n_genes, c(2L, 1L))
    j8 <- jsonlite::fromJSON(file.path(out, "chr8.json"), simplifyVector = FALSE)
    expect_equal(vapply(j8$genes, `[[`, "", "name"), "CCC")
    expect_equal(j8$n_transcripts, 1L)
})
