test_that("sc_export_gene_spans writes sorted gos-gene-spans/1 from a synthetic GENCODE", {
    skip_if_not_installed("GenomicRanges")
    gr <- GenomicRanges::GRanges(c("2", "1", "1"), IRanges::IRanges(c(500, 300, 100), c(900, 400, 200)),
                                 strand = c("-", "+", "+"), type = "gene",
                                 gene_name = c("GENEB", "GENEA", "PSEUDO1"),
                                 gene_type = c("lncRNA", "protein_coding", "unprocessed_pseudogene"))
    d <- tempfile(); dir.create(d)
    rds <- file.path(d, "gencode.test.rds"); saveRDS(gr, rds)
    cfg <- list(gencode = rds, cna_gencode = rds)
    out <- file.path(d, "_cohort", "gene_spans.json")
    res <- sc_export_gene_spans(out, cfg)
    expect_equal(res$n_genes, 2)
    expect_equal(res$n_dropped, 1)
    j <- jsonlite::fromJSON(out, simplifyVector = FALSE)
    expect_equal(j$format, "gos-gene-spans/1")
    expect_equal(j$source, "gencode.test.rds")
    expect_equal(vapply(j$genes, `[[`, "", "gene"), c("GENEA", "GENEB"))
    expect_true(j$genes[[1]]$coding)
    expect_false(j$genes[[2]]$coding)
    expect_equal(j$genes[[1]]$strand, "+")
    expect_equal(j$genes[[1]]$start, 300)
    sc_export_gene_spans(out, cfg, drop_pseudogenes = FALSE)       # second run backs up the first
    expect_true(length(list.files(dirname(out), "bak_")) == 1)
})
