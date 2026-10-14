suppressWarnings(devtools::load_all())
library(testthat)

## synthetic oncotable fusion rows, no patient data
fus_rows <- function() {
    data.table::data.table(
        gene = c("ZNF345", "EGFR", "AAA"),
        fusion_genes = c("ZNF345::MLLT3", "EGFR::SEPT14", "AAA::BBB"),
        variant.p = c("Exon -1 (p.5'UTR)::Exon 5 (p.140)", "Exon 24 (p.1000)::Exon 3 (p.20)", "Exon -1 (p.5'UTR)::Exon -1 (p.5'UTR)"),
        type = "fusion", vartype = "fusion", source = "oncokb_fusions", tier = c(1L, 1L, 3L),
        tier_description = c("Clinically Actionable", "Clinically Actionable", "VUS"),
        therapeutics = c("drug", "drug", NA))
}

test_that("fusions landing in the synthetic 5'UTR exon -1 are relabelled upstream and demoted", {
    r <- sc_relabel_upstream_fusions(fus_rows())
    expect_equal(r$fusion_genes, c("upstream of ZNF345::MLLT3", "EGFR::SEPT14", "upstream of AAA::upstream of BBB"))
    expect_equal(r$upstream_partner, c("ZNF345", NA, "AAA,BBB"))
    expect_equal(r$upstream_only, c(TRUE, FALSE, TRUE))
    expect_equal(r$tier, c(3L, 1L, 3L))
    expect_true(is.na(r$therapeutics[1]))
    expect_equal(r$therapeutics[2], "drug")
    expect_equal(nrow(r), 3L)
})

test_that("relabelling is idempotent and ignores CNA rows", {
    once <- sc_relabel_upstream_fusions(fus_rows())
    twice <- sc_relabel_upstream_fusions(data.table::copy(once))
    expect_equal(once$fusion_genes, twice$fusion_genes)
    cna <- data.table::data.table(gene = "EGFR", type = "amp", source = "oncokb_cna", fusion_genes = NA_character_, variant.p = NA_character_, tier = 1L)
    expect_false(sc_relabel_upstream_fusions(cna)$upstream_only %in% TRUE)
    expect_null(sc_relabel_upstream_fusions(NULL))
})
