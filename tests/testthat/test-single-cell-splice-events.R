suppressWarnings(devtools::load_all())
library(testthat)
library(data.table)

## synthetic gene G (+ strand): exons 100-200, 300-400, 500-600; introns 201-299, 401-499
synthetic_splice <- function() {
    set.seed(3)
    ex <- data.table(chromosome = "1", start = c(100L, 300L, 500L), end = c(200L, 400L, 600L), strand = "+", gene = "G", transcript = "T1",
                     exon = 1:3, gene_type = "protein_coding")
    introns <- data.table(chromosome = "1", start = c(201L, 401L), end = c(299L, 499L), strand = "+", gene = "G", transcript = "T1",
                          exon_left = 1:2, exon_right = 2:3)
    ## annotated 201-299 and a novel acceptor 201-305 (6 bp inside exon 2, frameshift 6 %% 3 == 0 -> in frame), used by clone A
    cells <- sprintf("c%02d", 1:30)
    clone <- setNames(rep(c("A", "B"), each = 15), cells)
    ann <- ifelse(clone == "A", 3L, 40L); nov <- ifelse(clone == "A", 40L, 0L)
    x <- rbind(data.table(rna_id = cells, chromosome = "1", start = 201L, end = 299L, count = ann),
               data.table(rna_id = cells, chromosome = "1", start = 201L, end = 305L, count = nov))[count > 0]
    x[, `:=`(cluster_id = "clu_1_201_305_+", cl_total = sum(count)), by = rna_id]
    cl <- data.table(chromosome = "1", start = c(201L, 201L), end = c(299L, 305L), strand = "+", gene = "G", annotated = c(TRUE, FALSE),
                     type = c("annotated", "novel_acceptor"), n_skipped = 0L, cluster_id = "clu_1_201_305_+")
    allj <- rbind(x[, .(reads = sum(count), cells = .N, patient = "P1"), by = .(chromosome, start, end)],
                  data.table(chromosome = "1", start = c(201L, 201L), end = c(299L, 305L), reads = c(500L, 0L), cells = c(50L, 0L), patient = "P2"))[reads > 0]
    list(ex = ex, introns = introns, x = x, cl = cl, allj = allj, clone = clone, cells = cells)
}

test_that("site offsets and regions follow the transcript direction", {
    expect_equal(sc_splice_site_offset(201L, 201L, 299L, "+", "donor"), 1L)
    expect_equal(sc_splice_site_offset(299L, 201L, 299L, "+", "acceptor"), 1L)
    expect_equal(sc_splice_site_offset(299L, 201L, 299L, "-", "donor"), 1L)
    expect_equal(sc_splice_site_offset(201L, 201L, 299L, "-", "acceptor"), 1L)
    expect_equal(sc_splice_site_region(c(1L, 2L, 3L, 6L, 7L, 0L, -2L), "donor"), c("canonical", "canonical", "motif", "motif", "near", "motif", "motif"))
    expect_equal(sc_splice_site_region(c(1L, 3L, 20L, 21L, 0L), "acceptor"), c("canonical", "motif", "motif", "near", "motif"))
})

test_that("events summarise the novel acceptor with PSI, cells and clone enrichment", {
    s <- synthetic_splice()
    tx_rank <- c(T1 = 1)
    r <- sc_splice_events("P1", s$x, s$cl, s$allj, s$introns, s$ex, tx_rank, coding = "G", cancer = "G", cg = s$clone, n_boot = 50L)
    e <- r$events[type == "novel_acceptor"]
    expect_equal(nrow(e), 1L)
    expect_equal(e$event_type, "alt_3")
    expect_equal(e$cells_covering, 30L)
    expect_equal(e$cells_using, 15L)
    expect_equal(e$reads, 600L)
    expect_equal(e$psi, 600 / (600 + 15 * 3 + 15 * 40))
    expect_true(e$psi_lo <= e$psi && e$psi <= e$psi_hi)
    expect_equal(e$clone, "A")
    expect_true(e$clone_p < 0.001)
    expect_equal(e$frame, "in_frame")
    expect_match(e$summary, "^Novel acceptor 6 bp inside exon 2 of G; used by \\d+% of reads in 15 of 30 cells; in frame$")
    expect_true(e$recommended)
    expect_match(e$reason, "novel in cancer gene")
    ## a cancer-gene-less call is not recommended without another reason
    r2 <- sc_splice_events("P1", s$x, s$cl, s$allj, s$introns, s$ex, tx_rank, coding = "G", cancer = "OTHER", cg = NULL, n_boot = 20L)
    expect_false(r2$events[type == "novel_acceptor"]$recommended)
    ## the same junction at a similar PSI in >= 3 other patients is flagged common
    expect_equal(r$cells[[e$event_id]][, sum(k)], 600L)
})

test_that("SNV test compares carriers and non-carriers", {
    k <- c(rep(10L, 8), rep(0L, 8)); n <- rep(20L, 16)
    t1 <- sc_splice_snv_test(k, n, c(rep(TRUE, 8), rep(FALSE, 8)))
    expect_equal(t1$delta_psi, 0.5)
    expect_lt(t1$p, 0.01)
    expect_null(sc_splice_snv_test(k, n, c(TRUE, rep(FALSE, 15))))
})

test_that("SNV links flag canonical splice-site positions and pair carriers with usage", {
    s <- synthetic_splice(); tx_rank <- c(T1 = 1)
    r <- sc_splice_events("P1", s$x, s$cl, s$allj, s$introns, s$ex, tx_rank, coding = "G", cancer = "G", cg = s$clone, n_boot = 20L)
    ## variants: 0 = donor +1 (canonical), 1 = donor +4 (extended), 2 = 5 kb away
    vid <- c("chr1_201_G_T", "chr1_204_A_C", "chr1_5000_A_C")
    cells <- setNames(lapply(s$cells, function(cid) {
        carrier <- s$clone[[cid]] == "A"
        list(list(0L, if (carrier) 0L else 20L, if (carrier) 10L else 0L, NULL), list(1L, 12L, 0L, NULL), list(2L, 5L, 0L, NULL))
    }), s$cells)
    f <- tempfile(fileext = ".json")
    jsonlite::write_json(list(schemaVersion = 1, format = "compact", variants = lapply(vid, function(i) list(id = i, category = "private", gene = "G")), cells = cells),
                         f, auto_unbox = TRUE, null = "null")
    cm <- setNames(s$cells, s$cells)
    l_all <- sc_splice_snv_links(r$events, r$cells, f, cm, cg = s$clone)
    l <- l_all[event_id == r$events[type == "novel_acceptor"]$event_id]
    expect_equal(sort(l$id), sort(c("chr1_201_G_T", "chr1_204_A_C")) )
    a <- l[id == "chr1_201_G_T"]
    expect_equal(a$region, "canonical"); expect_true(a$canonical); expect_equal(a$site, "donor")
    expect_equal(a$n_carriers, 15L); expect_equal(a$n_noncarriers, 15L)
    expect_gt(a$delta_psi, 0.3)
    expect_lt(a$p, 0.001)
    expect_equal(a$carrier_clone, "A")
    expect_equal(l[id == "chr1_204_A_C"]$region, "motif")
    j <- sc_splice_events_json(r, l_all)
    expect_true(any(vapply(j, function(e) length(e$snv_links) > 0, TRUE)))
    expect_true(all(c("summary", "psi_lo", "cells_covering") %in% names(j[[1]])))
})

test_that("the segment between a nearby new site and the annotated site is searched, a distant one is not", {
    ev <- data.table(event_id = c("near", "far"), chromosome = "1", strand = "+", start = 201L, end = c(305L, 5299L), ref_start = 201L, ref_end = 299L)
    v <- data.table(index = 0:1, chromosome = "1", position = c(302L, 3000L))
    h <- sc_splice_snv_hits(ev, v, window = 0L)
    expect_true(any(h$event_id == "near" & h$index == 0L))
    expect_false(any(h$event_id == "far" & h$index == 1L))
})

test_that("repeat flag only tests novel splice sites", {
    track <- GenomicRanges::GRanges("1", IRanges::IRanges(c(1000, 5000), c(1030, 5040)))
    ev <- data.table::data.table(chromosome = "1", start = c(1010L, 1010L, 3000L), end = c(9000L, 9000L, 5020L))
    ## 1st: novel start in the repeat; 2nd: the same start is annotated; 3rd: novel end in the repeat
    f1 <- sc_splice_repeat_flag(ev[c(1, 3)], track, ann_starts = character(0), ann_ends = "1:9000")
    expect_equal(f1, c(TRUE, TRUE))
    f2 <- sc_splice_repeat_flag(ev[2], track, ann_starts = "1:1010", ann_ends = "1:9000")
    expect_false(f2)
})

test_that("normal-panel flag needs reads, cells and a comparable PSI", {
    expect_true(sc_splice_normal_flag(0.4, 0.3, 5L, 3L))
    expect_false(sc_splice_normal_flag(0.4, 0.05, 5L, 3L))
    expect_false(sc_splice_normal_flag(0.4, 0.4, 5L, 1L))
    expect_false(sc_splice_normal_flag(0.4, NA_real_, NA_integer_, NA_integer_))
})

test_that("read cells: most supporting cells first, then controls that cover the cluster without the junction", {
    ev <- data.table(event_id = c("e1", "e2"), chromosome = "1", start = c(201L, 5000L), end = c(305L, 5100L))
    cd <- data.table(rna_id = sprintf("c%02d", 1:14), k = c(0L, 0L, 0L, 12L, 9L, 8L, 7L, 6L, 5L, 4L, 3L, 2L, 1L, 0L),
                     n = c(50L, 90L, 70L, rep(20L, 10), 0L))
    evr <- list(events = ev, cells = list(e1 = cd, e2 = cd[0]))
    rc <- sc_splice_read_cells(evr, n_use = 8L, n_control = 2L)
    expect_equal(rc$e1[role == "using"]$rna_id, sprintf("c%02d", 4:11))
    expect_equal(rc$e1[role == "control"]$rna_id, c("c02", "c03"))
    expect_null(rc$e2)
    rg <- sc_splice_event_regions(ev, rc)
    expect_length(rg, 1L)
    expect_equal(rg[[1]]$regions, "chr1:1-606")
    expect_equal(length(rg[[1]]$rna_id), 10L)    ## no cell without reads on the event: the lowest-share cells are the comparators
    cd2 <- data.table(rna_id = c("a", "b", "c", "d"), k = c(10L, 9L, 1L, 2L), n = c(10L, 10L, 10L, 10L))
    rc2 <- sc_splice_read_cells(list(events = ev[1], cells = list(e1 = cd2)), n_use = 2L, n_control = 2L)
    expect_equal(rc2$e1[role == "control"]$rna_id, c("c", "d"))
})
