## Synthetic tests for the per-cell allelic CN / het-SNP reconciliation (no patient data).
library(data.table)

sim_sites <- function(chr, start, n, f, r = 0.1, spacing = 1000L, lambda = 3, seed = 1) {
    set.seed(seed)
    d <- 2L + rpois(n, lambda)
    drop <- runif(n) < r
    ## minor-copy allele fraction for amplified sites, dropout sites show one allele
    p <- ifelse(drop, ifelse(runif(n) < 0.5, 0.005, 0.995), ifelse(runif(n) < 0.5, f, 1 - f))
    p <- pmin(pmax(p, 0.005), 0.995)
    alt <- rbinom(n, d, p)
    data.table(chr = chr, pos = as.integer(start + seq_len(n) * spacing), ref = d - alt, alt = alt)
}

test_that("sc_ado_loglik separates LOH from balanced sites", {
    h <- sim_sites("1", 0, 2000, f = 0.5)
    k <- pmin(h$ref, h$alt); d <- h$ref + h$alt
    expect_gt(sc_ado_loglik(k, d, f = 0.5), sc_ado_loglik(k, d, f = 0))
    l <- sim_sites("1", 0, 2000, f = 0)
    k <- pmin(l$ref, l$alt); d <- l$ref + l$alt
    expect_gt(sc_ado_loglik(k, d, f = 0), sc_ado_loglik(k, d, f = 0.5))
})

test_that("sc_fit_ado recovers a dropout rate", {
    h <- sim_sites("1", 0, 20000, f = 0.5, r = 0.3, seed = 3)
    kd <- data.table(k = pmin(h$ref, h$alt), d = h$ref + h$alt)[, .(n = .N), by = .(k, d)]
    fit <- sc_fit_ado(kd, rho_grid = 0, w = 0)
    expect_equal(fit$r, 0.3, tolerance = 0.06)
})

test_that("sc_allelic_baf_fix rejects false LOH, keeps true LOH and calls CN-LOH", {
    segs <- data.table(seg_id = 1:8,
                       chr = c("1", "19", "10", "7", "4", "5", "X", "13"),
                       start = 1L, end = 2e8L,
                       major = c(1L, 1L, 1L, 3L, 1L, 1L, 1L, 1L),
                       minor = c(1L, 0L, 0L, 0L, 1L, 0L, 0L, 0L))
    hets <- rbind(
        sim_sites("1", 0, 6000, f = 0.5, seed = 11),      # 1:1, calibration
        sim_sites("19", 0, 3000, f = 0.5, seed = 12),     # called 1:0, truly balanced (under-called total)
        sim_sites("10", 0, 3000, f = 0, seed = 13),       # true LOH
        sim_sites("7", 0, 4000, f = 1/3, seed = 14),      # called 3:0, truly 2:1
        sim_sites("4", 0, 3000, f = 0, seed = 15),        # called 1:1, truly CN-LOH
        sim_sites("5", 0, 10, f = 0.5, seed = 16),        # too few sites
        sim_sites("X", 3e6, 3000, f = 0, seed = 17),      # male X (sites from artefacts only)
        sim_sites("13", 0, 1500, f = 0, seed = 18),       # half LOH ...
        sim_sites("13", 2e6, 1500, f = 0.5, seed = 19))   # ... half het: spans a breakpoint
    res <- sc_allelic_baf_fix(segs, hets, sex = "male")
    s <- res$segs[order(seg_id)]
    expect_equal(s$status, c("het_supported", "loh_rejected_total_raised", "loh_supported", "loh_rejected",
                             "loh_by_baf", "few_sites", "hemizygous_x", "mixed_breakpoint"))
    expect_equal(s[seg_id == 2, c(new_major, new_minor)], c(1L, 1L))
    expect_equal(s[seg_id == 4, c(new_major, new_minor)], c(2L, 1L))
    expect_equal(s[seg_id == 5, c(new_major, new_minor)], c(2L, 0L))
    expect_equal(s[seg_id == 8, c(new_major, new_minor)], c(1L, 0L))
    ## fix_het = FALSE leaves het calls alone
    res2 <- sc_allelic_baf_fix(segs, hets, sex = "male", fix_het = FALSE)
    expect_equal(res2$segs[seg_id == 5, c(status, new_minor)], c("loh_by_baf_not_applied", "1"))
})

test_that("sex = auto: unknown without chrX sites, female with dense chrX hets", {
    segs <- data.table(seg_id = 1:2, chr = c("1", "X"), start = 1L, end = 2e8L, major = 1L, minor = c(1L, 0L))
    h1 <- sim_sites("1", 0, 5000, f = 0.5, seed = 21)
    expect_equal(sc_allelic_baf_fix(segs, h1)$params$sex, "unknown")
    hx <- sim_sites("X", 3e6, 3000, f = 0.5, seed = 22)
    res <- sc_allelic_baf_fix(segs, rbind(h1, hx))
    expect_equal(res$params$sex, "female")
    expect_equal(res$segs[seg_id == 2, status], "loh_rejected_total_raised")
})
