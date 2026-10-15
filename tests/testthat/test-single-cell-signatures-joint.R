## synthetic signatures only: columns are made-up smooth profiles, not COSMIC
make_S <- function() {
    set.seed(11)
    m <- vapply(1:8, function(i) { x <- rgamma(96, 0.3 + 0.2 * i)^2; x / sum(x) }, numeric(96))
    dimnames(m) <- list(SC_SBS96, paste0("S", 1:8))
    m
}
draw <- function(S, w, n) {
    p <- as.vector(base::`%*%`(S, w)); p <- p / sum(p)
    ctx <- sample(SC_SBS96, n, TRUE, prob = p)
    ctx
}

test_that("NNLS refit recovers exposures and bootstrap CI covers truth", {
    skip_if_not_installed("nnls")
    S <- make_S(); set.seed(1)
    w <- c(S1 = .5, S2 = .3, S3 = .2)
    fit <- sc_sig_refit_bootstrap(as.numeric(table(factor(draw(S, c(.5, .3, .2, 0, 0, 0, 0, 0), 3000), SC_SBS96))), S, c("S1", "S2", "S3"), n_boot = 200)
    expect_equal(unname(round(fit$fraction, 1)), unname(w), tolerance = 0.06)
    expect_true(all(fit$fraction >= apply(fit$boot_fraction, 2, quantile, 0.005) & fit$fraction <= apply(fit$boot_fraction, 2, quantile, 0.995)))
})

test_that("a real subset-only signature is added, a null subset adds nothing", {
    skip_if_not_installed("nnls")
    S <- make_S(); set.seed(2)
    sub_real <- as.numeric(table(factor(draw(S, c(.4, .2, 0, 0, .4, 0, 0, 0), 1500), SC_SBS96)))
    sub_null <- as.numeric(table(factor(draw(S, c(.5, .5, 0, 0, 0, 0, 0, 0), 1500), SC_SBS96)))
    ps <- sc_sig_patient_set(list(real = sub_real, null = sub_null), S, c("S1", "S2"), n_boot = 40)
    expect_true("S5" %in% ps$signatures)
    expect_equal(ps$reasons[signature == "S5"]$set, "real")
    expect_setequal(ps$signatures, c("S1", "S2", "S5"))
})

test_that("joint refit finds a subclonal-only exposure difference", {
    skip_if_not_installed("nnls")
    S <- make_S(); set.seed(3)
    sets <- list(truncal = draw(S, c(.6, .4, 0, 0, 0, 0, 0, 0), 800), subclonal = draw(S, c(.4, .3, 0, 0, .3, 0, 0, 0), 800),
                 private = draw(S, c(.5, .5, 0, 0, 0, 0, 0, 0), 800))
    sets$all <- unlist(sets, use.names = FALSE)
    out <- sc_refit_signatures_joint(sets, S, c("S1", "S2"), n_boot = 200, n_boot_select = 30)
    expect_true("S5" %in% vapply(out$signature_set, function(s) s$signature, ""))
    sc <- out$contrasts[[1]]
    expect_equal(sc$contrast, "subclonal vs truncal")
    r5 <- Filter(function(r) r$signature == "S5", sc$rows)[[1]]
    expect_true(r5$diff > 0.15 && r5$lo > 0 && r5$p < 0.05)
    pr <- Filter(function(r) r$signature == "S5", out$contrasts[[2]]$rows)[[1]]
    expect_true(pr$p > 0.05)
    expect_true(is.list(jsonlite::fromJSON(jsonlite::toJSON(out, auto_unbox = TRUE), simplifyVector = FALSE)))
})
