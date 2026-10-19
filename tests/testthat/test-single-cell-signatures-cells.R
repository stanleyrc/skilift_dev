testthat::test_that("root-to-tip catalogs include ancestors and private sites while excluding other clades", {
    tree <- ape::read.tree(text = "((a:1,b:1):1,c:2);")
    variants <- data.table::data.table(id = paste0("s", 1:5), category = c("truncal", "subclonal", "private", "ambiguous", "outside_tumor"),
        anchor = c("a|b", "a", "c", "a", "a"), context = SC_SBS96[1:5])
    counts <- sc_cell_signature_counts(variants, tree, c("a", "b", "c", "unplaced"))
    testthat::expect_equal(unname(colSums(counts)), c(2, 1, 1, 0))
    testthat::expect_equal(counts[1, c("a", "b", "c")], c(a = 1L, b = 1L, c = 0L))
    testthat::expect_equal(sum(counts[4:5, ]), 0)
})

testthat::test_that("export uses a synthetic counts matrix and flags low counts without fitting them", {
    counts <- matrix(0L, 96, 3, dimnames = list(SC_SBS96, c("a", "b", "small")))
    counts[1:2, "a"] <- c(90L, 30L); counts[1:2, "b"] <- c(30L, 120L); counts[1, "small"] <- 99L
    calls <- NULL
    fitter <- function(sets, work_dir, min_mutations, cfg) {
        calls <<- sets
        data.table::data.table(set = c("a", "a", "b", "b"), signature = c("SBS1", "SBS5", "SBS1", "SBS5"), activity = c(90, 30, 30, 120))
    }
    quality <- function(work_dir) list(sets = list(a = list(stats = list(cosine = 1)), b = list(stats = list(cosine = 1))))
    out_dir <- tempfile("cell_signatures_"); on.exit(unlink(out_dir, recursive = TRUE))
    out <- sc_export_cell_signatures("unused", out_dir = out_dir, counts = counts, fit_fun = fitter, quality_fun = quality, cfg = list())
    testthat::expect_named(calls, c("a", "b"))
    testthat::expect_equal(length(calls$a), 120L)
    testthat::expect_equal(out$cells[[1]]$exposures$SBS1, 0.75)
    testthat::expect_equal(out$cells[[2]]$n, 150)
    testthat::expect_true(out$cells[[3]]$too_few)
    testthat::expect_length(out$cells[[3]]$exposures, 0)
    parsed <- jsonlite::read_json(file.path(out_dir, "signatures_cells.json"))
    testthat::expect_equal(parsed$cells[[1]]$cosine, 1)
    sc_export_cell_signatures("unused", out_dir = out_dir, counts = counts, fit_fun = fitter, quality_fun = quality, cfg = list())
    testthat::expect_length(list.files(out_dir, pattern = "[.]bak_"), 1)
})
