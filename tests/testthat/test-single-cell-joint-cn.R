## sc_merge_joint_genome_json / sc_export_joint_cn on synthetic graphs (no patient data)
if (!exists("sc_merge_joint_genome_json", mode = "function")) source(file.path(testthat::test_path(), "..", "..", "R", "single-cell-joint-cn.R"))

iv = function(chr, s, e, iid, y) list(chromosome = chr, startPoint = s, endPoint = e, iid = iid, y = y, title = as.character(iid), type = "interval", strand = "*")
cn = function(cid, so, si, type = "ALT") list(cid = cid, source = so, sink = si, title = type, type = type, weight = 1)

test_that("joint chromosomes come from the joint graph, the rest from regular, ids renumbered", {
  joint = list(settings = list(y_axis = list(title = "copy number")), intervals = list(iv("7", 1, 50, 1, 4), iv("7", 51, 100, 2, 2)),
               connections = list(cn(1, 1, -2, "REF")))
  regular = list(settings = list(description = "<h3>JaBbA</h3>"),
                 intervals = list(iv("7", 1, 100, 1, 3), iv("1", 1, 100, 2, 2), iv("1", 101, 200, 3, 1)),
                 connections = list(cn(1, 2, -3, "REF"), cn(2, -1, 2), cn(3, NULL, 3, "LOOSE")))
  m = sc_merge_joint_genome_json(joint, regular, "7")
  expect_equal(length(m$intervals), 4)
  expect_equal(vapply(m$intervals, `[[`, "", "cn_source"), c("joint", "joint", "regular", "regular"))
  expect_equal(vapply(m$intervals, `[[`, 0, "iid"), 1:4)
  ## regular REF 2->-3 becomes 3->-4; the chr7-chr1 junction is dropped; the loose end keeps its kept side
  expect_equal(attr(m, "dropped"), 1L)
  srcs = lapply(m$connections, `[[`, "source"); sinks = vapply(m$connections, function(x) x$sink, 0)
  expect_equal(sinks, c(-2, -4, 4)); expect_null(srcs[[3]])
  expect_equal(vapply(m$connections, `[[`, 0, "cid"), 1:3)
  expect_equal(m$settings$description, "<h3>Joint junction balance</h3>")
})

test_that("export writes per-cell complex.json and a manifest with the verdict", {
  skip_if_not_installed("gGnome")
  suppressMessages(library(gGnome))
  d = tempfile(); dir.create(d)
  gr = GenomicRanges::GRanges(c("7", "7"), IRanges::IRanges(c(1, 1e6 + 1), c(1e6, 2e6)), cn = c(2, 5))
  gg = gGnome::gG(nodes = gr); gg$nodes$mark(cn = c(2, 5))
  saveRDS(list(gg = gg), file.path(d, "c1.rds"))
  reg = list(settings = list(), intervals = list(iv("1", 1, 100, 1, 2)), connections = list())
  jsonlite::write_json(reg, file.path(d, "c1.json"), auto_unbox = TRUE)
  v = list(much_better = FALSE, default_joint = FALSE, criteria = list(list(key = "c2_fit", label = "fit", pass = FALSE, detail = "x")))
  man = suppressWarnings(sc_export_joint_cn("P", c("c1", "c2"), file.path(d, c("c1.rds", "nope.rds")), file.path(d, c("c1.json", "c1.json")),
                           file.path(d, "out"), "7", v))
  expect_equal(man$default_source, "regular")
  expect_equal(as.character(man$cells), "c1")
  expect_equal(man$cells_missing[[1]]$cell, "c2")
  j = jsonlite::fromJSON(file.path(d, "out", "cells", "c1", "complex.json"), simplifyVector = FALSE)
  expect_true(all(c("joint", "regular") %in% vapply(j$intervals, `[[`, "", "cn_source")))
  mj = jsonlite::fromJSON(file.path(d, "out", "manifest.json"))
  expect_equal(mj$format, "gos-sc-cn-joint/1"); expect_equal(mj$cells, "c1")
})
