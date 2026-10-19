test_that("sc_neftel_quadrant: halves, signs and arg-max state (synthetic scores)", {
  q <- sc_neftel_quadrant(MES = c(0.1, 0.9, -0.5, 0), AC = c(0.2, 0.1, 0.4, 0), OPC = c(0.8, 0.1, -1, 0.2), NPC = c(0.3, -0.2, -1, 0.1))
  expect_equal(q$state, c("OPC", "MES", "AC", "OPC"))
  # cell 1: D = 0.8 - 0.2 = 0.6; OPC > NPC -> x = -log2(0.5 + 1)
  expect_equal(q$D[1], 0.6)
  expect_equal(q$x[1], -log2(1.5))
  # cell 2: lower half (D = 0.1 - 0.9 < 0), MES > AC -> positive
  expect_lt(q$D[2], 0)
  expect_equal(q$x[2], log2(0.8 + 1))
  # cell 3: AC > MES in the lower half -> negative x
  expect_lt(q$x[3], 0)
})

test_that("sc_neftel_quadrant: missing scores give NA state", {
  q <- sc_neftel_quadrant(NA_real_, 0, 0, 0)
  expect_true(is.na(q$state))
})

test_that("sc_export_state_scores writes columns, state == arg-max (synthetic Seurat-like object)", {
  skip_if_not_installed("jsonlite")
  tmp <- tempfile("scores")
  dir.create(file.path(tmp, "_cohort", "rna"), recursive = TRUE)
  meta <- data.frame(MES1 = c(0.5, 0.1, NA), AC2 = c(0, 0.6, NA), OPC3 = c(0.1, 0.2, NA), NPC4 = c(0.2, 0.3, NA),
                     state = c("MES", "AC", NA), row.names = c("c1", "c2", "c3"), stringsAsFactors = FALSE)
  setClass("FakeSeuratForTest", representation(meta.data = "data.frame"))
  rds <- file.path(tmp, "P1.rds")
  saveRDS(new("FakeSeuratForTest", meta.data = meta), rds)
  out <- suppressMessages(sc_export_state_scores(c(P1 = rds), tmp))
  expect_equal(out$state_equals_argmax_fraction, 1)
  j <- jsonlite::read_json(file.path(tmp, "_cohort", "rna", "state_scores.json"), simplifyVector = TRUE)
  expect_equal(nrow(j$cells), 2)
  expect_equal(j$cells$state, c("MES", "AC"))
})
