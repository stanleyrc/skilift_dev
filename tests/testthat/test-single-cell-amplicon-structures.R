library(data.table)

test_that("sc_graph_states: identical junction sets share a state, distant sets split (synthetic)", {
  x <- data.table(pair = c("a", "a", "b", "b", "c", "c", "d"), coord = c("j1", "j2", "j1", "j2", "j3", "j4", "j1"))
  s <- sc_graph_states(x, 0.25)
  expect_equal(s[pair == "a", graph_state], s[pair == "b", graph_state])
  expect_false(s[pair == "a", graph_state] == s[pair == "c", graph_state])
  # d shares 1 of 2 junctions with a (Jaccard distance 0.5 > 0.25): its own state
  expect_false(s[pair == "d", graph_state] == s[pair == "a", graph_state])
  expect_equal(sc_graph_states(data.table(pair = "a", coord = "j1"))$graph_state, 1L)
})

jj_fixture <- function() {
  data.table(
    patient = "P1",
    pair = c("c1", "c1", "c2", "c2", "c3", "c3", "c4", "c5"),
    coord = c("j1", "j2", "j1", "j2", "j3", "j9", "j1", "j8"),
    cn = c(5, 3, 4, 2, 6, 1, 5, 2),
    span = c(2e6, 3e6, 2e6, 3e6, 4e6, 5e4, 2e6, 1e6),
    genes_in_cluster = c("EGFR", "EGFR", "EGFR", "EGFR", "EGFR", "EGFR", "EGFR", "Not_in_subgraph"),
    alt_junctions = c("j1", "j2", "j1", "j2", "j3", "j9", "j1", "j8"),
    type = "ALT")
}

test_that("sc_amplicon_structure_input drops short junctions and keeps junction-less cells as one placeholder", {
  jj <- jj_fixture()
  jj <- rbind(jj, data.table(patient = "P1", pair = "c6", coord = "j7", cn = 1, span = 10, genes_in_cluster = "EGFR", alt_junctions = "j7", type = "ALT"))
  x <- sc_amplicon_structure_input(jj, "P1", min_span = 1e5)
  expect_false("j9" %in% x$coord)
  expect_true(x[pair == "c6", coord] == "No_ALT_junctions")
})

test_that("sc_amplicon_structures: PASS classes with >= min_cells, structures per class, per-cell ids", {
  clusters <- data.table(pair = c("c1", "c2", "c3", "c4", "c1", "c2"), genes_in_cluster = c("EGFR", "EGFR", "EGFR", "EGFR", "MDM2,CDK4", "MDM2,CDK4"),
                         mean_cn = c(40, 50, 60, 30, 20, 12), filter = c("PASS", "PASS", "PASS", "PASS", "PASS", "FAIL"))
  a <- sc_amplicon_structures(clusters, jj_fixture(), "P1", min_cells = 3)
  expect_equal(a$classes$cls, "EGFR")
  expect_equal(a$classes$n_cells, 4L)
  expect_equal(a$dropped$cls, "CDK4+MDM2")
  ids <- setNames(a$cells$structure, a$cells$cell)
  expect_equal(ids[["c1"]], ids[["c2"]])
  expect_false(ids[["c1"]] == ids[["c3"]])
  expect_equal(a$classes$n_structures, length(unique(ids)))
  # duplicate clusters of one class in a cell: one entry, CN summed
  dup <- rbind(clusters, data.table(pair = "c1", genes_in_cluster = "EGFR", mean_cn = 5, filter = "PASS"))
  b <- sc_amplicon_structures(dup, jj_fixture(), "P1", min_cells = 3)
  expect_equal(nrow(b$cells[cls == "EGFR"]), 4L)
  expect_equal(b$cells[cell == "c1" & cls == "EGFR", cn], 45)
})

test_that("sc_export_amplicon_structures writes the cohort JSON (synthetic rds inputs)", {
  skip_if_not_installed("jsonlite")
  tmp <- tempfile("amp"); dir.create(tmp)
  clusters <- data.table(pair = c("c1", "c2", "c3", "c4"), genes_in_cluster = "EGFR", mean_cn = c(40, 50, 60, 30), filter = "PASS")
  f1 <- file.path(tmp, "gc.rds"); saveRDS(clusters, f1)
  f2 <- file.path(tmp, "jj.rds"); saveRDS(list(jj.dt = jj_fixture()), f2)
  db <- data.table(patient = c("P1", "P2"), gene_clusters_amp_thresh4 = c(f1, NA), gene_cluster_sets_thresh4_lst = c(f2, NA))
  sc_export_amplicon_structures(db, tmp, n_tumour = c(P1 = 10), min_cells = 3)
  j <- jsonlite::read_json(file.path(tmp, "_cohort", "amplicon_structures.json"))
  expect_equal(j$format, "gos-sc-amplicon-structures/1")
  expect_equal(j$patients$P1$n_tumour_cells, 10)
  expect_equal(j$patients$P1$classes[[1]]$cls, "EGFR")
  expect_length(j$patients$P1$classes[[1]]$structure, 4)
  expect_equal(j$not_available$P2, "no gene-cluster table")
})

test_that("sc_ecdna_snv_sites / sc_export_ecdna_snvs keep sites with no normal alt reads (synthetic)", {
  pile <- data.table(cell = c("a", "b", "a", "b", "a"), site = c("100_A_T", "100_A_T", "200_G_C", "200_G_C", "300_C_A"),
                     ref_count = c(3, 2, 1, 4, 0), alt_count = c(2, 1, 1, 1, 3))
  normal <- data.table(pos = c(100L, 200L, 300L), normal_alt = c(0L, 2L, 0L), normal_depth = c(8L, 9L, 9L))
  expect_equal(sc_ecdna_snv_sites(pile, normal), "100_A_T")
  tmp <- tempfile("esnv"); dir.create(tmp)
  sc_export_ecdna_snvs("P1", pile, normal, tmp, chromosome = "7")
  j <- jsonlite::read_json(file.path(tmp, "P1", "ecdna_snvs.json"))
  expect_equal(j$n_sites, 1)
  expect_equal(j$sites[[1]]$id, "7_100_A_T")
  expect_equal(unlist(j$sites[[1]]$alt_count), c(2, 1))
})
