test_that("sc_export_walk_set writes a versioned gos-sc-walks/1 file with method fields (synthetic)", {
    d <- tempfile(); dir.create(d)
    w <- list(list(id = "pv4_sp1_w1", label = "GENEA peel v4 sp1", circular = TRUE, species = "sp1",
                   nodes = data.frame(chromosome = c("7", "7", "7"), start = c(100, 2001, 50000), end = c(2000, 9000, 60000), strand = "+"),
                   cells = c(cellA = 4.2, cellB = 0.2), cells_cplex = c(cellA = 4, cellB = 0),
                   evidence = list(closed = TRUE, n_junctions = 1L), curated_match = list(walk_id = "12", segment_cosine = 0.98, status = "reproduces")))
    p <- sc_export_walk_set(w, d, "PX", method = "test method", version = "v_test", fit = list(R2_free = 0.9))
    expect_equal(basename(p), "walks_peel_v4.json")
    j <- jsonlite::fromJSON(p, simplifyVector = FALSE)
    expect_equal(j$format, "gos-sc-walks/1"); expect_equal(j$version, "v_test")
    x <- j$walks[[1]]
    expect_equal(x$n_nodes, 2)                         # 100-2000 and 2001-9000 collapse into one node
    expect_equal(x$n_nodes_raw, 3)
    expect_equal(vapply(x$junctions, `[[`, "", "type"), c("ALT", "ALT"))
    expect_equal(x$ncells, 1); expect_false(x$curated)
    expect_equal(x$curated_match$status, "reproduces"); expect_equal(x$species, "sp1")
    expect_equal(x$cells_cplex$cellA, 4)
    expect_error(sc_export_walk_set(w, d, "PX", file = "walks.json", method = "m", version = "v"), "never writes walks.json")
})
