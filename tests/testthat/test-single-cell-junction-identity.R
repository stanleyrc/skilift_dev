library(data.table)
test_that("master ids are ordered by genome position and strand-carrying", {
  bp = data.table(chrom1 = c("chr7", "chr5"), pos1 = c(100L, 5624784L), strand1 = c("-", "+"),
                  chrom2 = c("chr1", "chr5"), pos2 = c(50L, 6302857L), strand2 = c("+", "+"),
                  merged_id = c("JUNC_MERGED_1", "JUNC_MERGED_2"))
  m = sc_junction_master_from_bp(bp)
  expect_equal(m$id, c("1:50+|7:100-", "5:5624784+|5:6302857+"))
  expect_equal(m$chrom1, c("chr1", "chr5"))
})
test_that("support is per cell, max over merged members, zeros dropped, 3 synthetic cells", {
  bp = data.table(chrom1 = "chr5", pos1 = 10L, strand1 = "+", chrom2 = "chr5", pos2 = 20L, strand2 = "-", merged_id = "JUNC_MERGED_1")
  m = sc_junction_master_from_bp(bp)
  counts = data.table(pair = c("c1", "c2", "c3", "c1"), merged_id = "JUNC_MERGED_1", max_reads = c(3, 0, 5, 2))
  s = sc_junction_support_from_counts(counts, m)
  expect_equal(nrow(s), 2L)
  expect_equal(s[cell == "c1"]$reads, 3)
  expect_false("c2" %in% s$cell)
})
