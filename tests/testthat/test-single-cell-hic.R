## Synthetic matrices only.

.read_hic_bin <- function(path, sec) {
    con <- file(path, "rb"); on.exit(close(con))
    raw_all <- readBin(con, "raw", n = file.size(path))
    take <- function(byte, n, what) {
        if (n == 0) return(if (what == "int") integer(0) else numeric(0))
        r <- raw_all[byte + seq_len(4 * n)]
        readBin(r, what, n = n, size = 4L, endian = "little")
    }
    off <- take(sec$offsets$byte, sec$offsets$n, "int")
    w <- take(sec$weights$byte, sec$weights$n, "numeric")
    np <- sec$pixels$n
    px <- if (np) matrix(raw_all[sec$pixels$byte + seq_len(8 * np)], nrow = 8L) else matrix(raw(), nrow = 8L)
    list(offsets = off, weights = w,
         bin2 = if (np) readBin(as.vector(px[1:4, ]), "int", n = np, size = 4L, endian = "little") else integer(0),
         count = if (np) readBin(as.vector(px[5:8, ]), "numeric", n = np, size = 4L, endian = "little") else numeric(0))
}

.write_synthetic_cool <- function(path, group, chroms, lengths, res, pix) {
    nb_c <- ceiling(lengths / res)
    coff <- c(0, cumsum(nb_c))
    nb <- sum(nb_c)
    rhdf5::h5createGroup(path, group)
    for (g in c("chroms", "bins", "pixels", "indexes")) rhdf5::h5createGroup(path, paste0(group, "/", g))
    rhdf5::h5write(chroms, path, paste0(group, "/chroms/name"))
    rhdf5::h5write(as.integer(lengths), path, paste0(group, "/chroms/length"))
    bc <- rep(seq_along(chroms) - 1L, nb_c)
    st <- unlist(lapply(nb_c, function(n) (seq_len(n) - 1L) * res))
    rhdf5::h5write(as.integer(bc), path, paste0(group, "/bins/chrom"))
    rhdf5::h5write(as.integer(st), path, paste0(group, "/bins/start"))
    pix <- pix[order(pix$b1, pix$b2), ]
    rhdf5::h5write(as.integer(pix$b1), path, paste0(group, "/pixels/bin1_id"))
    rhdf5::h5write(as.integer(pix$b2), path, paste0(group, "/pixels/bin2_id"))
    rhdf5::h5write(as.integer(pix$v), path, paste0(group, "/pixels/count"))
    rhdf5::h5write(as.integer(coff), path, paste0(group, "/indexes/chrom_offset"))
    rhdf5::h5write(as.integer(c(0, cumsum(tabulate(pix$b1 + 1, nb)))), path, paste0(group, "/indexes/bin1_offset"))
    nb
}

test_that("sc_hic_write_resolution round-trips pixels, offsets and weights", {
    d <- tempfile(); dir.create(d)
    b1 <- c(0, 0, 1, 3, 3); b2 <- c(0, 2, 1, 3, 4); ct <- c(5, 1, 2, 7, 3)
    wr <- sc_hic_write_resolution(file.path(d, "r10.bin"), b1, b2, ct, 5, weights = c(1, NA, 0.5, 2, 1))
    expect_equal(wr$nnz, 5)
    expect_equal(wr$bytes, file.size(file.path(d, "r10.bin")))
    r <- .read_hic_bin(file.path(d, "r10.bin"), wr$sections)
    expect_equal(r$offsets, c(0L, 2L, 3L, 3L, 5L, 5L))
    expect_equal(r$bin2, as.integer(b2))
    expect_equal(r$count, ct)
    expect_true(is.nan(r$weights[2]))
    expect_equal(r$weights[c(1, 3, 4, 5)], c(1, 0.5, 2, 1))
    expect_error(sc_hic_write_resolution(file.path(d, "bad.bin"), c(1, 0), c(1, 0), c(1, 1), 2), "sorted")
})

test_that("sc_hic_ice equalizes the marginals of a synthetic matrix", {
    skip_if_not_installed("Matrix")
    set.seed(1)
    n <- 60
    bias <- runif(n, 0.5, 2)
    pr <- expand.grid(b1 = 0:(n - 1), b2 = 0:(n - 1)); pr <- pr[pr$b1 <= pr$b2, ]
    pr$v <- round(1000 * bias[pr$b1 + 1] * bias[pr$b2 + 1] / (1 + abs(pr$b2 - pr$b1)))
    pr <- pr[pr$v > 0, ]
    ice <- sc_hic_ice(pr$b1, pr$b2, pr$v, n)
    expect_true(ice$converged)
    keep <- (pr$b2 - pr$b1) >= 2
    bal <- pr$v[keep] * ice$weights[pr$b1[keep] + 1] * ice$weights[pr$b2[keep] + 1]
    marg <- as.numeric(rowsum(c(bal, bal), c(pr$b1[keep], pr$b2[keep]))[, 1])
    expect_lt(stats::sd(marg) / mean(marg), 0.01)
    ## weights undo the bias up to a constant
    expect_gt(stats::cor(1 / ice$weights, bias), 0.9)   # up to edge effects of the distance decay
})

test_that("sc_export_hic_tiles exports a synthetic mcool, drops other chromosomes and renumbers bins", {
    skip_if_not_installed("rhdf5")
    d <- tempfile(); dir.create(d)
    f <- file.path(d, "toy.mcool")
    rhdf5::h5createFile(f)
    rhdf5::h5createGroup(f, "resolutions")
    ## chr1 (3 bins), chrUn_x (1 bin, dropped), chr2 (2 bins) at res 100
    pix <- data.frame(b1 = c(0, 0, 1, 2, 3, 4, 0, 5), b2 = c(0, 4, 2, 5, 3, 5, 3, 5), v = c(9, 4, 3, 2, 8, 1, 6, 5))
    .write_synthetic_cool(f, "resolutions/100", c("chr1", "chrUn_x", "chr2"), c(300, 50, 200), 100, pix)
    pix2 <- data.frame(b1 = c(0, 0, 1), b2 = c(0, 2, 2), v = c(20, 3, 7))
    .write_synthetic_cool(f, "resolutions/200", c("chr1", "chrUn_x", "chr2"), c(300, 50, 200), 200, pix2)
    idx <- sc_export_hic_tiles(f, file.path(d, "hic"), unit = "U1", chroms = c("1", "2"), balance = FALSE, data_label = "synthetic")
    expect_equal(vapply(idx$resolutions, `[[`, 0, "res"), c(100, 200))
    r100 <- idx$resolutions[[1]]
    expect_equal(r100$nbins, 5)
    expect_equal(vapply(r100$chroms, `[[`, "", "name"), c("1", "2"))
    expect_equal(vapply(r100$chroms, `[[`, 0, "bin_offset"), c(0, 3))
    rd <- .read_hic_bin(file.path(d, "hic", "r100.bin"), r100$sections)
    ## old bins 0,1,2 -> 0,1,2; old 4,5 -> 3,4; old 3 (chrUn_x) dropped with its pixels
    expect_equal(rd$offsets, c(0L, 2L, 3L, 4L, 5L, 6L))
    expect_equal(rd$bin2, c(0L, 3L, 2L, 4L, 4L, 4L))
    expect_equal(rd$count, c(9, 4, 3, 2, 1, 5))
    expect_true(all(is.nan(rd$weights)))
    j <- jsonlite::fromJSON(file.path(d, "hic", "index.json"), simplifyVector = FALSE)
    expect_equal(j$format, "gos-hic-tiles/1")
    expect_equal(j$data_label, "synthetic")
    expect_equal(j$resolutions[[1]]$sections$pixels$byte, 4 * 6 + 4 * 5)
    sc_export_hic_tiles(f, file.path(d, "hic"), unit = "U1", resolutions = 200, chroms = c("1", "2"), balance = FALSE)
    expect_length(list.files(file.path(d, "hic"), "bak_"), 1)
})
