## Hi-C contact maps for the gOS track browser ("gos-hic-tiles/1").
##
## The browser cannot read .mcool (HDF5), so each unit's .mcool is re-written as one small JSON index plus one
## flat little-endian binary per resolution that the browser reads with HTTP range requests:
##
##   <out_dir>/index.json      format, unit, genome, chroms, resolutions[{res, file, nbins, nnz, sum, chroms
##                             (name, length, bin_offset), sections (byte offsets), balance}]
##   <out_dir>/r<res>.bin      [bin1 offsets u32 x (nbins + 1)] [balancing weights f32 x nbins, NaN = filtered]
##                             [pixels: (bin2 u32, count f32) x nnz]
##
## Pixels are the cooler upper triangle (bin1 <= bin2), sorted by bin1 then bin2, exactly as in the .mcool, so the
## rows of a window are one contiguous byte range: offsets[a] .. offsets[b + 1]. Contacts between two windows are
## found in the rows of the earlier one. Bins are renumbered over the kept chromosomes only (default 1-22, X, Y,
## in file order), chromosome names lose the "chr" prefix. Balanced value = count * w[bin1] * w[bin2].
##
## Balancing: cooler's iterative correction (ICE) defaults (ignore 2 diagonals, min 10 non-zero pixels per bin,
## MAD-max 5 on log marginals, tol 1e-5, at most 200 iterations), computed here because the source .mcool files
## carry no weights. Amplified bins (ecDNA) are kept unless the MAD filter drops them; balancing flattens copy
## number by design, so raw counts are the default view in gOS.

#' @name sc_hic_ice
#' @title sc_hic_ice
#' @description Iterative correction (ICE) weights for an upper-triangle sparse contact matrix, cooler defaults.
#' @param bin1 0-based row bin ids (bin1 <= bin2)
#' @param bin2 0-based column bin ids
#' @param count contact counts
#' @param nbins number of bins
#' @param ignore_diags diagonals ignored (0 = main diagonal); cooler default 2
#' @param min_nnz bins with fewer non-zero pixels are filtered
#' @param mad_max bins whose log marginal is more than mad_max MADs below the median are filtered
#' @param max_iter maximum iterations
#' @param tol stop when the variance of the marginals is below tol
#' @return list(weights = numeric nbins (NA = filtered), converged, iterations, n_filtered, scale)
#' @export
sc_hic_ice <- function(bin1, bin2, count, nbins, ignore_diags = 2L, min_nnz = 10L, mad_max = 5, max_iter = 200L, tol = 1e-5) {
    if (!requireNamespace("Matrix", quietly = TRUE)) stop("sc_hic_ice needs the Matrix package")
    keep <- (bin2 - bin1) >= ignore_diags & count > 0
    b1 <- as.integer(bin1[keep]) + 1L
    b2 <- as.integer(bin2[keep]) + 1L
    v <- as.numeric(count[keep])
    nnz <- tabulate(c(b1, b2[b2 != b1]), nbins)
    bad <- nnz < min_nnz
    A <- Matrix::sparseMatrix(i = b1, j = b2, x = v, dims = c(nbins, nbins), symmetric = TRUE)
    marg0 <- as.numeric(A %*% rep(1, nbins))
    lm <- log(marg0[marg0 > 0 & !bad])
    if (length(lm)) {
        med <- stats::median(lm)
        dev <- stats::median(abs(lm - med))
        bad <- bad | marg0 < exp(med - mad_max * dev)
    }
    bad <- bad | marg0 == 0
    w <- ifelse(bad, 0, 1)
    converged <- FALSE
    it <- 0L
    for (it in seq_len(max_iter)) {
        marg <- w * as.numeric(A %*% w)
        mm <- mean(marg[!bad])
        if (!is.finite(mm) || mm == 0) break
        marg <- marg / mm
        marg[bad] <- 1
        w <- w / marg
        vr <- stats::var(marg[!bad])
        if (is.na(vr) || vr < tol) { converged <- TRUE; break }
    }
    ## rescale so balanced marginals average 1 (cooler rescale_marginals)
    marg <- w * as.numeric(A %*% w)
    scale <- mean(marg[!bad])
    if (is.finite(scale) && scale > 0) w <- w / sqrt(scale)
    w[bad] <- NA_real_
    list(weights = w, converged = converged, iterations = it, n_filtered = sum(bad), scale = scale)
}

## little-endian 4-byte writers
.hic_u32 <- function(x) writeBin(as.integer(x), raw(), size = 4L, endian = "little")
.hic_f32 <- function(x) writeBin(as.numeric(x), raw(), size = 4L, endian = "little")

#' @name sc_hic_write_resolution
#' @title sc_hic_write_resolution
#' @description Writes one resolution of a gos-hic-tiles/1 export from in-memory pixels (used by
#' sc_export_hic_tiles; also the entry point for tests with synthetic matrices).
#' @param path output .bin path
#' @param bin1,bin2,count 0-based upper-triangle pixels over the kept bins, sorted by bin1 then bin2
#' @param nbins number of kept bins
#' @param weights balancing weights (length nbins, NA = filtered) or NULL
#' @return list(nnz, sum, bytes, sections)
#' @export
sc_hic_write_resolution <- function(path, bin1, bin2, count, nbins, weights = NULL) {
    stopifnot(length(bin1) == length(bin2), length(bin1) == length(count), all(bin1 <= bin2))
    if (length(bin1) > 1 && is.unsorted(as.numeric(bin1) * nbins + as.numeric(bin2))) stop("pixels must be sorted by bin1, bin2")
    offsets <- c(0L, cumsum(tabulate(as.integer(bin1) + 1L, nbins)))
    w <- if (is.null(weights)) rep(NaN, nbins) else ifelse(is.na(weights), NaN, weights)
    px <- if (length(bin2)) as.vector(rbind(matrix(.hic_u32(bin2), nrow = 4L), matrix(.hic_f32(count), nrow = 4L))) else raw()
    con <- file(path, "wb")
    on.exit(close(con))
    writeBin(.hic_u32(offsets), con)
    writeBin(.hic_f32(w), con)
    writeBin(px, con)
    n_off <- 4 * (nbins + 1)
    list(nnz = length(bin1), sum = sum(as.numeric(count)), bytes = n_off + 4 * nbins + 8 * length(bin1),
         sections = list(offsets = list(byte = 0, n = nbins + 1, dtype = "u32"),
                         weights = list(byte = n_off, n = nbins, dtype = "f32"),
                         pixels = list(byte = n_off + 4 * nbins, n = length(bin1), dtype = "u32 bin2 + f32 count", stride = 8)))
}

#' @name sc_export_hic_tiles
#' @title sc_export_hic_tiles
#' @description
#'
#' Exports one unit's .mcool (cooler multi-resolution file) as gos-hic-tiles/1 for the gOS track browser `hic`
#' track: index.json plus one range-readable binary per resolution (layout in the file header comment of
#' R/single-cell-hic.R). Raw counts are copied as stored; ICE weights are computed per resolution (cooler
#' defaults) unless balance = FALSE. An existing index.json is backed up as index.json.bak_<yymmdd>.
#'
#' @param mcool path to the .mcool (or a single-resolution .cool when resolutions has one entry)
#' @param out_dir output folder (e.g. <data>/<unit>/hic)
#' @param unit unit id written into the index
#' @param resolutions bin sizes to export (default: all in the file)
#' @param chroms chromosomes to keep, without "chr" (default 1-22, X, Y)
#' @param balance compute ICE weights
#' @param genome genome label for the index
#' @param data_label provenance label shown by gOS (e.g. "public data ...")
#' @param max_balance_bins skip balancing (weights NaN, balance$method "none") above this many bins
#' @return the index (list), invisibly
#' @export
sc_export_hic_tiles <- function(mcool, out_dir, unit, resolutions = NULL, chroms = c(1:22, "X", "Y"), balance = TRUE,
                                genome = "hg38", data_label = NULL, max_balance_bins = 4e6) {
    if (!requireNamespace("rhdf5", quietly = TRUE)) stop("sc_export_hic_tiles needs rhdf5")
    groups <- tryCatch(rhdf5::h5ls(mcool, recursive = 2), error = function(e) NULL)
    is_multi <- !is.null(groups) && any(groups$group == "/resolutions")
    avail <- if (is_multi) sort(as.integer(groups$name[groups$group == "/resolutions"])) else NA_integer_
    if (is.null(resolutions)) resolutions <- avail
    if (is_multi && !all(resolutions %in% avail)) stop("resolutions not in ", mcool, ": ", paste(setdiff(resolutions, avail), collapse = ","))
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    chroms <- sub("^chr", "", as.character(chroms))
    res_out <- list()
    for (res in resolutions) {
        g <- if (is_multi) paste0("resolutions/", res, "/") else ""
        rd <- function(p, ...) as.vector(rhdf5::h5read(mcool, paste0(g, p), ...))
        cn <- sub("^chr", "", rd("chroms/name"))
        cl <- as.numeric(rd("chroms/length"))
        coff <- as.numeric(rd("indexes/chrom_offset"))           # nchroms + 1, 0-based bin offsets
        b1off <- as.numeric(rd("indexes/bin1_offset"))           # nbins + 1
        binsize <- if (is_multi) as.integer(res) else as.integer(rhdf5::h5readAttributes(mcool, "/")$`bin-size`)
        kc <- which(cn %in% chroms)                              # file order, so the bin remap is monotonic
        if (!length(kc)) stop("none of the chromosomes are in ", mcool)
        nb_c <- coff[kc + 1] - coff[kc]
        new_off <- c(0, cumsum(nb_c))
        nb <- sum(nb_c)
        map <- rep(NA_real_, coff[length(coff)])                # old bin -> new bin (0-based)
        for (k in seq_along(kc)) map[(coff[kc[k]] + 1):coff[kc[k] + 1]] <- new_off[k] + seq_len(nb_c[k]) - 1
        ## read only the pixel rows whose bin1 is on a kept chromosome
        rows <- range(unlist(lapply(kc, function(i) c(coff[i], coff[i + 1]))))
        p0 <- b1off[rows[1] + 1]; p1 <- b1off[rows[2] + 1]
        n <- p1 - p0
        if (n > 0) {
            b1 <- rd("pixels/bin1_id", start = p0 + 1, count = n)
            b2 <- rd("pixels/bin2_id", start = p0 + 1, count = n)
            ct <- rd("pixels/count", start = p0 + 1, count = n)
            nb1 <- map[b1 + 1]; nb2 <- map[b2 + 1]
            ok <- !is.na(nb1) & !is.na(nb2)
            nb1 <- nb1[ok]; nb2 <- nb2[ok]; ct <- ct[ok]
            rm(b1, b2, ok)
        } else nb1 <- nb2 <- ct <- numeric(0)
        bal <- list(method = "none")
        w <- NULL
        if (balance && nb <= max_balance_bins && length(ct)) {
            ice <- sc_hic_ice(nb1, nb2, ct, nb)
            w <- ice$weights
            bal <- list(method = "ICE", ignore_diags = 2, min_nnz = 10, mad_max = 5, converged = ice$converged,
                        iterations = ice$iterations, n_filtered = ice$n_filtered)
        }
        f <- paste0("r", binsize, ".bin")
        wr <- sc_hic_write_resolution(file.path(out_dir, f), nb1, nb2, ct, nb, w)
        res_out[[length(res_out) + 1]] <- list(
            res = binsize, file = f, nbins = nb, nnz = wr$nnz, sum = wr$sum, bytes = wr$bytes,
            chroms = lapply(seq_along(kc), function(k) list(name = cn[kc[k]], length = cl[kc[k]], bin_offset = new_off[k], nbins = nb_c[k])),
            sections = wr$sections, balance = bal)
        rm(nb1, nb2, ct, w); gc()
    }
    idx <- list(format = "gos-hic-tiles/1", unit = unit, genome = genome, source = basename(mcool),
                storage = "symmetric-upper", created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
                data_label = data_label, resolutions = res_out)
    out <- file.path(out_dir, "index.json")
    if (file.exists(out)) file.copy(out, paste0(out, ".bak_", format(Sys.Date(), "%y%m%d")), overwrite = FALSE)
    jsonlite::write_json(idx, out, auto_unbox = TRUE, digits = NA, null = "null", na = "null", pretty = FALSE)
    invisible(idx)
}
