## Per-cell allelic copy number: reconcile the allelic JaBbA graph with the
## cell's own germline het-SNP reads (2026-10-20).
##
## Why: the allelic graph (LPPhasedBalance_blocks) splits the *total* JaBbA CN
## between the two haplotypes. The total-CN marginal carries the same weight as
## the het evidence (both = n het SNPs), so wherever the total CN is under-called
## (GC-rich, PTA-depleted chr16/17/19/22) a node with total CN 1 is forced to
## 1:0 even when most het sites show both alleles. These functions test LOH
## against the reads directly, with an explicit PTA allelic-dropout model
## fitted per cell on its 1:1 segments, and correct calls the reads reject.
##
## Site model (folded count k = min(ref, alt), depth d >= 2), per het site:
##   with prob r      allelic dropout: one allele only (error e)
##   with prob w      artefact site: behaves balanced whatever the CN (paralog,
##                    false germline het); w is a small fixed floor
##   otherwise        both alleles amplified, minor copy fraction f, PTA
##                    amplification bias as beta-binomial overdispersion rho
## LOH is f = 0. r and rho are fitted per cell (f = 1/2 on its 1:1 autosomal
## segments); the LOH-vs-het log-likelihood ratio of a segment is scaled by
## (independent 50 kb bins / sites) because neighbouring sites share PTA
## amplicons.

## log beta-binomial density, vectorised over k, d (p, rho scalars); rho = 0 is binomial
.sc_lbb <- function(k, d, p, rho) {
    if (rho <= 0) return(stats::dbinom(k, d, p, log = TRUE))
    a <- p * (1 - rho) / rho
    b <- (1 - p) * (1 - rho) / rho
    lchoose(d, k) + lbeta(k + a, d - k + b) - lbeta(a, b)
}

#' @name sc_ado_loglik
#' @title sc_ado_loglik
#' @description
#'
#' Log-likelihood of tabulated folded het-SNP counts under the single-cell
#' allelic-dropout model (see the header of single-cell-allelic.R).
#'
#' @param k folded minor read count per site (or per table row)
#' @param d depth per site (>= 2)
#' @param n number of sites with this (k, d) (default 1 each)
#' @param f minor copy fraction (0 = LOH, 0.5 = balanced)
#' @param r allelic-dropout probability
#' @param rho beta-binomial overdispersion of amplified sites
#' @param w artefact-site fraction (balanced whatever the CN)
#' @param err per-read error rate
#' @return summed log-likelihood
#' @export
#' @author Stanley Clarke
sc_ado_loglik <- function(k, d, n = rep(1, length(k)), f, r = 0.1, rho = 0.05, w = 0.02, err = 0.01) {
    if (!length(k)) return(0)
    sum(n * .sc_ado_rows(k, d, f, r, rho, w, err))
}

## per-row log-likelihood (unsummed)
.sc_ado_rows <- function(k, d, f, r, rho, w, err) {
    p <- f * (1 - 2 * err) + err
    ldrop <- log(0.5) + matrixStats_logsumexp2(stats::dbinom(k, d, err, log = TRUE), stats::dbinom(k, d, 1 - err, log = TRUE))
    lamp <- log(0.5) + matrixStats_logsumexp2(.sc_lbb(k, d, p, rho), .sc_lbb(k, d, 1 - p, rho))
    lart <- .sc_lbb(k, d, 0.5, rho)
    q_amp <- max(1 - r - w, 1e-9)
    matrixStats_logsumexp3(log(max(r, 1e-9)) + ldrop, log(q_amp) + lamp, log(max(w, 1e-9)) + lart)
}

matrixStats_logsumexp2 <- function(a, b) {
    m <- pmax(a, b)
    m + log(exp(a - m) + exp(b - m))
}
matrixStats_logsumexp3 <- function(a, b, c) {
    m <- pmax(a, b, c)
    m + log(exp(a - m) + exp(b - m) + exp(c - m))
}

#' @name sc_fit_ado
#' @title sc_fit_ado
#' @description
#'
#' Fit the per-cell dropout rate r and overdispersion rho by grid search on
#' sites from segments believed balanced (f = 1/2).
#'
#' @param kd data.table with k, d, n (tabulated folded counts)
#' @param r_grid,rho_grid search grids
#' @param w,err fixed artefact fraction and error rate
#' @return list(r, rho, n_sites)
#' @export
#' @author Stanley Clarke
sc_fit_ado <- function(kd, r_grid = seq(0, 0.9, by = 0.025), rho_grid = c(0, 0.01, 0.03, 0.06, 0.1, 0.2, 0.3),
                       w = 0.02, err = 0.01) {
    if (!nrow(kd)) return(list(r = NA_real_, rho = NA_real_, n_sites = 0))
    best <- c(-Inf, NA, NA)
    for (rho in rho_grid) for (r in r_grid) {
        ll <- sc_ado_loglik(kd$k, kd$d, kd$n, f = 0.5, r = r, rho = rho, w = w, err = err)
        if (ll > best[1]) best <- c(ll, r, rho)
    }
    list(r = best[2], rho = best[3], n_sites = sum(kd$n))
}

#' @name sc_read_het_counts
#' @title sc_read_het_counts
#' @description
#'
#' Per-cell germline het-SNP read counts as data.table(chr, pos, ref, alt), chr
#' without the "chr" prefix. Reads either the genotyping-run counts written by
#' the gOS het-SNP exporter (rds with mutation = chr_pos_ref_alt, ref_count,
#' alt_count) or an Eagle hets_with_phase table (seqnames, start, ref.count.t,
#' alt.count.t; sites with a normal allele fraction outside 0.2-0.8 are dropped).
#'
#' @param path rds or txt/tsv path
#' @return data.table(chr, pos, ref, alt)
#' @export
#' @author Stanley Clarke
sc_read_het_counts <- function(path) {
    if (grepl("\\.rds$", path, ignore.case = TRUE)) {
        x <- data.table::as.data.table(readRDS(path))
        if (!"mutation" %in% names(x)) stop("het count rds needs mutation, ref_count, alt_count")
        sp <- data.table::tstrsplit(x$mutation, "_", fixed = TRUE, keep = 1:2)
        out <- data.table::data.table(chr = sub("^chr", "", sp[[1]]), pos = as.integer(sp[[2]]),
                                      ref = as.integer(x$ref_count), alt = as.integer(x$alt_count))
    } else {
        x <- data.table::fread(path)
        if ("alt.count.n" %in% names(x)) {
            fn <- x$alt.count.n / (x$alt.count.n + x$ref.count.n)
            x <- x[is.finite(fn) & fn > 0.2 & fn < 0.8]
        }
        out <- data.table::data.table(chr = sub("^chr", "", as.character(x$seqnames)), pos = as.integer(x$start),
                                      ref = as.integer(x$ref.count.t), alt = as.integer(x$alt.count.t))
    }
    out[!is.na(ref) & !is.na(alt)]
}

## hg38 pseudo-autosomal regions (sites there are diploid in males)
.SC_PAR <- data.table::data.table(chr = "X", start = c(10001L, 155701383L), end = c(2781479L, 156030895L))

#' @name sc_allelic_baf_fix
#' @title sc_allelic_baf_fix
#' @description
#'
#' Test every allelic segment of a cell against its het-SNP reads and correct
#' the calls the reads reject:
#'  * called LOH (minor 0) but the reads support both alleles (scaled LLR >
#'    llr_min): minor is set to the best-supported value in 1..floor(total/2);
#'    a total of 1 becomes 1:1 (total raised to 2, status
#'    loh_rejected_total_raised: the total CN graph is under-called there);
#'  * called heterozygous but the reads support LOH (scaled LLR < -llr_min):
#'    copy-neutral LOH, major = total, minor = 0 (status loh_by_baf; only when
#'    fix_het = TRUE);
#'  * male chrX / chrY outside the PARs: hemizygous, minor = 0, never tested;
#'  * fewer than min_sites informative sites: unchanged, status few_sites;
#'  * part LOH, part het (chunk check) or weak minor fraction: unchanged,
#'    status mixed_breakpoint / loh_partial.
#'
#' @param segs data.table(seg_id, chr, start, end, major, minor); major/minor are CN magnitudes
#' @param hets data.table(chr, pos, ref, alt) from sc_read_het_counts()
#' @param sex "female", "male", "unknown" or "auto" (female when non-PAR chrX het
#'   sites with reads are at least x_ratio of the autosomal density per Mb,
#'   unknown when there are none); with unknown, chrX is tested like an
#'   autosome only when it has sites (else few_sites)
#' @param min_sites minimum sites with depth >= 2 for a segment to be tested
#' @param llr_min threshold on the scaled log-likelihood ratio
#' @param bin_width bin size used to scale the LLR for correlated sites
#' @param fix_het also apply the het -> LOH direction
#' @param w,err artefact fraction and error rate
#' @param x_ratio sex call threshold for sex = "auto"
#' @param chunk_sites sites per chunk for the within-segment heterogeneity check;
#'   a segment with a clearly LOH chunk and a clearly het chunk is left unchanged
#'   (status mixed_breakpoint: the segment spans an allelic breakpoint)
#' @param max_total_het_test het calls with a larger total CN are never turned
#'   into LOH (status high_cn_untested): one minor copy among many is not
#'   separable from LOH with a few reads per site
#' @param dcap depth cap per site
#' @param mixed_frac minimum fraction of chunks of each sign for mixed_breakpoint
#' @param min_f_reject an LOH call is only rejected when the best minor fraction
#'   is at least this (else status loh_partial, unchanged)
#' @param min_sites_reject an LOH call on a segment with fewer than min_sites
#'   (but at least min_sites_reject) sites is still tested, in the LOH -> het
#'   direction only (status loh_rejected / loh_rejected_total_raised /
#'   loh_partial; otherwise few_sites). A fine segmentation splits arms into
#'   segments with a few dozen sites each, which min_sites alone leaves untested.
#'   Default min_sites (off).
#' @return list(segs = segs with new_major, new_minor, status, n_sites, llr, f_hat; params = list(r, rho, sex, ...))
#' @export
#' @author Stanley Clarke
sc_allelic_baf_fix <- function(segs, hets, sex = "auto", min_sites = 30L, llr_min = 10, bin_width = 5e4,
                               fix_het = TRUE, w = 0.02, err = 0.01, x_ratio = 0.2,
                               chunk_sites = 300L, mixed_frac = 0.2, min_f_reject = 0.2, dcap = 40L,
                               max_total_het_test = 4L, min_sites_reject = min_sites) {
    segs <- data.table::copy(data.table::as.data.table(segs))
    segs[, chr := sub("^chr", "", as.character(chr))]
    hets <- data.table::as.data.table(hets)[(ref + alt) >= 2]
    hets[, chr := sub("^chr", "", as.character(chr))]
    ## PAR flag and sex
    hets[, par := chr == "X" & ((pos >= .SC_PAR$start[1] & pos <= .SC_PAR$end[1]) | (pos >= .SC_PAR$start[2] & pos <= .SC_PAR$end[2]))]
    if (identical(sex, "auto")) {
        auto_len <- 2.875e3  # Mb, hg38 autosomes
        x_len <- 153.4       # Mb, chrX outside the PARs
        dens_a <- hets[chr %in% as.character(1:22), .N] / auto_len
        dens_x <- hets[chr == "X" & !par, .N] / x_len
        sex <- if (hets[chr == "X" & !par, .N] == 0) "unknown" else if (dens_a > 0 && dens_x >= x_ratio * dens_a) "female" else "male"
    }
    if (!sex %in% c("female", "male", "unknown")) stop("sex must be female, male, unknown or auto")
    ## site -> segment
    hets[, `:=`(start = pos, end = pos)]
    data.table::setkey(segs, chr, start, end)
    ov <- data.table::foverlaps(hets, segs[, .(seg_id, chr, start, end)], by.x = c("chr", "start", "end"), type = "within", nomatch = NULL)
    ov[, bin := pos %/% bin_width]
    ## chunks of ~chunk_sites consecutive sites (heterogeneity check: LOH over part of a segment only)
    data.table::setorder(ov, seg_id, pos)
    ov[, chunk := {
        nck <- max(1L, .N %/% chunk_sites)
        pmin(nck, (seq_len(.N) - 1L) %/% chunk_sites + 1L)
    }, by = seg_id]
    ## folded counts, depth capped at dcap (deep sites add little once both alleles are seen)
    ov[, d := ref + alt]
    ov[, k := pmin(ref, alt)]
    ov[d > dcap, `:=`(k = as.integer(round(k * dcap / d)), d = dcap)]
    kdc <- ov[, .(n = .N), by = .(seg_id, chunk, k, d)]
    kd <- kdc[, .(n = sum(n)), by = .(seg_id, k, d)]
    cbin <- ov[, .(s = min(1, data.table::uniqueN(bin) / .N)), by = .(seg_id, chunk)]
    nsite <- ov[, .(n_sites = .N, n_bins = data.table::uniqueN(bin), n_both = sum(ref > 0 & alt > 0)), by = seg_id]
    segs <- merge(segs, nsite, by = "seg_id", all.x = TRUE)
    segs[is.na(n_sites), `:=`(n_sites = 0L, n_bins = 0L, n_both = 0L)]
    ## male X/Y are hemizygous, except segments lying entirely inside a PAR (X+Y diploid)
    segs[, hemi := sex == "male" & chr %in% c("X", "Y") &
             !((start >= .SC_PAR$start[1] & end <= .SC_PAR$end[1]) | (start >= .SC_PAR$start[2] & end <= .SC_PAR$end[2]))]
    ## calibration: autosomal 1:1 segments (fall back to all autosomal non-LOH)
    ref_ids <- segs[chr %in% as.character(1:22) & major == 1 & minor == 1 & n_sites >= min_sites, seg_id]
    if (sum(segs[seg_id %in% ref_ids, n_sites]) < 500)
        ref_ids <- segs[chr %in% as.character(1:22) & minor >= 1 & n_sites >= min_sites, seg_id]
    fit <- sc_fit_ado(kd[seg_id %in% ref_ids, .(n = sum(n)), by = .(k, d)], w = w, err = err)
    ## one refinement: drop calibration segments that themselves look LOH (copy-neutral LOH called 1:1) and refit
    if (!is.na(fit$r) && length(ref_ids)) {
        lr <- vapply(ref_ids, function(id) {
            z <- kd[seg_id == id]
            sc_ado_loglik(z$k, z$d, z$n, f = 0.5, r = fit$r, rho = fit$rho, w = w, err = err) -
                sc_ado_loglik(z$k, z$d, z$n, f = 0, r = fit$r, rho = fit$rho, w = w, err = err)
        }, 0)
        if (any(lr < 0)) {
            keep_ids <- ref_ids[lr >= 0]
            if (sum(segs[seg_id %in% keep_ids, n_sites]) >= 500)
                fit <- sc_fit_ado(kd[seg_id %in% keep_ids, .(n = sum(n)), by = .(k, d)], w = w, err = err)
        }
    }
    segs[, `:=`(new_major = major, new_minor = minor, status = "unchanged", llr = NA_real_, f_hat = NA_real_, chunk_llr_range = NA_character_)]
    segs[(hemi), `:=`(new_major = major + minor, new_minor = 0L, status = "hemizygous_x")]
    if (is.na(fit$r)) {
        segs[!(hemi), status := "no_calibration"]
        return(list(segs = segs[], params = list(r = NA, rho = NA, sex = sex, n_ref_sites = 0)))
    }
    fgrid <- c(1/2, 2/5, 1/3, 1/4, 1/5, 1/6)
    ## log-likelihood of every (k, d) row under each minor fraction, then summed per chunk / segment
    fall <- c(0, fgrid)
    ukd <- unique(kdc[, .(k, d)])
    llm <- vapply(fall, function(f) .sc_ado_rows(ukd$k, ukd$d, f, fit$r, fit$rho, w, err), numeric(nrow(ukd)))
    if (is.null(dim(llm))) llm <- matrix(llm, nrow = 1)
    ukd[, row := .I]
    kdc <- merge(kdc, ukd, by = c("k", "d"))
    fcols <- paste0("f", seq_along(fall))
    kdc[, (fcols) := lapply(seq_along(fall), function(j) n * llm[row, j])]
    seg_ll <- kdc[, lapply(.SD, sum), by = seg_id, .SDcols = fcols]
    chunk_ll <- merge(kdc[, lapply(.SD, sum), by = .(seg_id, chunk), .SDcols = fcols], cbin, by = c("seg_id", "chunk"))
    j0 <- which(fall == 0); jg <- match(fgrid, fall)
    best_het <- function(m) apply(m[, jg, drop = FALSE], 1, max)
    chunk_ll[, l := s * (best_het(as.matrix(.SD)) - .SD[[j0]]), .SDcols = fcols]
    cmix <- chunk_ll[, .(nck = .N, lmin = min(l), lmax = max(l),
                         fpos = mean(l > llr_min / 4), fneg = mean(l < -llr_min / 4)), by = seg_id]
    seg_mat <- as.matrix(seg_ll[, fcols, with = FALSE])
    seg_ll[, `:=`(ll0 = seg_mat[, j0], llhet = best_het(seg_mat), fbest = fgrid[apply(seg_mat[, jg, drop = FALSE], 1, which.max)])]
    rowsel <- match(segs$seg_id, seg_ll$seg_id)
    cm <- cmix[match(segs$seg_id, cmix$seg_id)]
    for (i in which(!segs$hemi)) {
        tot <- segs$major[i] + segs$minor[i]
        few <- segs$n_sites[i] < min_sites
        if (is.na(rowsel[i]) || (few && !(segs$minor[i] == 0 && tot > 0 && segs$n_sites[i] >= min_sites_reject))) {
            data.table::set(segs, i, "status", "few_sites"); next
        }
        if (tot == 0) { data.table::set(segs, i, "status", "cn0"); next }
        r_i <- rowsel[i]
        s <- min(1, segs$n_bins[i] / segs$n_sites[i])
        llr <- s * (seg_ll$llhet[r_i] - seg_ll$ll0[r_i])
        fhat <- if (llr > 0) seg_ll$fbest[r_i] else 0
        data.table::set(segs, i, c("llr", "f_hat"), list(llr, fhat))
        ## chunk LLRs: a clearly LOH part and a clearly het part -> the segment spans an allelic breakpoint
        ## (artefact clusters such as segmental duplications make isolated het-looking chunks inside true LOH,
        ##  so both signs must cover at least mixed_frac of the chunks)
        mixed <- cm$nck[i] > 1 && cm$fpos[i] >= mixed_frac && cm$fneg[i] >= mixed_frac
        if (cm$nck[i] > 1) data.table::set(segs, i, "chunk_llr_range", sprintf("%.1f..%.1f (het %.2f, loh %.2f)", cm$lmin[i], cm$lmax[i], cm$fpos[i], cm$fneg[i]))
        mn <- segs$minor[i]
        if (few && !(llr > llr_min)) {
            ## below min_sites only the LOH -> het direction is applied
            data.table::set(segs, i, "status", "few_sites")
        } else if (mixed) {
            data.table::set(segs, i, "status", "mixed_breakpoint")
        } else if (mn == 0 && llr > llr_min && fhat < min_f_reject) {
            data.table::set(segs, i, "status", "loh_partial")
        } else if (mn == 0 && llr > llr_min) {
            if (tot >= 2) {
                m <- seq_len(floor(tot / 2))
                z <- kd[seg_id == segs$seg_id[i]]
                llm <- vapply(m / tot, function(f) sc_ado_loglik(z$k, z$d, z$n, f = f, r = fit$r, rho = fit$rho, w = w, err = err), 0)
                mbest <- m[which.max(llm)]
                data.table::set(segs, i, c("new_major", "new_minor", "status"), list(as.integer(tot - mbest), as.integer(mbest), "loh_rejected"))
            } else {
                data.table::set(segs, i, c("new_major", "new_minor", "status"), list(1L, 1L, "loh_rejected_total_raised"))
            }
        } else if (mn > 0 && llr < -llr_min && tot > max_total_het_test) {
            ## amplified segments: a single minor copy (f <= 1/(max_total_het_test + 1)) is not separable from LOH
            data.table::set(segs, i, "status", "high_cn_untested")
        } else if (mn > 0 && llr < -llr_min) {
            if (fix_het) data.table::set(segs, i, c("new_major", "new_minor", "status"), list(as.integer(tot), 0L, "loh_by_baf"))
            else data.table::set(segs, i, "status", "loh_by_baf_not_applied")
        } else {
            st <- if (mn == 0) { if (llr <= -llr_min) "loh_supported" else "inconclusive" } else { if (llr >= llr_min) "het_supported" else "inconclusive" }
            data.table::set(segs, i, "status", st)
        }
    }
    list(segs = segs[], params = list(r = fit$r, rho = fit$rho, sex = sex, n_ref_sites = fit$n_sites, w = w, err = err,
                                      min_sites = min_sites, min_sites_reject = min_sites_reject, llr_min = llr_min, bin_width = bin_width))
}

#' @name sc_allelic_gg_segments
#' @title sc_allelic_gg_segments
#' @description
#'
#' One row per original node of a phased (haploid) allelic gGraph: chr, start,
#' end, the two haploid node ids and their CN as major (larger) / minor.
#'
#' @param gg allelic gGraph (nodes with og.node.id, allele, cn) or rds path
#' @return data.table(seg_id, chr, start, end, major, minor, hi_node, lo_node)
#' @export
#' @author Stanley Clarke
sc_allelic_gg_segments <- function(gg) {
    if (is.character(gg)) gg <- readRDS(gg)
    nd <- data.table::as.data.table(gg$nodes$dt)
    if (!all(c("og.node.id", "cn") %in% names(nd))) stop("allelic gGraph needs og.node.id and cn on nodes")
    nd[, cn := as.numeric(cn)]
    nd[is.na(cn), cn := 0]
    nd[, .(chr = as.character(seqnames[1]), start = start[1], end = end[1],
           major = as.integer(round(max(cn))), minor = as.integer(round(min(cn))),
           hi_node = node.id[which.max(cn)], lo_node = node.id[if (.N > 1) setdiff(seq_len(.N), which.max(cn))[1] else 1L]),
       by = .(seg_id = og.node.id)]
}

#' @name sc_allelic_fix_gg
#' @title sc_allelic_fix_gg
#' @description
#'
#' Apply sc_allelic_baf_fix() to an allelic gGraph: the haploid node CNs of
#' corrected segments are replaced and every node gets a baf_fix annotation
#' (status, scaled LLR, sites) and a logical hemizygous flag (male chrX/Y
#' outside the PARs: major = total, minor = 0, not LOH). Edge CNs are left as
#' solved by the LP.
#'
#' @param gg allelic gGraph or rds path
#' @param hets data.table(chr, pos, ref, alt) or a path for sc_read_het_counts()
#' @param ... passed to sc_allelic_baf_fix()
#' @return list(gg, segs, params)
#' @export
#' @author Stanley Clarke
sc_allelic_fix_gg <- function(gg, hets, ...) {
    if (is.character(gg)) gg <- readRDS(gg)
    if (is.character(hets)) hets <- sc_read_het_counts(hets)
    gg <- gg$copy
    segs <- sc_allelic_gg_segments(gg)
    res <- sc_allelic_baf_fix(segs[, .(seg_id, chr, start, end, major, minor)], hets, ...)
    s <- merge(segs[, .(seg_id, hi_node, lo_node)], res$segs, by = "seg_id")
    cn <- as.numeric(gg$nodes$dt$cn)
    ids <- gg$nodes$dt$node.id
    ch <- s[new_major != major | new_minor != minor]
    if (nrow(ch)) {
        cn[match(ch$hi_node, ids)] <- ch$new_major
        cn[match(ch$lo_node, ids)] <- ch$new_minor
    }
    lab <- s[, sprintf("%s; llr=%s; sites=%d", status, ifelse(is.na(llr), "NA", sprintf("%.1f", llr)), n_sites)]
    ann <- lab[match(gg$nodes$dt$og.node.id, s$seg_id)]
    hemi <- s$hemi[match(gg$nodes$dt$og.node.id, s$seg_id)]
    gg$nodes$mark(cn = cn, baf_fix = ann, hemizygous = !is.na(hemi) & hemi)
    list(gg = gg, segs = res$segs, params = res$params)
}
