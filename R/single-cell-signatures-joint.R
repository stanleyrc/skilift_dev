## Joint (patient-level) COSMIC SBS refit of single-cell SNV site sets.
##
## sc_fit_sbs_signatures() fits every site set independently against all of COSMIC,
## so small sets (subclonal, private, one clone) pick spurious signatures. The
## functions here choose ONE signature set per patient, refit every set on it
## with NNLS, and put bootstrap intervals and a subclonal-vs-truncal test on the
## exposure fractions. Base R + nnls only (no SigProfiler container needed).

#' @name sc_read_cosmic_reference
#' @title sc_read_cosmic_reference
#' @description
#'
#' Reads a COSMIC SBS96 signature table (MutationType column then one column per
#' signature, e.g. the Assignment_Solution_Signatures.txt that SigProfilerAssignment
#' cosmic_fit writes, which is the COSMIC v3.4 GRCh38 reference it used) as a
#' 96 x n matrix in SC_SBS96 order with columns summing to 1.
#'
#' @param path tab-separated signature table
#' @return numeric matrix, rownames SC_SBS96, colnames signature names
#' @export
#' @author Stanley Clarke
sc_read_cosmic_reference <- function(path) {
    d <- data.table::fread(path)
    m <- as.matrix(d[, -1, with = FALSE])
    rownames(m) <- d[[1]]
    if (!all(SC_SBS96 %in% rownames(m))) stop("reference lacks SBS96 channels: ", path)
    m <- m[SC_SBS96, , drop = FALSE]
    sweep(m, 2, colSums(m), "/")
}

## base matrix product: gUtils masks %*% when skilift is attached
.sj_mm <- function(a, b) base::`%*%`(a, b)

.sj_cosine <- function(a, b) {
    d <- sqrt(sum(a^2)) * sqrt(sum(b^2))
    if (d == 0) NA_real_ else sum(a * b) / d
}

## NNLS exposures (in mutations) of a 96-vector of counts on the columns idx of S
.sj_nnls <- function(counts, S, idx) {
    e <- numeric(ncol(S))
    if (!length(idx)) return(e)
    e[idx] <- nnls::nnls(S[, idx, drop = FALSE], counts)$x
    e
}

## multinomial log-likelihood of counts under the normalised NNLS reconstruction
.sj_loglik <- function(counts, S, e) {
    rec <- as.vector(.sj_mm(S, e))
    if (sum(rec) <= 0) return(-Inf)
    p <- pmax(rec / sum(rec), 1e-12)
    sum(counts * log(p))
}

#' @name sc_sig_forward_select
#' @title sc_sig_forward_select
#' @description
#'
#' Greedy forward selection of signatures on top of a base set for ONE catalog.
#' Each round adds the candidate signature with the largest likelihood gain
#' (multinomial log-likelihood of the NNLS reconstruction) if its likelihood-ratio
#' test (1 df) passes alpha / number of candidates and it takes at least min_frac
#' of the exposure.
#'
#' @param counts 96-vector of channel counts
#' @param S 96 x n signature matrix
#' @param base column indices (or names) of the base signature set
#' @param alpha family-wise level of the likelihood-ratio test
#' @param min_frac smallest exposure fraction of an added signature
#' @param max_add most signatures added
#' @return list(added = integer indices in selection order, gain = data.table of the
#'   steps: signature, delta_loglik, p, frac, cosine_gain)
#' @export
#' @author Stanley Clarke
sc_sig_forward_select <- function(counts, S, base, alpha = 0.05, min_frac = 0.05, max_add = 4L) {
    if (is.character(base)) base <- match(base, colnames(S))
    cur <- base
    steps <- list()
    for (r in seq_len(max_add)) {
        cands <- setdiff(seq_len(ncol(S)), cur)
        e0 <- .sj_nnls(counts, S, cur)
        ll0 <- .sj_loglik(counts, S, e0)
        cos0 <- .sj_cosine(counts, as.vector(.sj_mm(S, e0)))
        best <- NULL
        for (cd in cands) {
            e1 <- .sj_nnls(counts, S, c(cur, cd))
            if (sum(e1) <= 0 || e1[cd] / sum(e1) < min_frac) next
            ll1 <- .sj_loglik(counts, S, e1)
            if (is.null(best) || ll1 > best$ll) best <- list(cd = cd, ll = ll1, frac = e1[cd] / sum(e1),
                                                           cos = .sj_cosine(counts, as.vector(.sj_mm(S, e1))))
        }
        if (is.null(best)) break
        delta <- 2 * (best$ll - ll0)
        p <- stats::pchisq(delta, 1, lower.tail = FALSE)
        if (!is.finite(delta) || p >= alpha / length(cands)) break
        cur <- c(cur, best$cd)
        steps[[length(steps) + 1L]] <- data.table::data.table(
            signature = colnames(S)[best$cd], index = best$cd, delta_loglik = delta / 2, p = p,
            frac = best$frac, cosine_gain = best$cos - cos0)
    }
    list(added = if (length(steps)) data.table::rbindlist(steps)$index else integer(0),
         gain = if (length(steps)) data.table::rbindlist(steps) else NULL)
}

## bootstrap resample of a catalog = resampling its mutations
.sj_resample <- function(counts) {
    n <- sum(counts)
    as.numeric(stats::rmultinom(1L, n, counts / n))
}

#' @name sc_sig_patient_set
#' @title sc_sig_patient_set
#' @description
#'
#' Chooses the patient's active COSMIC signatures: the base set (SigProfiler
#' cosmic_fit of all sites) plus signatures that are strongly supported in a subset
#' of sites, improve that subset's NNLS fit significantly (forward selection, see
#' sc_sig_forward_select) and are selected in at least min_stability of n_boot
#' bootstrap resamples of the subset. Records why each signature is in the set.
#'
#' @param catalogs named list of 96-vectors (the subsets that may add signatures)
#' @param S COSMIC matrix (sc_read_cosmic_reference())
#' @param base_signatures names of the base set
#' @param all_catalog 96-vector of the all-sites catalog (bootstrap support of the base set)
#' @param n_boot bootstrap resamples per subset
#' @param min_stability fraction of resamples that must select an added signature
#' @param min_select smallest subset (mutations) allowed to add signatures
#' @param alpha,min_frac,max_add see sc_sig_forward_select()
#' @param seed random seed
#' @return list(signatures = character, reasons = data.table(signature, source, set, ...))
#' @export
#' @author Stanley Clarke
sc_sig_patient_set <- function(catalogs, S, base_signatures, all_catalog = NULL, n_boot = 200L,
                               min_stability = 0.8, min_select = 100L, alpha = 0.05, min_frac = 0.05,
                               max_add = 4L, seed = 1L) {
    set.seed(seed)
    base_signatures <- intersect(base_signatures, colnames(S))
    base <- match(base_signatures, colnames(S))
    reasons <- data.table::data.table(
        signature = base_signatures, source = "base", set = "all",
        detail = "active in the SigProfilerAssignment cosmic_fit of all sites", stability = NA_real_,
        delta_loglik = NA_real_, p = NA_real_, frac = NA_real_, cosine_gain = NA_real_)
    if (!is.null(all_catalog) && sum(all_catalog) >= min_select && length(base)) {
        hit <- matrix(FALSE, n_boot, length(base))
        for (b in seq_len(n_boot)) {
            e <- .sj_nnls(.sj_resample(all_catalog), S, base)
            hit[b, ] <- e[base] / max(sum(e), 1e-12) >= 0.01
        }
        reasons[source == "base", stability := colMeans(hit)]
    }
    for (nm in names(catalogs)) {
        ct <- catalogs[[nm]]
        if (sum(ct) < min_select) next
        sel <- sc_sig_forward_select(ct, S, base, alpha, min_frac, max_add)
        if (!length(sel$added)) next
        count <- integer(ncol(S))
        for (b in seq_len(n_boot)) {
            rs <- sc_sig_forward_select(.sj_resample(ct), S, base, alpha, min_frac, max_add)
            count[rs$added] <- count[rs$added] + 1L
        }
        for (k in seq_len(nrow(sel$gain))) {
            g <- sel$gain[k]
            stab <- count[g$index] / n_boot
            if (stab >= min_stability) {
                reasons <- rbind(reasons, data.table::data.table(
                    signature = g$signature, source = "subset", set = nm,
                    detail = sprintf("%s: %.0f%% of %d mutations, delta logLik %.1f, LRT p %.1e, selected in %.0f%% of %d bootstraps",
                                     nm, 100 * g$frac, sum(ct), g$delta_loglik, g$p, 100 * stab, n_boot),
                    stability = stab, delta_loglik = g$delta_loglik, p = g$p, frac = g$frac, cosine_gain = g$cosine_gain))
            }
        }
    }
    reasons <- reasons[!duplicated(signature)]
    list(signatures = reasons$signature, reasons = reasons[])
}

#' @name sc_sig_refit_bootstrap
#' @title sc_sig_refit_bootstrap
#' @description
#'
#' NNLS refit of a catalog on a signature set with bootstrap (mutation resampling)
#' exposure fractions and reconstruction cosine.
#'
#' @param counts 96-vector
#' @param S 96 x n COSMIC matrix
#' @param signatures names of the patient signature set
#' @param n_boot resamples
#' @return list(exposure, fraction, cosine, boot_fraction matrix (n_boot x signatures), boot_cosine, reconstruction)
#' @export
#' @author Stanley Clarke
sc_sig_refit_bootstrap <- function(counts, S, signatures, n_boot = 1000L) {
    idx <- match(signatures, colnames(S))
    S1 <- S[, idx, drop = FALSE]
    fit <- function(ct) {
        e <- nnls::nnls(S1, ct)$x
        list(e = e, cos = .sj_cosine(ct, as.vector(.sj_mm(S1, e))))
    }
    f0 <- fit(counts)
    bf <- matrix(NA_real_, n_boot, length(idx), dimnames = list(NULL, signatures))
    bc <- numeric(n_boot)
    for (b in seq_len(n_boot)) {
        fb <- fit(.sj_resample(counts))
        bf[b, ] <- if (sum(fb$e) > 0) fb$e / sum(fb$e) else NA_real_
        bc[b] <- fb$cos
    }
    list(exposure = stats::setNames(f0$e, signatures),
         fraction = stats::setNames(if (sum(f0$e) > 0) f0$e / sum(f0$e) else rep(NA_real_, length(idx)), signatures),
         cosine = f0$cos, boot_fraction = bf, boot_cosine = bc,
         reconstruction = as.vector(.sj_mm(S1, f0$e)))
}

#' @name sc_sig_contrast
#' @title sc_sig_contrast
#' @description
#'
#' Per-signature difference of exposure fraction between two refit sets from their
#' independent bootstrap replicates: difference, percentile 95\% CI and a two-sided
#' bootstrap p (with +1 correction), BH-adjusted across signatures.
#'
#' @param a,b outputs of sc_sig_refit_bootstrap() (a minus b)
#' @return data.table (signature, frac_a, frac_b, diff, lo, hi, p, q)
#' @export
#' @author Stanley Clarke
sc_sig_contrast <- function(a, b) {
    sigs <- colnames(a$boot_fraction)
    nb <- min(nrow(a$boot_fraction), nrow(b$boot_fraction))
    out <- data.table::rbindlist(lapply(sigs, function(s) {
        d <- a$boot_fraction[seq_len(nb), s] - b$boot_fraction[seq_len(nb), s]
        d <- d[is.finite(d)]
        qs <- stats::quantile(d, c(0.025, 0.975), names = FALSE)
        p <- min(1, 2 * (min(sum(d <= 0), sum(d >= 0)) + 1) / (length(d) + 1))
        data.table::data.table(signature = s, frac_a = a$fraction[[s]], frac_b = b$fraction[[s]],
                               diff = a$fraction[[s]] - b$fraction[[s]], lo = qs[1], hi = qs[2], p = p)
    }))
    out[, q := stats::p.adjust(p, "BH")]
    out[]
}

.sj_stats <- function(obs, rec) {
    p <- obs / sum(obs); q <- pmax(rec / sum(rec), 1e-12)
    list(cosine = round(.sj_cosine(obs, rec), 4),
         l1_pct = round(100 * sum(abs(obs - rec)) / sum(obs), 2),
         l2_pct = round(100 * sqrt(sum((obs - rec)^2)) / sqrt(sum(obs^2)), 2),
         kl = round(sum(ifelse(p > 0, p * log(p / q), 0)), 4),
         correlation = round(suppressWarnings(stats::cor(obs, rec)), 4))
}

#' @name sc_refit_signatures_joint
#' @title sc_refit_signatures_joint
#' @description
#'
#' Joint signature refit of a patient's SNV site sets (see the file header):
#' patient signature set (sc_sig_patient_set), NNLS refit + bootstrap of every
#' set (sc_sig_refit_bootstrap), subclonal / private / clone-specific / clade vs
#' truncal contrasts (sc_sig_contrast). The returned list is written as
#' signatures_joint.json; sets[[i]]$activities has the same shape as the legacy
#' signatures.json (signature, activity = exposure in mutations) so the browser's
#' per-site assignment can use it unchanged.
#'
#' @param sets named list of SBS96 context vectors (as sc_fit_sbs_signatures())
#' @param S COSMIC matrix (sc_read_cosmic_reference())
#' @param base_signatures signatures of the cosmic_fit of the "all" set
#' @param set_info optional named list (per set) of list(group, definition, note) labels
#' @param subset_sets names of sets that may add signatures to the patient set
#' @param legacy optional legacy signatures.json (list) to report dropped signatures
#' @param min_mutations smallest set refitted
#' @param n_boot bootstrap resamples of the refits and tests
#' @param n_boot_select bootstrap resamples of the signature selection
#' @param seed seed
#' @return list ready for jsonlite::write_json(auto_unbox = TRUE)
#' @export
#' @author Stanley Clarke
sc_refit_signatures_joint <- function(sets, S, base_signatures, set_info = list(), subset_sets = NULL,
                                      legacy = NULL, min_mutations = 20L, n_boot = 1000L,
                                      n_boot_select = 200L, seed = 1L) {
    count_of <- function(x) as.numeric(table(factor(x, levels = SC_SBS96)))
    catalogs <- lapply(sets, count_of)
    n <- vapply(catalogs, sum, numeric(1))
    catalogs <- catalogs[n >= min_mutations]
    if (!("all" %in% names(catalogs))) stop("sets needs an 'all' set")
    subset_sets <- intersect(if (is.null(subset_sets)) setdiff(names(catalogs), "all") else subset_sets, names(catalogs))
    ps <- sc_sig_patient_set(catalogs[subset_sets], S, base_signatures, all_catalog = catalogs$all,
                             n_boot = n_boot_select, seed = seed)
    sigs <- ps$signatures
    set.seed(seed + 1L)
    fits <- lapply(catalogs, sc_sig_refit_bootstrap, S = S, signatures = sigs, n_boot = n_boot)
    ci <- function(x) stats::quantile(x, c(0.025, 0.975), names = FALSE, na.rm = TRUE)
    num <- function(x, d = 4) round(as.numeric(x), d)
    out_sets <- lapply(names(fits), function(nm) {
        f <- fits[[nm]]; ct <- catalogs[[nm]]
        cq <- ci(f$boot_cosine)
        acts <- lapply(order(-f$exposure), function(i) {
            s <- sigs[i]
            if (f$exposure[i] <= 0) return(NULL)
            q <- ci(f$boot_fraction[, s])
            list(signature = s, activity = num(f$exposure[i], 2), fraction = num(f$fraction[i]), lo = num(q[1]), hi = num(q[2]))
        })
        info <- set_info[[nm]]
        list(name = nm, n = as.integer(sum(ct)), group = info$group, definition = info$definition, note = info$note,
             activities = Filter(Negate(is.null), acts),
             cosine = list(est = num(f$cosine), lo = num(cq[1]), hi = num(cq[2])),
             counts = as.integer(ct), reconstruction = round(f$reconstruction, 2),
             stats = .sj_stats(ct, f$reconstruction))
    })
    names(out_sets) <- names(fits)
    contrast_table <- function(a, b, label, kind) {
        if (!(a %in% names(fits)) || !(b %in% names(fits))) return(NULL)
        tb <- sc_sig_contrast(fits[[a]], fits[[b]])
        list(contrast = label, kind = kind, a = a, b = b, n_a = as.integer(sum(catalogs[[a]])), n_b = as.integer(sum(catalogs[[b]])),
             rows = lapply(seq_len(nrow(tb)), function(i) as.list(tb[i, .(signature, frac_a = num(frac_a), frac_b = num(frac_b),
                                                                         diff = num(diff), lo = num(lo), hi = num(hi),
                                                                         p = signif(p, 3), q = signif(q, 3))])))
    }
    contrasts <- Filter(Negate(is.null), c(
        list(contrast_table("subclonal", "truncal", "subclonal vs truncal", "tree"),
             contrast_table("private", "truncal", "private vs truncal", "tree")),
        lapply(grep("^clone-specific: |^clade: ", names(fits), value = TRUE), function(nm)
            contrast_table(nm, "truncal", paste(nm, "vs truncal"), if (grepl("^clade", nm)) "clade" else "clone"))))
    legacy_dropped <- NULL
    if (!is.null(legacy)) {
        seen <- list()
        for (s in legacy$sets) for (a in s$activities) {
            sg <- a$signature
            if (!(sg %in% sigs)) seen[[sg]] <- c(seen[[sg]], s$name)
        }
        legacy_dropped <- lapply(names(seen), function(sg) list(signature = sg, legacy_sets = unique(seen[[sg]])))
    }
    list(cosmic_version = "3.4", genome = "GRCh38",
         method = "patient signature set (cosmic_fit of all sites + bootstrap-stable subset additions), NNLS refit of every set, bootstrap CIs",
         n_boot = n_boot, n_boot_select = n_boot_select, seed = seed,
         signature_set = lapply(seq_len(nrow(ps$reasons)), function(i) as.list(ps$reasons[i])),
         sets = unname(out_sets), channels = SC_SBS96, contrasts = contrasts, legacy_dropped = legacy_dropped)
}

#' @name sc_signatures_joint_cohort
#' @title sc_signatures_joint_cohort
#' @description
#'
#' Cohort summary of per-patient signatures_joint.json lists: for every patient and
#' signature the subclonal-vs-truncal and private-vs-truncal difference with CI, p and q,
#' plus each patient's signature set.
#'
#' @param per_patient named list patient -> signatures_joint list
#' @return list(patients, contrasts = list of rows with patient)
#' @export
#' @author Stanley Clarke
sc_signatures_joint_cohort <- function(per_patient) {
    rows <- list()
    for (P in names(per_patient)) for (ct in per_patient[[P]]$contrasts) {
        if (!(ct$kind %in% "tree")) next
        for (r in ct$rows) rows[[length(rows) + 1L]] <- c(list(patient = P, contrast = ct$contrast, n_a = ct$n_a, n_b = ct$n_b), r)
    }
    list(method = "per-signature exposure-fraction difference, bootstrap CI and p (see each patient's signatures_joint.json)",
         patients = lapply(names(per_patient), function(P) list(
             patient = P, signatures = vapply(per_patient[[P]]$signature_set, function(s) s$signature, ""),
             legacy_dropped = vapply(if (is.null(per_patient[[P]]$legacy_dropped)) list() else per_patient[[P]]$legacy_dropped, function(s) s$signature, ""))),
         contrasts = rows)
}
