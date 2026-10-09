#' @name sc_export_rna_fusions
#' @title sc_export_rna_fusions
#' @description
#'
#' RNA fusions of one patient for gOS: the per-cell Arriba calls (STAR chimeric
#' alignment of each cell, fusions.tsv) with confidence in `keep_confidence`,
#' aggregated per fusion (gene pair + breakpoints) with the read support of every
#' cell, the gOS cell id of cells that also have DNA, and whether the patient's
#' DNA filtered events hold the same gene pair. Most single-cell calls are
#' one-cell artefacts, so a fusion is exported when it matches a DNA fusion, is
#' an Arriba known / recurrent fusion (tags), is seen in `min_cells` cells, or is
#' a high-confidence call of one cell with `single_min_reads` reads. With `slices`, every carrier
#' cell gets an IGV slice of its RNA reads within `pad` bp of its fusion
#' breakpoints (rna/reads/<rna_id>.bam; supporting reads tagged ZF:i:1).
#' Every fusion gets a tier (see sc_rna_fusion_tier: 1 known / actionable
#' driver fusion, 2 cancer-gene fusion with functional evidence, 3 other) with
#' its reasons, the cancer genes involved, Arriba's retained protein domains
#' and transcripts. Writes data/<P>/rna/fusions.json (gos-sc-rna-fusions/1).
#'
#' @param patient patient id
#' @param cell_dirs named character: rna_id -> per-cell folder (fusions.tsv, Aligned.sorted.bam)
#' @param out_dir patient folder of the gOS dataset (data/<P>)
#' @param cell_map named character rna_id -> gOS cell id (NA when the cell has no DNA)
#' @param dna_events optional path to the patient's filtered.events.json
#' @param keep_confidence Arriba confidence levels kept
#' @param slices write per-cell RNA BAM slices
#' @param min_cells fusions seen in at least this many cells are kept
#' @param single_min_reads a fusion of one cell is kept when high confidence with at least this many reads
#' @param max_fusions cap on exported fusions (DNA-matched ones always kept)
#' @param pad bp around each breakpoint in the slices
#' @param cores parallel cells
#' @param refs reference gene lists for the tiers (sc_rna_fusion_references())
#' @return the fusions data.table (invisibly)
#' @export
#' @author Stanley Clarke
sc_export_rna_fusions <- function(patient, cell_dirs, out_dir, cell_map = NULL, dna_events = NULL,
                                  keep_confidence = c("high", "medium"), min_cells = 2, single_min_reads = 10,
                                  max_fusions = 2000, slices = TRUE, pad = 300, cores = 4,
                                  refs = sc_rna_fusion_references()) {
    files <- file.path(cell_dirs, "fusions.tsv")
    keep <- file.exists(files)
    calls <- data.table::rbindlist(parallel::mclapply(which(keep), function(i) {
        x <- data.table::fread(files[i], sep = "\t", quote = "", colClasses = "character", showProgress = FALSE)
        if (!nrow(x)) return(NULL)
        data.table::setnames(x, sub("^#", "", names(x)))
        x[, rna_id := names(cell_dirs)[i]]
        x
    }, mc.cores = cores), fill = TRUE)
    rna_dir <- file.path(out_dir, "rna")
    dir.create(rna_dir, showWarnings = FALSE, recursive = TRUE)
    empty <- list(format = "gos-sc-rna-fusions/1", patient = patient, n_cells_rna = sum(keep), fusions = list())
    if (!nrow(calls)) {
        jsonlite::write_json(empty, file.path(rna_dir, "fusions.json"), auto_unbox = TRUE)
        return(invisible(data.table::data.table()))
    }
    calls <- calls[confidence %in% keep_confidence]
    num <- function(v) suppressWarnings(as.integer(v))
    calls[, `:=`(split1 = num(split_reads1), split2 = num(split_reads2), discordant = num(discordant_mates))]
    calls[, fusion_id := paste0(gene1, "::", gene2, "|", breakpoint1, "|", breakpoint2)]
    calls[, cell_id := if (is.null(cell_map)) NA_character_ else unname(cell_map[rna_id])]
    rank <- c(high = 3L, medium = 2L, low = 1L)

    ## DNA fusions of the patient (gene pairs, either order)
    dna_pairs <- character(0)
    dna_tier_of <- integer(0)
    if (!is.null(dna_events) && file.exists(dna_events)) {
        ev <- data.table::as.data.table(jsonlite::fromJSON(dna_events))
        if ("vartype" %in% names(ev)) {
            dfx <- ev[vartype %in% c("fusion", "outframe_fusion")]
            dna_pairs <- unique(dfx$gene)
            if ("Tier" %in% names(dfx) && nrow(dfx)) {
                tt <- dfx[, .(tier = suppressWarnings(min(as.integer(Tier), na.rm = TRUE))), by = gene]
                dna_tier_of <- stats::setNames(ifelse(is.finite(tt$tier), tt$tier, NA_integer_), tt$gene)
            }
        }
    }
    pair_key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "::")
    dna_keys <- if (length(dna_pairs)) {
        g <- data.table::tstrsplit(dna_pairs, "::", fixed = TRUE)
        stats::setNames(dna_pairs, pair_key(g[[1]], g[[2]]))
    } else character(0)
    clean_gene <- function(g) sub("\\(.*$", "", sub(",.*$", "", g))   ## intergenic "A(123),B(456)" -> A
    ## most frequent informative value of an Arriba column across the cells of a fusion ("." when none)
    top_value <- function(v) {
        v <- v[!is.na(v) & v != "." & nzchar(v)]
        if (!length(v)) return(".")
        names(sort(table(v), decreasing = TRUE))[1]
    }
    for (col in c("retained_protein_domains", "transcript_id1", "transcript_id2"))
        if (!col %in% names(calls)) calls[, (col) := "."]

    fus <- calls[, .(gene1 = gene1[1], gene2 = gene2[1], breakpoint1 = breakpoint1[1], breakpoint2 = breakpoint2[1],
                     strand1 = `strand1(gene/fusion)`[1], strand2 = `strand2(gene/fusion)`[1],
                     site1 = site1[1], site2 = site2[1], type = type[1], reading_frame = top_value(reading_frame),
                     confidence = names(rank)[match(max(rank[confidence]), rank)],
                     n_cells = data.table::uniqueN(rna_id), n_cells_dna = data.table::uniqueN(cell_id[!is.na(cell_id)]),
                     split_reads = sum(split1 + split2, na.rm = TRUE), discordant_mates = sum(discordant, na.rm = TRUE),
                     retained_domains = top_value(retained_protein_domains),
                     transcript1 = top_value(transcript_id1), transcript2 = top_value(transcript_id2)),
                 by = .(id = fusion_id)]
    fus[, dna_gene := unname(dna_keys[pair_key(clean_gene(gene1), clean_gene(gene2))])]
    fus[, dna_tier := if (length(dna_tier_of)) unname(dna_tier_of[dna_gene]) else NA_integer_]
    known <- unique(calls[!is.na(tags) & tags != "." & grepl("Mitelman|known|COSMIC|CCLE", tags, ignore.case = TRUE)]$fusion_id)
    fus[, known := id %in% known | paste(clean_gene(gene1), clean_gene(gene2), sep = "::") %in% refs$known_pairs]
    n_all <- nrow(fus)
    fus <- fus[!is.na(dna_gene) | known | n_cells >= min_cells |
               (confidence == "high" & split_reads + discordant_mates >= single_min_reads)]
    data.table::setorder(fus, -n_cells, -split_reads)
    fus <- fus[!is.na(dna_gene) | seq_len(.N) <= max_fusions]
    message(sprintf("sc_export_rna_fusions %s: %d of %d fusions kept (%d DNA-matched, %d known)", patient, nrow(fus), n_all,
                    sum(!is.na(fus$dna_gene)), sum(fus$known)))
    calls <- calls[fusion_id %in% fus$id]
    tiers <- sc_rna_fusion_tier(fus, refs, n_cells_rna = sum(keep))
    message(sprintf("sc_export_rna_fusions %s: tier 1 %d, tier 2 %d, tier 3 %d", patient,
                    sum(tiers$tier == 1L), sum(tiers$tier == 2L), sum(tiers$tier == 3L)))

    ## per-cell RNA slices around the breakpoints, supporting reads tagged ZF:i:1
    if (slices) {
        reads_dir <- file.path(rna_dir, "reads")
        dir.create(reads_dir, showWarnings = FALSE)
        by_cell <- split(calls, calls$rna_id)
        script <- system.file("extdata", "scripts", "sc_rna_fusion_slice.py", package = "Skilift")
        if (!nzchar(script)) script <- file.path(getNamespaceInfo("Skilift", "path"), "inst", "extdata", "scripts", "sc_rna_fusion_slice.py")
        py <- Sys.getenv("SC_PYTHON", "/gpfs/commons/groups/imielinski_lab/Software/miniforge3/envs/mskilab_ne1/bin/python")
        ok <- parallel::mclapply(names(by_cell), function(r) {
            bam <- file.path(cell_dirs[[r]], "Aligned.sorted.bam")
            if (!file.exists(paste0(bam, ".bai"))) return(FALSE)
            x <- by_cell[[r]]
            regions <- unique(unlist(lapply(c(x$breakpoint1, x$breakpoint2), function(b) {
                p <- strsplit(b, ":", fixed = TRUE)[[1]]
                pos <- as.integer(p[2])
                paste0(p[1], ":", max(1, pos - pad), "-", pos + pad)
            })))
            ids <- unique(unlist(strsplit(x$read_identifiers[!is.na(x$read_identifiers) & x$read_identifiers != "."], ",", fixed = TRUE)))
            idf <- tempfile(fileext = ".txt"); on.exit(unlink(idf))
            writeLines(ids, idf)
            out <- file.path(reads_dir, paste0(r, ".bam"))
            rc <- system2(py, c(shQuote(script), shQuote(bam), shQuote(out), shQuote(idf), shQuote(paste(regions, collapse = ","))),
                          stdout = FALSE, stderr = FALSE)
            rc == 0 && file.exists(paste0(out, ".bai"))
        }, mc.cores = cores, mc.preschedule = FALSE)
        sliced <- names(by_cell)[unlist(ok) %in% TRUE]
        message(sprintf("sc_export_rna_fusions %s: RNA slices for %d of %d carrier cells", patient, length(sliced), length(by_cell)))
    } else sliced <- character(0)

    cells_of <- split(calls, calls$fusion_id)
    out <- lapply(seq_len(nrow(fus)), function(i) {
        f <- fus[i]
        cc <- cells_of[[f$id]][order(-(split1 + split2 + discordant))]
        tr <- tiers[i]
        c(as.list(f[, !c("dna_gene", "dna_tier")]),
          list(tier = tr$tier, tier_label = tr$tier_label, tier_reasons = as.list(tr$tier_reasons[[1]]),
               cancer_genes = tr$cancer_genes[[1]], oncokb_level = if (is.na(tr$oncokb_level)) NULL else tr$oncokb_level,
               productive = tr$productive, cell_fraction = tr$cell_fraction),
          list(dna_match = if (is.na(f$dna_gene)) NULL else list(kind = "gene_pair", event_gene = f$dna_gene,
                                                                 tier = if (is.na(f$dna_tier)) NULL else f$dna_tier),
               cells = lapply(seq_len(nrow(cc)), function(k) list(
                   rna_id = cc$rna_id[k], cell_id = if (is.na(cc$cell_id[k])) NULL else cc$cell_id[k],
                   split1 = cc$split1[k], split2 = cc$split2[k], discordant = cc$discordant[k], confidence = cc$confidence[k],
                   bam = if (cc$rna_id[k] %in% sliced) paste0("rna/reads/", cc$rna_id[k], ".bam") else NULL))))
    })
    res <- c(empty[c("format", "patient", "n_cells_rna")], list(fusions = out))
    jsonlite::write_json(res, file.path(rna_dir, "fusions.json"), auto_unbox = TRUE, null = "null", na = "null", digits = NA)
    invisible(fus)
}

#' @name sc_rna_fusion_references
#' @title sc_rna_fusion_references
#' @description
#'
#' Reference lists for the RNA fusion tiers: Arriba's known fusion gene pairs
#' (5'::3', from the known_fusions database), cancer genes with their role
#' (OncoKB cancer gene list, then the COSMIC Cancer Gene Census) and OncoKB's
#' fusion biomarkers (genes whose fusions carry a level, e.g. NTRK1 "Fusions",
#' and named pairs such as FGFR3-TACC3) with their best level (1-4 therapeutic,
#' then Dx, then Px). Missing files give empty lists.
#'
#' @param known_fusions Arriba known_fusions_*.tsv.gz
#' @param oncokb_genes OncoKB cancer gene list (rds data.table with Hugo_Symbol, Role)
#' @param oncokb_biomarkers OncoKB biomarker-drug associations tsv (Level, Gene, Alterations)
#' @param cgc COSMIC Cancer Gene Census csv
#' @return list(known_pairs, cancer_roles, oncokb_genes, oncokb_pairs)
#' @export
#' @author Stanley Clarke
sc_rna_fusion_references <- function(
    known_fusions = "/nfs/sw/easybuild/software/arriba/2.5.0/database/known_fusions_hg38_GRCh38_v2.5.0.tsv.gz",
    oncokb_genes = "/gpfs/commons/groups/imielinski_lab/DB/OncoKB/OncoKB_cancer_genes.rds",
    oncokb_biomarkers = "/gpfs/commons/groups/imielinski_lab/DB/OncoKB/oncokb_biomarker_drug_associations.tsv",
    cgc = "/gpfs/commons/groups/imielinski_lab/DB/COSMIC/v97_GRCh38/cancer_gene_census.csv") {
    ok <- function(f) is.character(f) && length(f) == 1 && !is.na(f) && file.exists(f)
    known_pairs <- character(0)
    if (ok(known_fusions)) {
        con <- gzfile(known_fusions); l <- readLines(con, warn = FALSE); close(con)
        p <- strsplit(sub("^#", "", grep("^#", l, value = TRUE)), "\t", fixed = TRUE)
        p <- p[vapply(p, length, 1L) >= 2]
        known_pairs <- unique(toupper(vapply(p, function(x) paste(x[1], x[2], sep = "::"), "")))
    }
    roles <- character(0)
    if (ok(oncokb_genes)) {
        x <- data.table::as.data.table(readRDS(oncokb_genes))
        r <- tolower(ifelse(is.na(x$Role) | !nzchar(x$Role), "cancer gene", x$Role))
        roles <- stats::setNames(r, toupper(x$Hugo_Symbol))
    }
    if (ok(cgc)) {
        y <- data.table::fread(cgc, showProgress = FALSE)
        g <- toupper(y[["Gene Symbol"]])
        r <- y[["Role in Cancer"]]
        r <- ifelse(is.na(r) | !nzchar(r), "cancer gene", r)
        new <- !g %in% names(roles)
        roles <- c(roles, stats::setNames(r[new], g[new]))
    }
    okb_genes <- okb_pairs <- character(0)
    if (ok(oncokb_biomarkers)) {
        b <- data.table::fread(oncokb_biomarkers, sep = "\t", showProgress = FALSE)
        data.table::setnames(b, 1:3, c("level", "gene", "alt"))
        b <- b[grepl("Fusion", alt) & !grepl("excluding Fusions", alt)]
        rank <- function(l) ifelse(grepl("^\\d$", l), suppressWarnings(as.integer(l)),
                            ifelse(grepl("^Dx", l), 4L + suppressWarnings(as.integer(sub("Dx", "", l))),
                                   8L + suppressWarnings(as.integer(sub("Px", "", l)))))
        best <- function(gene, level) {
            d <- data.table::data.table(gene = gene, level = level, r = rank(level))[order(r)]
            d <- d[!duplicated(gene)]
            stats::setNames(d$level, d$gene)
        }
        gl <- b[grepl("(^|, )Fusions", alt)]
        okb_genes <- best(toupper(gl$gene), as.character(gl$level))
        pr <- b[, .(pair = unlist(regmatches(alt, gregexpr("[A-Za-z0-9]+-[A-Za-z0-9]+(?= Fusion)", alt, perl = TRUE)))), by = .(level = as.character(level))]
        if (nrow(pr)) okb_pairs <- best(toupper(sub("-", "::", pr$pair, fixed = TRUE)), pr$level)
    }
    list(known_pairs = known_pairs, cancer_roles = roles, oncokb_genes = okb_genes, oncokb_pairs = okb_pairs)
}

#' @name sc_rna_fusion_tier
#' @title sc_rna_fusion_tier
#' @description
#'
#' Tier of each RNA fusion, on the scale of the DNA driver tiers (OncoKB:
#' 1 actionable, 2 significant, 3 VUS). Only productive fusions (not
#' read-through, not 5'-5' / 3'-3') rank above 3.
#'   Tier 1 (known / actionable driver): an OncoKB fusion biomarker pair
#'     (e.g. FGFR3::TACC3); a fusion of an OncoKB fusion gene (NTRK1-3, ALK,
#'     ROS1, RET, FGFR1-3, BRAF, MET, PDGFRA/B, ...) with a partner, in frame
#'     or (3' partner) keeping its kinase domain; an Arriba known (Mitelman /
#'     literature) pair in frame, or DNA-matched in >= 2 cells without a known
#'     frame shift; an in-frame EGFR intragenic
#'     deletion (EGFRvIII-like); or the RNA call of a tier-1 DNA fusion.
#'   Tier 2 (cancer-gene fusion with evidence): a cancer gene (OncoKB cancer
#'     genes, COSMIC CGC) at a genic breakpoint and in frame, DNA-matched,
#'     keeping a kinase domain, or in >= min_fraction of the RNA cells; any
#'     other known pair; an in-frame DNA-matched fusion; or the RNA call of a
#'     tier-2 DNA fusion.
#'   Tiers 1-2 also need support beyond one cell (>= 2 carrier cells or a
#'   DNA-matched gene pair).
#'   Tier 3: everything else exported.
#'
#' @param fus fusions (gene1, gene2, site1, site2, type, reading_frame, known, n_cells, retained_domains, dna_gene, dna_tier)
#' @param refs sc_rna_fusion_references()
#' @param n_cells_rna RNA cells of the patient (for the cell fraction)
#' @param min_fraction cell fraction making a cancer-gene fusion tier 2
#' @return data.table(tier, tier_label, tier_reasons (list), cancer_genes (list), oncokb_level, productive, cell_fraction)
#' @export
#' @author Stanley Clarke
sc_rna_fusion_tier <- function(fus, refs = sc_rna_fusion_references(), n_cells_rna = NA, min_fraction = 0.1) {
    n <- nrow(fus)
    col <- function(x, d) if (is.null(fus[[x]])) rep(d, n) else fus[[x]]
    clean <- function(g) toupper(sub("\\(.*$", "", sub(",.*$", "", g)))
    g1 <- clean(fus$gene1); g2 <- clean(fus$gene2)
    genic1 <- col("site1", ".") != "intergenic"; genic2 <- col("site2", ".") != "intergenic"
    type <- col("type", ".")
    readthrough <- grepl("read-through", type)
    nonprod <- grepl("5'-5'|3'-3'", type)
    productive <- !readthrough & !nonprod
    frame <- col("reading_frame", ".")
    inframe <- frame %in% "in-frame"
    dna_gene <- col("dna_gene", NA_character_); dna_tier <- col("dna_tier", NA_integer_)
    dna <- !is.na(dna_gene)
    known <- col("known", FALSE) %in% TRUE
    dom <- strsplit(ifelse(is.na(col("retained_domains", ".")), ".", col("retained_domains", ".")), "|", fixed = TRUE)
    side <- function(k) vapply(dom, function(d) if (length(d) >= k) d[k] else ".", "")
    ## a protein kinase domain (Pfam "Protein kinase domain", "Protein tyrosine (and serine/threonine) kinase"),
    ## at least half of it kept; metabolic kinases (PI3/4-kinase, NDK, ...) do not count
    kinase <- function(d) vapply(regmatches(d, gregexpr("(Protein_kinase_domain|Protein_tyrosine[^,(]*kinase)[^,]*\\((\\d+)%\\)", d, ignore.case = TRUE)), function(m)
        length(m) > 0 && any(as.integer(sub(".*\\((\\d+)%\\)$", "\\1", m)) >= 50), TRUE)
    kin1 <- genic1 & kinase(side(1)); kin2 <- genic2 & kinase(side(2))
    frac <- if (is.finite(n_cells_rna) && n_cells_rna > 0) col("n_cells", 0) / n_cells_rna else rep(NA_real_, n)
    role1 <- ifelse(genic1, unname(refs$cancer_roles[g1]), NA_character_)
    role2 <- ifelse(genic2, unname(refs$cancer_roles[g2]), NA_character_)
    cancer <- !is.na(role1) | !is.na(role2)
    okb1 <- ifelse(genic1, unname(refs$oncokb_genes[g1]), NA_character_)
    okb2 <- ifelse(genic2, unname(refs$oncokb_genes[g2]), NA_character_)
    pair_lvl <- unname(refs$oncokb_pairs[paste(g1, g2, sep = "::")])
    ## frame not out-of-frame / stop-codon: in frame, or not determined (".": e.g. intronic / UTR breakpoints)
    frame_ok <- inframe | frame %in% c(".", "", NA)
    ## OncoKB fusion gene: in frame, or the 3' partner keeping its kinase domain without a known frame shift
    okb_partner <- g1 != g2 & ((!is.na(okb1) & inframe) | (!is.na(okb2) & (inframe | (kin2 & frame_ok))))
    egfrviii <- g1 == "EGFR" & g2 == "EGFR" & inframe & grepl("deletion", type)
    recurrent <- (frac >= min_fraction) %in% TRUE
    ## support beyond one cell: >= 2 carrier cells or the same gene pair in the DNA
    supported <- col("n_cells", 1L) >= 2 | dna
    ## a known pair: in frame, or DNA-matched in >= 2 cells without a known frame shift
    known1 <- known & (inframe | (dna & col("n_cells", 1L) >= 2 & frame_ok))
    t1 <- productive & supported & (!is.na(pair_lvl) | okb_partner | known1 | egfrviii | (dna & dna_tier %in% 1L))
    t2 <- productive & supported & !t1 & ((cancer & (inframe | dna | (kin2 & frame_ok) | recurrent)) | known | (dna & inframe) | (dna & dna_tier %in% 2L))
    tier <- ifelse(t1, 1L, ifelse(t2, 2L, 3L))
    labels <- c("known / actionable driver", "cancer-gene fusion", "other")
    reasons <- lapply(seq_len(n), function(i) {
        r <- character(0)
        if (readthrough[i]) r <- c(r, "read-through (adjacent genes)")
        if (nonprod[i]) r <- c(r, "non-productive orientation (5'-5' / 3'-3')")
        if (!is.na(pair_lvl[i])) r <- c(r, sprintf("OncoKB fusion %s::%s (level %s)", g1[i], g2[i], pair_lvl[i]))
        if (!is.na(okb1[i])) r <- c(r, sprintf("OncoKB fusion gene %s (level %s)", g1[i], okb1[i]))
        if (!is.na(okb2[i]) && g2[i] != g1[i]) r <- c(r, sprintf("OncoKB fusion gene %s (level %s)", g2[i], okb2[i]))
        if (egfrviii[i]) r <- c(r, "in-frame EGFR intragenic deletion (EGFRvIII-like)")
        if (known[i]) r <- c(r, "known fusion (Arriba / Mitelman)")
        if (!is.na(role1[i])) r <- c(r, sprintf("cancer gene %s (%s)", g1[i], role1[i]))
        if (!is.na(role2[i]) && g2[i] != g1[i]) r <- c(r, sprintf("cancer gene %s (%s)", g2[i], role2[i]))
        if (inframe[i]) r <- c(r, "in frame") else if (!frame[i] %in% c(".", "")) r <- c(r, frame[i])
        if (kin1[i]) r <- c(r, sprintf("%s kinase domain retained", g1[i]))
        if (kin2[i] && g2[i] != g1[i]) r <- c(r, sprintf("%s kinase domain retained", g2[i]))
        if (dna[i]) r <- c(r, if (is.na(dna_tier[i])) sprintf("DNA fusion %s", dna_gene[i]) else sprintf("DNA fusion %s (DNA tier %d)", dna_gene[i], dna_tier[i]))
        if (recurrent[i]) r <- c(r, sprintf("in %.0f%% of RNA cells", 100 * frac[i]))
        if (!supported[i]) r <- c(r, "single cell, no DNA support")
        r
    })
    cg <- lapply(seq_len(n), function(i) {
        out <- list()
        if (!is.na(role1[i])) out <- c(out, list(list(gene = g1[i], role = role1[i])))
        if (!is.na(role2[i]) && g2[i] != g1[i]) out <- c(out, list(list(gene = g2[i], role = role2[i])))
        out
    })
    lvl <- ifelse(!is.na(pair_lvl), pair_lvl, ifelse(!is.na(okb1), okb1, okb2))
    data.table::data.table(tier = tier, tier_label = labels[tier], tier_reasons = reasons, cancer_genes = cg,
                           oncokb_level = lvl, productive = productive, cell_fraction = round(frac, 4))
}

#' @name sc_rna_fusion_recurrence
#' @title sc_rna_fusion_recurrence
#' @description
#'
#' Cohort recurrence of the RNA fusions: adds to every fusion of every
#' data/<P>/rna/fusions.json the other patients with the same gene pair
#' (either order, first gene of intergenic fields), as `recurrence`
#' [{patient, n_cells}] (most cells first) and `n_patients`. A tier 1-2
#' fusion seen in `artefact_patients` or more patients that is neither known
#' nor DNA-matched drops to tier 3 (recurrent single-cell calls without DNA
#' support across unrelated tumours are mostly alignment / transcriptional
#' artefacts, e.g. BOLA2B::SMG1); the export's tier is kept as `tier_base`.
#' Run after the per-patient exports; rewrites the files in place (fields only added).
#'
#' @param data_dir data folder of the gOS dataset
#' @param artefact_patients patients from which an unsupported recurrent fusion is demoted
#' @return data.table(patient, pair, n_cells) (invisibly)
#' @export
#' @author Stanley Clarke
sc_rna_fusion_recurrence <- function(data_dir, artefact_patients = 3) {
    files <- Sys.glob(file.path(data_dir, "*", "rna", "fusions.json"))
    if (!length(files)) return(invisible(data.table::data.table()))
    pats <- basename(dirname(dirname(files)))
    js <- lapply(files, jsonlite::fromJSON, simplifyVector = FALSE)
    clean <- function(g) toupper(sub("\\(.*$", "", sub(",.*$", "", g)))
    key <- function(f) { a <- clean(f$gene1); b <- clean(f$gene2); paste(pmin(a, b), pmax(a, b), sep = "::") }
    tab <- data.table::rbindlist(lapply(seq_along(js), function(i) {
        fs <- js[[i]]$fusions
        if (!length(fs)) return(NULL)
        data.table::data.table(patient = pats[i], pair = vapply(fs, key, ""), n_cells = vapply(fs, function(f) as.integer(f$n_cells), 1L))
    }))
    if (!nrow(tab)) return(invisible(tab))
    tab <- tab[, .(n_cells = max(n_cells)), by = .(patient, pair)][order(-n_cells)]
    by_pair <- split(tab, tab$pair)
    for (i in seq_along(js)) {
        if (!length(js[[i]]$fusions)) next
        js[[i]]$fusions <- lapply(js[[i]]$fusions, function(f) {
            o <- by_pair[[key(f)]]
            o <- o[patient != pats[i]]
            f$recurrence <- lapply(seq_len(nrow(o)), function(k) list(patient = o$patient[k], n_cells = o$n_cells[k]))
            f$n_patients <- nrow(o) + 1L
            if (!is.null(f$tier)) {
                if (is.null(f$tier_base)) { f$tier_base <- f$tier; f$tier_reasons_base <- f$tier_reasons }
                f$tier <- f$tier_base
                f$tier_reasons <- f$tier_reasons_base
                f$tier_label <- c("known / actionable driver", "cancer-gene fusion", "other")[f$tier]
                if (f$tier_base < 3 && f$n_patients >= artefact_patients && !isTRUE(f$known) && is.null(f$dna_match)) {
                    f$tier <- 3L
                    f$tier_label <- "other"
                    f$tier_reasons <- c(f$tier_reasons, list(sprintf("in %d patients without DNA support (likely artefact)", f$n_patients)))
                }
            }
            f
        })
        jsonlite::write_json(js[[i]], files[i], auto_unbox = TRUE, null = "null", na = "null", digits = NA)
    }
    message(sprintf("sc_rna_fusion_recurrence: %d patients, %d gene pairs in >1 patient", length(files),
                    sum(tab[, data.table::uniqueN(patient), by = pair]$V1 > 1)))
    invisible(tab)
}

## ------------------------------------------------------------------ splicing

#' @name sc_read_junctions
#' @title sc_read_junctions
#' @description
#'
#' Per-cell splice junctions from regtools `junctions extract` BED files:
#' intron coordinates (1-based first and last intronic base, chromosome
#' without "chr") and read counts.
#'
#' @param cell_dirs named character rna_id -> folder holding junctions.bed
#' @param cores parallel reads
#' @return data.table(rna_id, chromosome, start, end, count)
#' @export
#' @author Stanley Clarke
sc_read_junctions <- function(cell_dirs, cores = 4) {
    files <- file.path(cell_dirs, "junctions.bed")
    data.table::rbindlist(parallel::mclapply(which(file.exists(files)), function(i) {
        ## blockSizes "77,33" must stay text (fread would read it as the decimal 77.33)
        x <- data.table::fread(files[i], header = FALSE, select = c(1:3, 5, 11), showProgress = FALSE,
                               colClasses = list(character = c(1, 11)), col.names = c("chr", "s", "e", "count", "blocks"))
        if (!nrow(x)) return(NULL)
        b <- data.table::tstrsplit(x$blocks, ",", fixed = TRUE, type.convert = TRUE)
        ## unstranded regtools may list one intron more than once: sum
        x[, .(rna_id = names(cell_dirs)[i], chromosome = sub("^chr", "", chr),
              start = as.integer(s + b[[1]] + 1L), end = as.integer(e - b[[2]]), count = as.integer(count))
          ][, .(count = sum(count)), by = .(rna_id, chromosome, start, end)]
    }, mc.cores = cores, mc.preschedule = FALSE))
}

## introns of the GTF transcripts: chromosome, start, end, strand, gene, transcript, exon numbers
sc_gtf_introns <- function(gtf) {
    ex <- if (data.table::is.data.table(gtf)) gtf else sc_gtf_exons(gtf)
    data.table::setorder(ex, transcript, start)
    ex[, .(chromosome = chromosome[-.N], start = end[-.N] + 1L, end = start[-1] - 1L, strand = strand[1], gene = gene[1],
           exon_left = exon[-.N], exon_right = exon[-1]), by = transcript][start <= end]
}

## exons of the GTF: chromosome (no "chr"), start, end, strand, gene, transcript (no version), exon number
sc_gtf_exons <- function(gtf) {
    ex <- data.table::as.data.table(rtracklayer::import(gtf, feature.type = "exon"))
    ex <- ex[, .(chromosome = sub("^chr", "", as.character(seqnames)), start, end, strand = as.character(strand),
                 gene = gene_name, transcript = sub("\\..*$", "", transcript_id), exon = as.integer(exon_number))]
    data.table::setorder(ex, transcript, start)
    ex[]
}

## collapsed gene models: the union of each gene's exons (all transcripts)
sc_merged_exons <- function(ex) {
    m <- unique(ex[, .(gene, chromosome, strand, start, end)])
    data.table::setorder(m, gene, chromosome, start, end)
    m[, grp := cumsum(c(TRUE, start[-1] > cummax(end)[-.N])), by = .(gene, chromosome)]
    m[, .(strand = strand[1], start = min(start), end = max(end)), by = .(gene, chromosome, grp)][, grp := NULL][]
}

## LeafCutter-style clusters: junctions of one chromosome/strand sharing a donor or acceptor
sc_junction_clusters <- function(j) {
    j <- data.table::copy(j)
    j[, jid := .I]
    parent <- seq_len(nrow(j))
    find <- function(x) { while (parent[x] != x) { parent[x] <<- parent[parent[x]]; x <- parent[x] }; x }
    for (key in list(c("chromosome", "strand", "start"), c("chromosome", "strand", "end"))) {
        grp <- j[, .(ids = list(jid)), by = key]$ids
        for (g in grp) if (length(g) > 1) for (k in g[-1]) { a <- find(g[1]); b <- find(k); if (a != b) parent[b] <- a }
    }
    j[, cluster := vapply(jid, find, 1L)]
    j[, n_in_cluster := .N, by = cluster]
    j[n_in_cluster > 1][, `:=`(jid = NULL, n_in_cluster = NULL)][]
}

#' @name sc_junction_types
#' @title sc_junction_types
#' @description
#'
#' Junction classes against the GTF: "annotated" (an intron of a transcript),
#' "exon_skip" (annotated donor and acceptor, unannotated pair, whole exons
#' of the gene in between), "novel_combination" (annotated ends, no exon in
#' between), "novel_donor" / "novel_acceptor" (one end annotated; donor and
#' acceptor by strand) or "novel" (neither end annotated); plus the number of
#' (collapsed) exons of the gene inside the intron.
#'
#' @param j data.table(chromosome, start, end, strand, gene, annotated)
#' @param introns GTF introns (chromosome, start, end)
#' @param merged collapsed exons (sc_merged_exons)
#' @return j with type and n_skipped
#' @export
sc_junction_types <- function(j, introns, merged) {
    j <- data.table::copy(j)
    s_known <- paste(j$chromosome, j$start) %in% paste(introns$chromosome, introns$start)
    e_known <- paste(j$chromosome, j$end) %in% paste(introns$chromosome, introns$end)
    donor <- ifelse(j$strand == "-", e_known, s_known)
    acceptor <- ifelse(j$strand == "-", s_known, e_known)
    j[, jid := .I]
    inside <- merged[j, on = .(gene, chromosome, start > start, end < end), .(jid = i.jid), nomatch = 0L, allow.cartesian = TRUE]
    j[, n_skipped := 0L]
    if (nrow(inside)) j[inside[, .N, by = jid], on = "jid", n_skipped := i.N]
    j[, type := data.table::fifelse(annotated, "annotated",
                data.table::fifelse(donor & acceptor, data.table::fifelse(n_skipped > 0, "exon_skip", "novel_combination"),
                data.table::fifelse(donor, "novel_acceptor", data.table::fifelse(acceptor, "novel_donor", "novel"))))]
    j[, jid := NULL][]
}

## canonical-transcript exon junctions of known GBM splice variants (alt vs reference junction)
SC_SPLICE_VARIANTS <- list(
    list(id = "EGFRvIII", gene = "EGFR", transcript = "ENST00000275493", alt = c(1, 8), ref = c(1, 2),
         description = "EGFR exon 1 → exon 8 (Δ exons 2–7, EGFRvIII)"),
    list(id = "EGFRvII", gene = "EGFR", transcript = "ENST00000275493", alt = c(13, 16), ref = c(13, 14),
         description = "EGFR exon 13 → exon 16 (Δ exons 14–15, EGFRvII)"),
    list(id = "EGFR_d25_27", gene = "EGFR", transcript = "ENST00000275493", alt = c(24, 28), ref = c(24, 25),
         description = "EGFR exon 24 → exon 28 (Δ exons 25–27, C-terminal deletion)"),
    list(id = "EGFR_d25_26", gene = "EGFR", transcript = "ENST00000275493", alt = c(24, 27), ref = c(24, 25),
         description = "EGFR exon 24 → exon 27 (Δ exons 25–26, C-terminal deletion)"),
    list(id = "MET_ex14_skip", gene = "MET", transcript = "ENST00000397752", alt = c(13, 15), ref = c(13, 14),
         description = "MET exon 13 → exon 15 (exon 14 skipping)"),
    list(id = "PDGFRA_d8_9", gene = "PDGFRA", transcript = "ENST00000257290", alt = c(7, 10), ref = c(7, 8),
         description = "PDGFRA exon 7 → exon 10 (Δ exons 8–9)"))

## genes whose intron clusters are always exported per patient (when the patient has reads there)
SC_SPLICE_GENES <- c("EGFR", "PDGFRA", "MET", "PTPRZ1", "CD44", "PKM", "MKNK2", "BCL2L1", "MDM4", "TP53", "PTEN", "NF1",
                     "CDK4", "MDM2", "SOX2", "OLIG2", "NCAM1", "BIN1", "FN1", "TNC", "VEGFA", "KLHDC4", "SRSF3", "PTBP1")

sc_variant_junctions <- function(gtf_exons) {
    rbindlist_safe <- function(x) data.table::rbindlist(Filter(Negate(is.null), x))
    rbindlist_safe(lapply(SC_SPLICE_VARIANTS, function(v) {
        e <- gtf_exons[sub("\\..*$", "", transcript) == v$transcript]
        if (!nrow(e)) return(NULL)
        jn <- function(a, b) {
            l <- e[exon == a]; r <- e[exon == b]
            if (!nrow(l) || !nrow(r)) return(c(NA, NA))
            if (l$strand[1] == "-") c(r$end[1] + 1L, l$start[1] - 1L) else c(l$end[1] + 1L, r$start[1] - 1L)
        }
        a <- jn(v$alt[1], v$alt[2]); r <- jn(v$ref[1], v$ref[2])
        data.table::data.table(id = v$id, gene = v$gene, description = v$description, chromosome = e$chromosome[1],
                               strand = e$strand[1], alt_start = a[1], alt_end = a[2], ref_start = r[1], ref_end = r[2])
    }))
}

## IGV slices of per-cell RNA reads at the known-variant junctions (both exon sides of the alt and
## the reference junction): rna/splice_reads/<rna_id>.bam. Returns rna_id -> relative path.
sc_splice_slices <- function(rna_ids, cell_dirs, variants, rna_dir, pad = 300, cores = 4,
                             samtools = Sys.getenv("SC_SAMTOOLS", "/gpfs/commons/groups/imielinski_lab/Software/miniforge3/envs/mskilab_ne1/bin/samtools")) {
    if (!length(rna_ids) || !nrow(variants)) return(list())
    out_dir <- file.path(rna_dir, "splice_reads")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    v <- variants[!is.na(alt_start) & !is.na(ref_start)]
    pos <- unique(data.table::data.table(chromosome = rep(v$chromosome, 4),
                                         pos = c(v$alt_start - 1L, v$alt_end + 1L, v$ref_start - 1L, v$ref_end + 1L)))
    regions <- pos[, paste0("chr", chromosome, ":", pmax(1L, pos - pad), "-", pos + pad)]
    ok <- parallel::mclapply(rna_ids, function(r) {
        bam <- file.path(cell_dirs[[r]], "Aligned.sorted.bam")
        if (!file.exists(paste0(bam, ".bai"))) return(FALSE)
        out <- file.path(out_dir, paste0(r, ".bam"))
        rc <- system2(samtools, c("view", "-b", "-M", "-o", shQuote(out), shQuote(bam), regions), stdout = FALSE, stderr = FALSE)
        if (rc != 0) return(FALSE)
        system2(samtools, c("index", shQuote(out)), stdout = FALSE, stderr = FALSE) == 0 && file.exists(paste0(out, ".bai"))
    }, mc.cores = cores, mc.preschedule = FALSE)
    done <- rna_ids[unlist(ok) %in% TRUE]
    stats::setNames(as.list(paste0("rna/splice_reads/", done, ".bam")), done)
}

## per-cell overdispersion of a patient's clusters: for each junction, the variance of per-cell PSI
## (cells with >= min_total reads in the cluster) minus the binomial variance expected from read
## depth alone; a cluster scores its most overdispersed junction
sc_cluster_dispersion <- function(x, min_total = 5, min_cells = 10) {
    xc <- x[cl_total >= min_total]
    if (!nrow(xc)) return(data.table::data.table(cluster_id = character(0), n_cells = integer(0), excess_sd = numeric(0)))
    cc <- unique(xc[, .(rna_id, cluster_id, cl_total)])[, .(n_cells = .N, N = sum(cl_total), minv = mean(1 / cl_total)), by = cluster_id]
    jj <- xc[, .(s1 = sum(count / cl_total), s2 = sum((count / cl_total)^2), k = sum(count)), by = .(cluster_id, start, end)]
    jj <- merge(jj, cc, by = "cluster_id")
    jj[, `:=`(m = s1 / n_cells, pbar = k / N)]
    jj[, excess := (s2 / n_cells - m^2) - pbar * (1 - pbar) * minv]
    jj[n_cells >= min_cells, .(n_cells = n_cells[1], excess_sd = sqrt(max(0, max(excess)))), by = cluster_id]
}

## clusters whose per-cell junction usage differs between groups (DNA clones): per junction a
## Kruskal-Wallis test of per-cell PSI (cells with >= min_total reads, groups with >= min_cells
## such cells); cluster p = smallest junction p x junctions (Bonferroni), q = BH over clusters
sc_cluster_group_test <- function(x, group_of, min_total = 3, min_cells = 5) {
    xc <- x[cl_total >= min_total]
    xc[, group := group_of[rna_id]]
    xc <- xc[!is.na(group)]
    if (!nrow(xc)) return(data.table::data.table(cluster_id = character(0), group_p = numeric(0), group_q = numeric(0)))
    cells <- unique(xc[, .(rna_id, cluster_id, group)])
    jn <- unique(xc[, .(cluster_id, start, end)])
    ## zero-filled per-cell PSI of every junction of the cluster
    full <- merge(cells, jn, by = "cluster_id", allow.cartesian = TRUE)
    full <- merge(full, xc[, .(rna_id, cluster_id, start, end, psi = count / cl_total)], by = c("rna_id", "cluster_id", "start", "end"), all.x = TRUE)
    full[is.na(psi), psi := 0]
    full[, n_group := data.table::uniqueN(rna_id), by = .(cluster_id, group)]
    full <- full[n_group >= min_cells]
    res <- full[, if (data.table::uniqueN(group) >= 2) .(p = suppressWarnings(stats::kruskal.test(psi, factor(group))$p.value)), by = .(cluster_id, start, end)]
    if (!nrow(res)) return(data.table::data.table(cluster_id = character(0), group_p = numeric(0), group_q = numeric(0)))
    res <- res[!is.na(p), .(group_p = min(1, min(p) * .N)), by = cluster_id]
    res[, group_q := stats::p.adjust(group_p, "BH")][]
}

#' @name sc_export_rna_splicing
#' @title sc_export_rna_splicing
#' @description
#'
#' Splice junction usage of the single-cell RNA for gOS. Junctions of every
#' cell (regtools) are pooled per patient; junctions with at least `min_reads`
#' reads in `min_cells` cells (cohort) are grouped into LeafCutter-style
#' clusters (shared donor or acceptor, same strand). Strand and gene come from
#' the GTF introns (annotated) or, for novel junctions, the GTF gene they fall
#' in. Every junction is typed against the GTF (sc_junction_types) and every
#' cluster carries the collapsed exons of its gene around it (for sashimi
#' plots). Writes:
#'  - data/_cohort/rna/splicing.json: clusters whose junction usage differs
#'    between patients (chi-square on patient x junction pseudobulk counts,
#'    BH q, max delta-PSI between patients with >= `min_patient_reads`) and
#'    the known variants per patient,
#'  - data/<P>/rna/splicing.json per patient: per-cell counts of the clusters
#'    most overdispersed between the patient's cells (sc_cluster_dispersion) or
#'    differing between its DNA clones (sc_cluster_group_test, `clones`),
#'    the cohort clusters and those of SC_SPLICE_GENES, and per-cell alt /
#'    reference reads of known GBM splice variants (EGFRvIII, EGFRvII, EGFR
#'    C-terminal deletions, MET exon 14 skipping, PDGFRA delta 8-9), with IGV
#'    slices of the variant loci for carrier cells (`slices`).
#'
#' @param cell_dirs list patient -> named character (rna_id -> per-cell folder with junctions.bed)
#' @param data_dir data folder of the gOS dataset
#' @param gtf GTF of the alignment reference
#' @param cell_maps list patient -> named character rna_id -> gOS cell id; RNA ids of the
#'   junction folders (cell_dirs names) are matched to these case-insensitively and renamed to them
#' @param clones list patient -> named character rna_id -> DNA clone (clone-differential clusters)
#' @param min_reads,min_cells cohort junction filter
#' @param min_patient_reads patient cluster reads for delta-PSI
#' @param n_cohort,n_patient,n_novel clusters exported (between patients; per patient by dispersion / clone test; with novel junctions)
#' @param slices write per-cell RNA slices of the variant loci (carriers + `n_ref_slices` reference cells per variant)
#' @param cores parallel reads
#' @return list(cohort = cohort clusters, patients = per-patient cluster ids) invisibly
#' @export
#' @author Stanley Clarke
sc_export_rna_splicing <- function(cell_dirs, data_dir, gtf, cell_maps = list(), clones = list(), min_reads = 30, min_cells = 10,
                                   min_patient_reads = 30, n_cohort = 300, n_patient = 100, n_novel = 100, slices = TRUE, n_ref_slices = 6,
                                   cores = 8) {
    message("sc_export_rna_splicing: GTF introns")
    ex <- sc_gtf_exons(gtf)
    introns <- sc_gtf_introns(ex)
    ann <- unique(introns[, .(chromosome, start, end, strand, gene)], by = c("chromosome", "start", "end"))
    merged <- sc_merged_exons(ex)
    genes <- ex[, .(start = min(start), end = max(end), strand = strand[1], chromosome = chromosome[1]), by = gene]
    variants <- sc_variant_junctions(ex)

    ## per patient: per-cell junction counts
    per <- list()
    for (p in names(cell_dirs)) {
        ## folder names (sample sheet) may differ in case from the Seurat / gOS RNA ids: use the latter
        canon <- names(cell_maps[[p]])
        if (length(canon)) {
            m <- match(tolower(names(cell_dirs[[p]])), tolower(canon))
            names(cell_dirs[[p]])[!is.na(m)] <- canon[m[!is.na(m)]]
        }
        message("sc_export_rna_splicing: reading junctions of ", p, " (", length(cell_dirs[[p]]), " cells)")
        per[[p]] <- sc_read_junctions(cell_dirs[[p]], cores = cores)
        data.table::setkey(per[[p]], chromosome, start, end)
    }
    allj <- data.table::rbindlist(lapply(names(per), function(p) per[[p]][, .(reads = sum(count), cells = .N), by = .(chromosome, start, end)][, patient := p]))
    tot <- allj[, .(reads = sum(reads), cells = sum(cells)), by = .(chromosome, start, end)][reads >= min_reads & cells >= min_cells]
    ## strand / gene: annotated intron, else the gene the junction lies in (unique strand), else drop
    tot <- merge(tot, ann, by = c("chromosome", "start", "end"), all.x = TRUE)
    tot[, annotated := !is.na(strand)]
    nov <- tot[is.na(strand)]
    if (nrow(nov)) {
        g <- genes[, .(gchr = chromosome, gs = start, ge = end, gstrand = strand, gname = gene)]
        hit <- g[nov, on = .(gchr = chromosome, gs <= start, ge >= end), .(chromosome = i.chromosome, start = i.start, end = i.end, gstrand, gname),
                 allow.cartesian = TRUE, nomatch = 0L]
        hit <- hit[, .(strand = if (data.table::uniqueN(gstrand) == 1) gstrand[1] else NA_character_, gene = gname[1]), by = .(chromosome, start, end)]
        tot[hit, on = .(chromosome, start, end), `:=`(strand = i.strand, gene = i.gene)]
    }
    tot <- tot[!is.na(strand)]
    cl <- sc_junction_clusters(tot[, .(chromosome, strand, start, end, gene, annotated)])
    cl[, gene := { g <- gene[!is.na(gene)]; if (length(g)) names(sort(table(g), decreasing = TRUE))[1] else NA_character_ }, by = cluster]
    cl[, cluster_id := paste0("clu_", chromosome, "_", min(start), "_", max(end), "_", strand), by = cluster]
    cl <- sc_junction_types(cl, introns, merged)
    data.table::setorder(cl, cluster_id, start, end)
    message("sc_export_rna_splicing: ", data.table::uniqueN(cl$cluster_id), " clusters from ", nrow(cl), " junctions (",
            paste(names(table(cl$type)), table(cl$type), collapse = ", "), ")")
    cl_by <- split(cl, cl$cluster_id)

    ## JSON parts shared by the cohort and patient files
    junctions_of <- function(jj) lapply(seq_len(nrow(jj)), function(k) list(start = jj$start[k], end = jj$end[k], annotated = jj$annotated[k],
                                                                             type = jj$type[k], n_skipped = jj$n_skipped[k]))
    exons_of <- function(jj) {
        e <- merged[gene == jj$gene[1] & chromosome == jj$chromosome[1] & end >= min(jj$start) - 1L & start <= max(jj$end) + 1L]
        lapply(seq_len(nrow(e)), function(k) c(e$start[k], e$end[k]))
    }
    cluster_head <- function(id) {
        jj <- cl_by[[id]]
        list(id = id, gene = jj$gene[1], chromosome = jj$chromosome[1], strand = jj$strand[1], start = min(jj$start), end = max(jj$end),
             junctions = junctions_of(jj), exons = exons_of(jj))
    }

    ## cohort: patient x junction pseudobulk per cluster
    pb <- merge(allj, cl[, .(chromosome, start, end, cluster_id)], by = c("chromosome", "start", "end"))
    patients <- names(per)
    stats <- pb[, {
        m <- data.table::dcast(.SD, patient ~ start + end, value.var = "reads", fill = 0, fun.aggregate = sum)
        mat <- as.matrix(m[, -1, with = FALSE]); rownames(mat) <- m$patient
        mat <- mat[rowSums(mat) >= min_patient_reads, , drop = FALSE]
        if (nrow(mat) < 2 || ncol(mat) < 2) list(p = NA_real_, max_dpsi = NA_real_, n_patients = nrow(mat)) else {
            psi <- mat / rowSums(mat)
            list(p = suppressWarnings(stats::chisq.test(mat)$p.value), max_dpsi = max(apply(psi, 2, function(v) diff(range(v)))),
                 n_patients = nrow(mat))
        }
    }, by = cluster_id, .SDcols = c("patient", "start", "end", "reads")]
    stats <- stats[!is.na(p)]
    stats[, q := stats::p.adjust(p, "BH")]
    top <- stats[q < 0.05][order(-max_dpsi, q)][seq_len(min(.N, n_cohort))]
    ## cells with any read in a cluster, per patient
    cl_key <- cl[, .(chromosome, start, end, cluster_id)]
    cl_cells <- data.table::rbindlist(lapply(patients, function(p)
        merge(per[[p]], cl_key, by = c("chromosome", "start", "end"))[, .(n_cells = data.table::uniqueN(rna_id)), by = cluster_id][, patient := p]))
    cohort_clusters <- lapply(top$cluster_id, function(id) {
        jj <- cl_by[[id]]
        u <- pb[cluster_id == id]
        s <- top[cluster_id == id]
        c(cluster_head(id),
          list(p = s$p, q = s$q, max_dpsi = s$max_dpsi, n_patients = s$n_patients,
               usage = stats::setNames(lapply(patients, function(p) {
                   cnt <- vapply(seq_len(nrow(jj)), function(k) as.numeric(sum(u[patient == p & start == jj$start[k] & end == jj$end[k]]$reads)), 0)
                   nc <- cl_cells[cluster_id == id & patient == p]$n_cells
                   list(counts = cnt, total = sum(cnt), n_cells = if (length(nc)) nc else 0L)
               }), patients)))
    })

    ## known variants per patient (per-cell alt / ref reads)
    variant_cells <- function(p, v) {
        alt <- per[[p]][.(v$chromosome, v$alt_start, v$alt_end), .(rna_id, alt = count), nomatch = 0L]
        ref <- per[[p]][.(v$chromosome, v$ref_start, v$ref_end), .(rna_id, ref = count), nomatch = 0L]
        ar <- merge(alt, ref, by = "rna_id", all = TRUE)
        ar[is.na(alt), alt := 0L][is.na(ref), ref := 0L][]
    }
    var_tab <- list()
    for (p in patients) var_tab[[p]] <- lapply(seq_len(nrow(variants)), function(i) variant_cells(p, variants[i]))

    cdir <- file.path(data_dir, "_cohort", "rna")
    dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
    cohort_variants <- lapply(seq_len(nrow(variants)), function(i) {
        v <- variants[i]
        list(id = v$id, gene = v$gene, description = v$description,
             alt_junction = list(chromosome = v$chromosome, start = v$alt_start, end = v$alt_end),
             ref_junction = list(chromosome = v$chromosome, start = v$ref_start, end = v$ref_end),
             usage = stats::setNames(lapply(patients, function(p) {
                 ar <- var_tab[[p]][[i]]
                 list(alt = sum(ar$alt), ref = sum(ar$ref), n_cells_alt = sum(ar$alt > 0), n_cells = sum(ar$alt + ar$ref > 0))
             }), patients))
    })
    jsonlite::write_json(list(format = "gos-sc-splicing-cohort/1", patients = patients, n_clusters_tested = nrow(stats),
                              n_cells_rna = stats::setNames(lapply(patients, function(p) data.table::uniqueN(per[[p]]$rna_id)), patients),
                              variants = cohort_variants, clusters = cohort_clusters),
                         file.path(cdir, "splicing.json"), auto_unbox = TRUE, digits = NA, null = "null", na = "null")
    message("sc_export_rna_splicing: cohort: ", length(cohort_clusters), " clusters of ", nrow(stats), " tested")

    ## per patient: clusters most overdispersed between cells (+ cohort clusters + GBM genes), per-cell counts; known variants
    pat_ids <- list()
    for (p in patients) {
        x <- merge(per[[p]], cl_key, by = c("chromosome", "start", "end"))
        x[, cl_total := sum(count), by = .(rna_id, cluster_id)]
        ## overdispersion only for clusters covered in enough cells (few-cell clusters are all bursting noise)
        disp <- sc_cluster_dispersion(x, min_cells = max(20L, round(0.08 * data.table::uniqueN(per[[p]]$rna_id))))
        gbm <- unique(cl[gene %in% SC_SPLICE_GENES]$cluster_id)
        gbm <- intersect(gbm, x[, .(r = sum(count)), by = cluster_id][r >= 10]$cluster_id)
        ## clone-differential clusters (DNA clones of the linked cells; normal cells left out)
        cg <- clones[[p]]
        cg <- cg[!is.na(cg) & !grepl("^normal$", cg, ignore.case = TRUE)]
        gt <- if (length(unique(cg)) >= 2) sc_cluster_group_test(x, cg) else data.table::data.table(cluster_id = character(0), group_p = numeric(0), group_q = numeric(0))
        clone_ids <- gt[group_q < 0.1][order(group_p)][seq_len(min(.N, n_patient))]$cluster_id
        ## clusters with a well-used unannotated junction in this patient (exon skips, novel sites):
        ## >= 10 reads in >= 3 cells and >= 10% of the cluster's pooled reads
        ju <- x[, .(r = sum(count), nc = data.table::uniqueN(rna_id)), by = .(cluster_id, start, end)]
        ju[, psi := r / sum(r), by = cluster_id]
        ju <- merge(ju, cl[, .(cluster_id, start, end, type)], by = c("cluster_id", "start", "end"))
        novel_ids <- unique(ju[type != "annotated" & r >= 10 & nc >= 3 & psi >= 0.1][order(-nc, -r)]$cluster_id)
        novel_ids <- novel_ids[seq_len(min(length(novel_ids), n_novel))]
        disp_ids <- disp[order(-excess_sd)][seq_len(min(.N, n_patient))]$cluster_id
        cohort_ids <- intersect(top$cluster_id, unique(x$cluster_id))
        ids <- unique(c(cohort_ids, clone_ids, novel_ids, disp_ids, gbm))
        pat_ids[[p]] <- ids
        xs <- split(x[cluster_id %in% ids], by = "cluster_id")
        clusters <- lapply(ids, function(id) {
            jj <- cl_by[[id]]
            xx <- xs[[id]]
            w <- data.table::dcast(xx, rna_id ~ start + end, value.var = "count", fill = 0, fun.aggregate = sum)
            keyj <- paste(jj$start, jj$end, sep = "_")
            for (k in setdiff(keyj, names(w))) w[, (k) := 0L]
            m <- as.matrix(w[, keyj, with = FALSE])
            d <- disp[cluster_id == id]
            g <- gt[cluster_id == id]
            c(cluster_head(id),
              list(excess_sd = if (nrow(d)) d$excess_sd else NULL, clone_p = if (nrow(g)) g$group_p else NULL, clone_q = if (nrow(g)) g$group_q else NULL,
                   n_cells = nrow(w), cohort = id %in% top$cluster_id, gbm_gene = id %in% gbm,
                   selected_by = I(c("cohort", "clone", "novel", "dispersion", "gene")[c(id %in% cohort_ids, id %in% clone_ids, id %in% novel_ids, id %in% disp_ids, id %in% gbm)]),
                   cells = stats::setNames(lapply(seq_len(nrow(m)), function(i) unname(as.integer(m[i, ]))), w$rna_id)))
        })
        ## per-patient file only for patients already in the gOS dataset (cohort file covers all)
        if (!dir.exists(file.path(data_dir, p))) next
        pdir <- file.path(data_dir, p, "rna")
        dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
        bams <- list()
        if (slices && nrow(variants)) {
            pick <- unique(unlist(lapply(var_tab[[p]], function(ar)
                c(ar[alt > 0]$rna_id, ar[alt == 0 & ref > 0][order(-ref)][seq_len(min(.N, n_ref_slices))]$rna_id))))
            pick <- pick[!is.na(pick) & pick %in% names(cell_dirs[[p]])]
            bams <- sc_splice_slices(pick, cell_dirs[[p]], variants, pdir, cores = cores)
            message("sc_export_rna_splicing: ", p, ": variant-locus RNA slices for ", length(bams), " of ", length(pick), " cells")
        }
        vv <- lapply(seq_len(nrow(variants)), function(i) {
            v <- variants[i]
            ar <- var_tab[[p]][[i]]
            list(id = v$id, gene = v$gene, description = v$description, strand = v$strand,
                 alt_junction = list(chromosome = v$chromosome, start = v$alt_start, end = v$alt_end),
                 ref_junction = list(chromosome = v$chromosome, start = v$ref_start, end = v$ref_end),
                 n_cells_alt = sum(ar$alt > 0),
                 cells = stats::setNames(lapply(seq_len(nrow(ar)), function(k) c(ar$alt[k], ar$ref[k])), ar$rna_id))
        })
        cm <- cell_maps[[p]]
        empty_named <- structure(list(), names = character(0))
        jsonlite::write_json(list(format = "gos-sc-splicing/1", patient = p, n_cells = data.table::uniqueN(per[[p]]$rna_id),
                                  variants = vv, clusters = clusters,
                                  splice_reads = if (length(bams)) bams else empty_named,
                                  cell_map = if (length(cm)) as.list(cm) else empty_named),
                             file.path(pdir, "splicing.json"), auto_unbox = TRUE, digits = NA, null = "null", na = "null")
        message("sc_export_rna_splicing: ", p, ": ", length(ids), " clusters (", length(clone_ids), " clone-differential of ", nrow(gt), " tested, ",
                length(novel_ids), " with novel junctions, ",
                length(gbm), " GBM-gene), variants: ",
                paste(vapply(vv, function(v) paste0(v$id, " ", v$n_cells_alt), ""), collapse = ", "))
    }
    invisible(list(cohort = top, patients = pat_ids))
}
