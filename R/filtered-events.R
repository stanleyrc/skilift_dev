#' @title process_gencode
#' @description
#'
#' Helper script to process gencode parameter
#'
#' @param gencode path to gencode file. Gencode file must be either rds or some format accepted by rtracklayer::import (e.g. GTF) with just a single entry for each gene (so gencode entries for each gene are collapse to a single range). The input could be .gtf or .rds with GRanges object, or a GRanges object i.e. resulting from importing the (appropriate) GENCODE .gtf via rtracklayer, note: this input is only used in CNA to gene mapping.
#' @return gencode_gr GRanges
#' @author Marcin Imielinski
process_gencode <- function(gencode = NULL) {
  if (is.null(gencode)) {
    stop("gencode file must be provided")
  } else if (is.character(gencode)) {
    if (grepl(".rds$", gencode)) {
      gencode <- readRDS(gencode)
    } else {
      gencode <- rtracklayer::import(gencode)
    }
  }
  return(gencode)
}

#' @title process_cytoband
#' @description
#'
#' Helper script to process cytoband
#'
#' @param cytoband path to cytoband file provided by UCSC
#' @return GenomicRanges::GRanges object
#' @author Kevin Hadi
process_cytoband <- function(cyto = NULL, coarse = FALSE) {
  if (is.null(cyto)) {
    cyto <- system.file("extdata", "data", "cytoband.rds", package = "Skilift")
    warning("Using default cytoband, hg19")
    cyto <- readRDS(cyto)
    return(cyto)
  } else if (is.character(cyto)) {
    if (grepl(".rds$", cyto)) {
      cyto <- readRDS(cyto)
      return(cyto)
    } else {
      cyto <- data.table::fread(cyto)
    }
  }

  names(cyto) <- c("seqnames", "start", "end", "band", "stain")
  isZeroStart <- any(cyto[, length(intersect(start, end)), by = seqnames]$V1 > 0) || any(cyto$start == 0)
  if (isZeroStart) cyto$start <- cyto$start + 1
  cyto <- gUtils::dt2gr(cyto)
  GenomeInfoDb::seqlevelsStyle(cyto) <- "NCBI"
  cyto$chrom_name <- as.character(seqnames(cyto))
  pasteband <- cyto$band
  if (coarse) {
    pasteband <- gsub("\\.[0-9]+", "", cyto$band)
  }
  # cyto$rough_band = gsub("\\.[0-9]+", "", cyto$band)
  cyto$chromband <- paste(
    cyto$chrom_name,
    pasteband,
    sep = ""
  )
  return(cyto)
}



#' @title collect_gene_fusions
#' @description
#' Collects gene fusion data from a specified file and processes it.
#'
#' @param fusions Path to the fusions.rds file.
#' @param pge GRanges object with gencode annotations.
#' @param verbose Logical flag to indicate if messages should be printed.
#' @return A data.table containing processed gene fusion information.
collect_gene_fusions <- function(fusions, pge, verbose = TRUE) {
  if (is.null(fusions) || !file.exists(fusions)) {
    if (verbose) message("Fusions file is missing or does not exist.")
    return(data.table(type = NA, source = "fusions"))
  }

  if (verbose) message("pulling fusions")
  fus <- readRDS(fusions)$meta

  if (nrow(fus) == 0) {
    if (verbose) message("No fusions found in the file.")
    return(data.table(type = NA, source = "fusions"))
  }

  non_silent_fusions <- fus[silent == FALSE, ]
  unique_fusions <- non_silent_fusions[!duplicated(genes), ]
  unique_fusions[, vartype := ifelse(in.frame == TRUE, "fusion", "outframe_fusion")]

  split_genes <- function(genes) unlist(strsplit(genes, ","))
  gene_lengths <- function(genes) sapply(strsplit(genes, ","), length)

  fus <- unique_fusions[, .(
    gene = split_genes(genes),
    vartype = rep(vartype, gene_lengths(genes)),
    fusion_genes = rep(genes, gene_lengths(genes))
  )][, `:=`(track = "variants", type = vartype, source = "fusions")]

  get_gene_coords <- function(genes) {
    coords <- lapply(genes, function(gene) {
      gene_ranges <- pge[mcols(pge)$gene_name == gene]
      paste0(seqnames(gene_ranges), ":", start(gene_ranges), "-", end(gene_ranges))
    })
    paste(unlist(coords), collapse = ",")
  }

  fus[, fusion_gene_coords := unlist(lapply(strsplit(fusion_genes, ","), get_gene_coords))]

  return(fus)
}
#' @description
#' Collects complex events from a specified file and processes them.
#'
#' @param complex Path to the complex.rds file.
#' @param verbose Logical flag to indicate if messages should be printed.
#' @return A data.table containing processed complex event information.
collect_complex_events <- function(complex, verbose = TRUE) {
  if (is.null(complex) || !file.exists(complex)) {
    if (verbose) message("Complex events file is missing or does not exist.")
    return(data.table(type = NA, source = "complex"))
  }

  if (verbose) message("pulling complex events")
  sv <- readRDS(complex)$meta$events

  if (nrow(sv) == 0) {
    if (verbose) message("No complex events found in the file.")
    return(data.table(type = NA, source = "complex"))
  }

  sv_summary <- sv[, .(value = .N), by = type]
  simple_sv_types <- c("del", "dup", "invdup", "tra", "inv")
  sv_summary[, track := ifelse(type %in% simple_sv_types, "simple sv", "complex sv")]
  sv_summary[, source := "complex"]

  return(sv_summary)
}

check_GRanges_compatibility = function (gr1, gr2, name1 = "first", name2 = "second") 
{
    non_overlapping_seqnames1 = setdiff(seqlevels(gr1), seqlevels(gr2))
    non_overlapping_seqnames2 = setdiff(seqlevels(gr2), seqlevels(gr1))
    overlap = intersect(seqlevels(gr1), seqlevels(gr2))
    message("The following seqnames are only in the ", name1, 
        " GRanges, but not in the ", name2, " GRanges: ", paste(non_overlapping_seqnames1, 
            collapse = ", "))
    message("The following seqnames are only in the ", name2, 
        " GRanges, but not in the ", name1, " GRanges: ", paste(non_overlapping_seqnames2, 
            collapse = ", "))
    message("The follosing seqnames are in both GRanges objects: ", 
        paste(overlap, collapse = ", "))
    if (length(non_overlapping_seqnames1) > 0 | length(non_overlapping_seqnames2) > 
        0) {
        return(FALSE)
    }
    return(TRUE)
}

get_gene_copy_numbers <- function(
  gg, 
  gene_ranges, 
  nseg = NULL, 
  gene_id_col = "gene_id", 
  simplify_seqnames = FALSE, 
  mfields = c("gene_name", "source", "gene_id", "gene_type", "level", "hgnc_id", "havana_gene"), 
  output_type = "data.table",
  min_width = 1e3
) {
    if (is.character(gg)) {
        gg = readRDS(gg)
    }
    if (!inherits(gene_ranges, "GRanges")) {
        gene_ranges = rtracklayer::import(gene_ranges)
    }
    if (!(output_type %in% c("GRanges", "data.table"))) {
        stop("Invalid output_type: ", output_type, ". outputtype must be either \"GRanges\" or \"data.table\".")
    }
    ngr = gg$nodes$gr
    if (simplify_seqnames) {
        ngr = gr.sub(ngr)
        gene_ranges = gr.sub(gene_ranges)
    }
    GRanges_are_compatible = check_GRanges_compatibility(ngr, 
        gene_ranges, "gGraph segments", "genes")
    if (!is.null(nseg)) {
        ngr = ngr %$% nseg[, c("ncn")]
    }
    else {
        message("No normal copy number segmentation was provided so assuming CN = 2 for all seqnames.")
        ngr$ncn = 2
    }
    ndt = gr2dt(ngr)
    seq_widths = as.numeric(width(ngr))
    normal_ploidy = round(sum(seq_widths * ngr$ncn, na.rm = T)/sum(seq_widths, 
        na.rm = T))
    ndt[, `:=`(normalized_cn, cn * normal_ploidy/(gg$meta$ploidy * 
                                                    ncn))]
    gene_ranges$gene_width = width(gene_ranges)
    gene_cn_segments = dt2gr(ndt, seqlengths = seqlengths(gg)) %*% 
      gene_ranges %>% gr2dt
    
    # gene_cn_stats = gene_cn_segments[, .(
    #   max_normalized_cn = max(normalized_cn, na.rm = TRUE),
    #   max_cn = max(cn, na.rm = TRUE),
    #   min_normalized_cn = min(normalized_cn, na.rm = TRUE),
    #   min_cn = min(cn, na.rm = TRUE),
    #   avg_normalized_cn = sum(normalized_cn * width, na.rm = TRUE) / sum(width),
    #   avg_cn = sum(cn * width, na.rm = TRUE) / sum(width),
    #   total_node_width = sum(width, na.rm = TRUE),
    #   number_of_cn_segments = .N,
    #   ncn = ncn[1],
    #   list_of_segs = list(.SD[, list(gene_name = gene_name[1], normalized_cn, cn, width)])
    #   ## gene_width = gene_width[1]
    # ), by = gene_name]

    gene_cn_table = gene_cn_segments[, `:=`(
      max_normalized_cn = max(normalized_cn, na.rm = TRUE),
      max_cn = max(cn, na.rm = TRUE),
      min_normalized_cn = min(normalized_cn, na.rm = TRUE),
      min_cn = min(cn, na.rm = TRUE),
      avg_normalized_cn = sum(normalized_cn * width, na.rm = TRUE) / sum(width),
      avg_cn = sum(cn * width, na.rm = TRUE) / sum(width),
      # total_node_width = sum(width, na.rm = TRUE),
      number_of_cn_segments = .N,
      ncn = ncn[1],
      # list_of_segs = list(.SD[, list(gene_name = gene_name[1], normalized_cn, cn, width)])
      gene_width = gene_width[1]
    ), by = gene_name]
    
    # gene_cn_table = data.table::merge.data.table(
    #   gr2dt(gene_ranges),
    #   gene_cn_stats,
    #   by = "gene_name",
    #   suffixes = c("", "__DUPED")
    # )
    
    ## split_genes = gene_cn_segments[duplicated(get(gene_id_col)), 
    ##     get(gene_id_col)]
    ## gene_cn_non_split_genes = gene_cn_segments[!(get(gene_id_col) %in% 
    ##     split_genes)]
    ## gene_cn_non_split_genes[, `:=`(max_normalized_cn = normalized_cn, 
    ##     min_normalized_cn = normalized_cn, max_cn = cn, min_cn = cn, 
    ##     number_of_cn_segments = 1, cn = NULL, normalized_cn = NULL)]
    ## gene_cn_split_genes_min = gene_cn_segments[get(gene_id_col) %in% 
    ##     split_genes, .SD[which.min(cn)], by = gene_id_col]
    ## gene_cn_split_genes_min[, `:=`(min_normalized_cn = normalized_cn, 
    ##     min_cn = cn, cn = NULL, normalized_cn = NULL)]
    ## gene_cn_split_genes_max = gene_cn_segments[get(gene_id_col) %in% 
    ##     split_genes, .SD[which.max(cn)], by = gene_id_col][, 
    ##     .(get(gene_id_col), max_normalized_cn = normalized_cn, 
    ##         max_cn = cn)]
    ## setnames(gene_cn_split_genes_max, "V1", gene_id_col)
    ## number_of_segments_per_split_gene = gene_cn_segments[get(gene_id_col) %in% 
    ##     split_genes, .(number_of_cn_segments = .N), by = gene_id_col]
    ## gene_cn_split_genes = merge(gene_cn_split_genes_min, gene_cn_split_genes_max, 
    ##     by = gene_id_col)
    ## gene_cn_split_genes = merge(gene_cn_split_genes, number_of_segments_per_split_gene, 
    ##     by = gene_id_col)
    ## gene_cn_table = rbind(gene_cn_split_genes, gene_cn_non_split_genes)
    if (output_type == "data.table") {
        return(gene_cn_table)
    }
    return(dt2gr(gene_cn_table, seqlengths = seqlengths(gene_ranges)))
}

#' Get gene amplifications and deletions from jabba object
#'
#' Copied from github::mskilab-org/skitools.
#' This function takes in a jabba object and returns
#' amplifications and deletions
#'
#' @param jab character path to jabba file or gGraph object
get_gene_ampdels_from_jabba <- function(jab, pge, amp.thresh = 4, del.thresh = 0.5, nseg = NULL) {
  gg <- jab
  if (!inherits(gg, "gGraph")) {
    gg <- gG(jabba = jab)
  }
  gene_CN <- Skilift:::get_gene_copy_numbers(gg, gene_ranges = pge, nseg = nseg)
  gene_CN[, `:=`(type, NA_character_)]

  gene_CN[min_normalized_cn >= amp.thresh, `:=`(type, "amp")]
  gene_CN[min_cn > 1 & min_normalized_cn < del.thresh, `:=`(
    type,
    "del"
  )]
  gene_CN[min_cn == 1 & min_cn < ncn, `:=`(type, "hetdel")]
  gene_CN[min_cn == 0, `:=`(type, "homdel")]

  # scna_result = gene_CN[!is.na(type)]
  scna_result = (
      gene_CN[, 
      .(
        max_normalized_cn = max_normalized_cn[1],
        max_cn = max_cn[1],
        min_normalized_cn = min_normalized_cn[1],
        min_cn = min_cn[1],
        avg_normalized_cn = avg_normalized_cn[1],
        avg_cn = avg_cn[1],
        number_of_cn_segments = number_of_cn_segments[1],
        gene_width = gene_width[1],
        total_node_width = sum(width) # This must be calculated after nominating SCNA type, not before.,
      ), 
      by = .(gene_name, type)]
  )

  scna_result = base::subset(
    scna_result,
    !is.na(scna_result$type)
  )

  return(scna_result)
}

#' @title collect_copy_number_jabba
#' @description
#' Collects copy number and jabba data from a specified file and processes it.
#'
#' @param jabba_rds Path to the jabba.simple.rds file.
#' @param pge GRanges object with gencode annotations.
#' @param amp.thresh SCNA amplification threshold.
#' @param del.thresh SCNA deletion threshold.
#' @param verbose Logical flag to indicate if messages should be printed.
#' @param karyograph Optional path to the karyograph.rds file.
#' @return A data.table containing processed copy number and jabba information.
collect_copy_number_jabba <- function(
    jabba_rds,
    pge,
    amp.thresh,
    del.thresh,
    verbose = TRUE,
    karyograph = NULL) {
  if (is.null(jabba_rds) || !file.exists(jabba_rds)) {
    if (verbose) message("Jabba RDS file is missing or does not exist.")
    return(data.table(type = NA, source = "jabba_rds"))
  }

  if (verbose) message("pulling jabba_rds to get SCNA and purity / ploidy")
  jab <- readRDS(jabba_rds)
  jabpurity = (
    base::get0(
      "purity", 
      as.environment(jab$meta), 
      ifnotfound = base::get0(
        "purity", as.environment(jab), 
        ifnotfound = stop("purity not found in jabba object")
      )
    )
  )
  jabploidy = (
    base::get0(
      "ploidy", 
      as.environment(jab$meta), 
      ifnotfound = base::get0(
        "ploidy", as.environment(jab), 
        ifnotfound = stop("ploidy not found in jabba object")
      )
    )
  )
  result <- data.table(
    value = c(jabpurity, jabploidy),
    type = c("purity", "ploidy"),
    track = "pp"
  )

  # get the ncn data from jabba
  nseg <- NULL

  # Don't fail out if karyograph isn't found
  # Just provide nseg as NULL, and get_gene_ampdels_from_jabba
  # assumes ncn = 2 (pretty safe assumption)
  if (!is.null(karyograph) && file.exists(karyograph)) {
    nseg <- readRDS(karyograph)$segstats[, c("ncn")]
  }

  scna <- get_gene_ampdels_from_jabba(
    jab,
    amp.thresh = amp.thresh,
    del.thresh = del.thresh,
    pge = pge,
    nseg = nseg
  )

  if (nrow(scna)) {
    scna[, track := "variants"][, source := "jabba_rds"][, vartype := "scna"]
    result <- rbind(
      result,
      scna[, .(value = min_cn, type, track, gene = gene_name)],
      fill = TRUE,
      use.names = TRUE
    )
  }

  return(result)
}

#' @description
#' Collects signature data from a specified file and processes it.
#'
#' @param signature_counts Path to the signature_counts.txt file.
#' @param verbose Logical flag to indicate if messages should be printed.
#' @return A data.table containing processed signature information.
collect_signatures <- function(signature_counts, verbose = TRUE) {
  # TODO: update this function to use sigprofiler
  if (!is.null(signature_counts) && file.exists(signature_counts)) {
    if (verbose) message("pulling signature_counts")
    sig <- fread(signature_counts)
    sig <- sig[, .(value = num_events, type = Signature, etiology = Etiology, frac = frac.events, track = "signature", source = "signature_counts")]
    return(sig)
  } else {
    return(data.table(type = NA, source = "signature_counts"))
  }
}

#' @title collect_gene_mutations
#' @description
#' Collects gene mutation data from a specified file and processes it.
#'
#' @param annotated_bcf Path to the annotated.bcf file.
#' @param jabba_rds Path to the jabba.simple.rds file.
#' @param filter Filter to apply to the variants.
#' @param verbose Logical flag to indicate if messages should be printed.
#' @return A data.table containing processed gene mutation information.
collect_gene_mutations <- function(
    annotated_bcf,
    jabba_rds,
    filter,
    verbose = TRUE) {
  if (is.null(annotated_bcf) || !file.exists(annotated_bcf)) {
    if (verbose) message("Annotated BCF file is missing or does not exist.")
    return(data.table(type = NA, source = "annotated_bcf"))
  }

  if (verbose) message("pulling annotated_bcf using FILTER=", filter)

  local_bcftools_path <- Sys.which("bcftools")
  if (local_bcftools_path == "") {
    stop("bcftools not found in the system PATH. Please install or module load bcftools.")
  }
  if (verbose) message("bcftools found at: ", local_bcftools_path)

  bcf <- skitools::grok_bcf(
    annotated_bcf,
    label = "sample",
    long = TRUE,
    filter = filter,
    bpath = local_bcftools_path
  )
  if (verbose) message(length(bcf), " variants pass filter")

  genome.size <- sum(GenomeInfoDb::seqlengths(bcf), na.rm = TRUE) / 1e6
  if (is.na(genome.size)) {
    genome.size <- sum(
      GenomeInfoDb::seqlengths(gGnome::gG(jabba = jabba_rds)),
      na.rm = TRUE
    ) / 1e6
  }

  nmut <- unique(data.table(
    as.character(seqnames(bcf)),
    start(bcf),
    end(bcf),
    bcf$REF,
    bcf$ALT
  ))[, .N]
  mut.density <- data.table(
    value = c(nmut, nmut / genome.size),
    type = c("count", "density"),
    track = "tmb",
    source = "annotated_bcf"
  )

  keepeff <- c("trunc", "cnadel", "cnadup", "complexsv", "splice", "inframe_indel", "fusion", "missense", "promoter", "regulatory", "mir")
  bcf <- bcf[bcf$short %in% keepeff]
  if (verbose) message(length(bcf), " variants pass keepeff")

  if (length(bcf) == 0) {
    return(mut.density)
  }

  bcf$variant.g <- paste0(seqnames(bcf), ":", start(bcf), "-", end(bcf), " ", bcf$REF, ">", bcf$ALT)
  vars <- gr2dt(bcf)[, .(gene, vartype, variant.g, variant.p, distance, annotation, type = short, track = "variants", source = "annotated_bcf")]
  setkey(vars, variant.g)
  vars <- vars[, .SD[1], by = variant.g]

  return(rbind(mut.density, vars, fill = TRUE, use.names = TRUE))
}

#' @title collect_oncokb_cna
#' @description
#' Collects OncoKB CNA data from a specified file and processes it.
#'
#' @param oncokb_cna Path to the oncokb CNA file.
#' @param verbose Logical flag to indicate if messages should be printed.
#' @return A data.table containing processed OncoKB CNA information.
collect_oncokb_cna <- function(oncokb_cna, jabba_gg, pge, amp.thresh, del.thresh, karyograph = NULL, verbose = TRUE) {
  if (is.null(oncokb_cna) || !file.exists(oncokb_cna)) {
    if (verbose) message("OncoKB CNA file is missing or does not exist.")
    return(data.table(type = NA, source = "oncokb_cna"))
  }

  oncokb_cna <- data.table::fread(oncokb_cna)

  nseg = NULL
  if (!is.null(karyograph) && file.exists(karyograph)) {
      nseg <- readRDS(karyograph)$segstats[, c("ncn")]
  }

  scna <- get_gene_ampdels_from_jabba(
      jabba_gg,
      amp.thresh = amp.thresh,
      del.thresh = del.thresh,
      pge = pge,
      nseg = nseg
  ) 

  matches = list(
      c("Deletion", "homdel"),
      c("Amplification", "amp")
  )
  is_valid_oncokb_cna = rep(FALSE, NROW(oncokb_cna))
  for (match_lst in matches) {
      valid_gene = scna[
          scna$total_node_width > 1e3
          & scna$type == match_lst[2],

      ]$gene_name
      is_valid_oncokb_cna = (
          is_valid_oncokb_cna
          | (oncokb_cna$HUGO_SYMBOL %in% valid_gene & oncokb_cna$ALTERATION %in% match_lst[1])
      )
  }
  oncokb_cna = oncokb_cna[is_valid_oncokb_cna,]

  if (NROW(oncokb_cna) > 0) {
    oncokb_cna <- parse_oncokb_tier(
      oncokb_cna,
      tx_cols = c("LEVEL_1", "LEVEL_2"),
      rx_cols = c("LEVEL_R1"),
      dx_cols = c("LEVEL_Dx1"),
      px_cols = c("LEVEL_Px1")
    )
    return(oncokb_cna[, .(
      gene = HUGO_SYMBOL,
      gene_summary = GENE_SUMMARY,
      role = Role,
      value = min_cn,
      type = ifelse(
        ALTERATION == "Amplification",
        "amp",
        ifelse(ALTERATION == "Deletion", "homdel", NA_character_)
      ),
      tier = tier,
      tier_description = tier_factor,
      variant_summary = VARIANT_SUMMARY,
      therapeutics = tx_string,
      resistances = rx_string,
      diagnoses = dx_string,
      prognoses = px_string,
      effect = MUTATION_EFFECT,
      effect_description = MUTATION_EFFECT_DESCRIPTION,
      track = "scna",
      source = "oncokb_cna"
    )])
  }
  return(data.table(type = NA, source = "oncokb_cna"))
}


#' @title collect_oncokb_fusions
#' @description
#' Collects OncoKB Fusion data from a specified file and processes it.
#'
#' @param oncokb_fusions Path to the oncokb fusions file.
#' @param verbose Logical flag to indicate if messages should be printed.
#' @return A data.table containing processed OncoKB Fusion information.
#' @author Kevin Hadi
collect_oncokb_fusions <- function(oncokb_fusions, pge, cytoband, verbose = TRUE) {
  out <- data.table(vartype = NA, source = "oncokb_fusions")
  if (is.null(oncokb_fusions) || !file.exists(oncokb_fusions)) {
    if (verbose) message("OncoKB Fusions file is missing or does not exist.")
    return(out)
  }

  oncokb_fusions <- data.table::fread(oncokb_fusions)

  if (NROW(oncokb_fusions) > 0) {
    oncokb_fusions <- parse_oncokb_tier(
      oncokb_fusions,
      tx_cols = c("LEVEL_1", "LEVEL_2"),
      rx_cols = c("LEVEL_R1"),
      dx_cols = c("LEVEL_Dx1"),
      px_cols = c("LEVEL_Px1")
    )

    non_silent_fusions <- oncokb_fusions[silent == FALSE, ] # already de-duped
    non_silent_fusions[, vartype := ifelse(in.frame == TRUE, "fusion", "outframe_fusion")]

    if (!NROW(non_silent_fusions) > 0) {
      return(out)
    }

    genes <- strsplit(non_silent_fusions$FUSION, "-")
    genes_matrix <- do.call(rbind, genes)
    ixA <- match(genes_matrix[, 1], pge$gene_name)
    ixB <- match(genes_matrix[, 2], pge$gene_name)
    # remove NA from ixA and ixB (needed for tests that use subset of gencode genes)
    na.index <- which(is.na(ixA))
    ixA <- ixA[!is.na(ixA)]
    ixB <- ixB[!is.na(ixB)]
    grA <- pge[ixA]
    grB <- pge[ixB]
    coordA <- gUtils::gr.string(grA)
    coordB <- gUtils::gr.string(grB)
    grovA <- gUtils::gr.findoverlaps(grA, cytoband, scol = "chromband")
    grovB <- gUtils::gr.findoverlaps(grB, cytoband, scol = "chromband")

    if (length(grovA) > 0) {
      grovA <- GenomicRanges::sort(grovA) %Q% (order(query.id, ifelse(grepl("p", chromband), -1, 1) * start))
      non_silent_fusions$cytoA <- ifelse(!1:nrow(non_silent_fusions) %in% na.index,
        gUtils::gr2dt(grovA)[, paste(unique(chromband[chromband %in% c(chromband[1], tail(chromband, 1))]), collapse = "-"), by = query.id]$V1, NA
      )
    } else {
      non_silent_fusions$cytoA <- ""
    }

    if (length(grovB) > 0) {
      grovB <- GenomicRanges::sort(grovB) %Q% (order(query.id, ifelse(grepl("p", chromband), -1, 1) * start))
      non_silent_fusions$cytoB <- ifelse(!1:nrow(non_silent_fusions) %in% na.index,
        gUtils::gr2dt(grovB)[, paste(unique(chromband[chromband %in% c(chromband[1], tail(chromband, 1))]), collapse = "-"), by = query.id]$V1, NA
      )
    } else {
      non_silent_fusions$cytoB <- ""
    }

    non_silent_fusions$fusion_genes <- paste0(genes_matrix[, 1], "(", non_silent_fusions$exonA, ")::", genes_matrix[, 2], "(", non_silent_fusions$exonB, ")@", non_silent_fusions$cytoA, "::", non_silent_fusions$cytoB)
    non_silent_fusions$fusion_gene_coords <- ifelse(!1:nrow(non_silent_fusions) %in% na.index,
      paste(coordA, coordB, sep = ","),
      NA
    )
    out <- non_silent_fusions[, .(
      gene = Hugo_Symbol,
      gene_summary = GENE_SUMMARY,
      role = Role,
      value = min_cn,
      vartype,
      type = "fusion",
      tier = tier,
      tier_description = tier_factor,
      variant_summary = VARIANT_SUMMARY,
      therapeutics = tx_string,
      resistances = rx_string,
      diagnoses = dx_string,
      prognoses = px_string,
      effect = MUTATION_EFFECT,
      effect_description = MUTATION_EFFECT_DESCRIPTION,
      fusion_genes,
      fusion_gene_coords,
      track = "variants",
      source = "oncokb_fusions"
    )]
  }

  return(out)
}

#'
#' Helper function to parse oncokb outputs
#' levels of evidence and assign tier.
#' Tiering is simply:
#' 1 = Clinically Actionable:
#' Therapeutic sensitivity, resistance,
#' diagnostic, or prognostic information is assigned
#' to the variant.
#' 2 = Clinically Significant:
#' Clinically significant variant entails that the variant
#' does not have therapeutic, diagnostic, or prognostic effect
#' validated at FDA level, but is oncogenic.
#' 3 = All others (VUS)
#'
#' @param oncokb data.table object holding oncokb outputs
#' @author Kevin Hadi
parse_oncokb_tier <- function(
    oncokb,
    tx_cols = c("LEVEL_1", "LEVEL_2"),
    rx_cols = c("LEVEL_R1"),
    dx_cols = c("LEVEL_Dx1"),
    px_cols = c("LEVEL_Px1")) {
  .concat_string <- function(oncokb, cols) {
    out_string <- lapply(base::subset(oncokb, select = cols), function(y) (strsplit(y, ",")))
    concat_out <- IRanges::CharacterList(Reduce(concat_vectors, out_string))
    concat_out <- S4Vectors::unique(concat_out)
    concat_out <- concat_out[!is.na(concat_out)]
    concat_out[S4Vectors::elementNROWS(concat_out) == 0] <- NA_character_
    concat_out <- stringi::stri_c_list(as.list(concat_out), sep = ",")
    return(concat_out)
  }
  is_actionable <- logical(NROW(oncokb))
  for (col in c(tx_cols, rx_cols, dx_cols, px_cols)) {
    oncokb[[col]] <- as.character(oncokb[[col]])
    is_actionable <- is_actionable | (!is.na(oncokb[[col]]) & base::nzchar(oncokb[[col]]))
  }
  oncokb$is_actionable <- is_actionable
  oncokb$is_oncogenic <- oncokb$ONCOGENIC %in% c("Likely Oncogenic", "Oncogenic")
  tier_factor <- ifelse(
    oncokb$is_actionable, "Clinically Actionable",
    ifelse(oncokb$is_oncogenic, "Clinically Significant", "VUS")
  ) %>% factor(c("Clinically Actionable", "Clinically Significant", "VUS"))

  oncokb$tier_factor <- tier_factor
  oncokb$tier <- as.integer(tier_factor)
  oncokb$tx_string <- .concat_string(oncokb, tx_cols)
  oncokb$rx_string <- .concat_string(oncokb, rx_cols)
  oncokb$dx_string <- .concat_string(oncokb, dx_cols)
  oncokb$px_string <- .concat_string(oncokb, px_cols)

  return(oncokb)
}

#' @title collect_oncokb
#' @description
#' Collects OncoKB mutation data from a specified file and processes it.
#'
#' @param oncokb_maf Path to the oncokb MAF file.
#' @param multiplicity Path to the multiplicity file.
#' @param verbose Logical flag to indicate if messages should be printed.
#' @return A data.table containing processed OncoKB mutation information.
collect_oncokb <- function(oncokb_maf, multiplicity = NA_character_, verbose = TRUE) {
  if (is.null(oncokb_maf) || !file.exists(oncokb_maf)) {
    if (verbose) message("OncoKB MAF file is missing or does not exist.")
    return(data.table(type = NA, source = "oncokb_maf"))
  }

  # snpeff_ontology = readRDS(system.file("extdata", "data", "snpeff_ontology.rds", package = "Skilift"))
  oncokb <- data.table::fread(oncokb_maf)
  if (
    is.character(multiplicity) &&
      NROW(multiplicity) == 1 &&
      file.exists(multiplicity)
  ) {
    multiplicity <- readRDS(multiplicity)
    oncokb <- merge_oncokb_multiplicity(oncokb, multiplicity, overwrite = TRUE)
  }

  if (NROW(oncokb) > 0) {
    ## oncokb$snpeff_ontology <- snpeff_ontology$short[match(oncokb$Consequence, snpeff_ontology$eff)]
    oncokb$short <- dplyr::case_when(
      grepl("frameshift", oncokb$Consequence) & grepl("fs$", oncokb$HGVSp) ~ "trunc",
      grepl("stop", oncokb$Consequence) & grepl("^p\\.", oncokb$HGVSp) ~ "trunc",
      grepl("lost", oncokb$Consequence) ~ "trunc",
      grepl("missense", oncokb$Consequence) & grepl("^p\\.", oncokb$HGVSp) ~ "missense",
      grepl("splice", oncokb$Consequence) ~ "splice",
      grepl("(5|3)_prime_UTR_variant", oncokb$Consequence) ~ "UTR",
      TRUE ~ NA_character_
    )

    oncokb <- parse_oncokb_tier(
      oncokb,
      tx_cols = c("LEVEL_1", "LEVEL_2"),
      rx_cols = c("LEVEL_R1"),
      dx_cols = c("LEVEL_Dx1"),
      px_cols = c("LEVEL_Px1")
    )

    # coerce_column_tuples <- list(
    #   c("altered_copies", "total_copies")
    # )
    # current_nms <- names(oncokb)
    # for (col in coerce_column_tuples) {
    #   # if (col[1] %in% current_nms && !col[2] %in% current_nms) {
    #   current_nms[current_nms == col[1]] <- col[2]
    #   # }
    # }
    # names(oncokb) <- current_nms
    return(oncokb[, .(
      gene = Hugo_Symbol,
      gene_summary = GENE_SUMMARY,
      role = Role,
      variant.g = paste(Chromosome, ":", Start_Position, "-", End_Position, " ", variant.g, sep = ""),
      variant.c = HGVSc,
      variant.p = HGVSp,
      annotation = Consequence,
      type = short,
      tier = tier,
      tier_description = tier_factor,
      variant_summary = VARIANT_SUMMARY,
      therapeutics = tx_string,
      resistances = rx_string,
      diagnoses = dx_string,
      prognoses = px_string,
      distance = NA_integer_,
      effect = MUTATION_EFFECT,
      effect_description = MUTATION_EFFECT_DESCRIPTION,
      major_count = major.count,
      minor_count = minor.count,
      major_snv_copies,
      minor_snv_copies,
      # total_copies, ## total_copies gets converted to estimated_altered_copies downstream
      altered_copies,
      segment_cn,
      ref,
      alt,
      VAF,
      vartype = "SNV",
      track = "variants",
      source = "oncokb_maf"
    )])
  }
  return(data.table(type = NA, source = "oncokb_maf"))
}


#' @title oncotable
#' @description
#'
#' @param somatic_variant_annotations Path to annotated.bcf file
#' @param fusions Path to fusion.rds file
#' @param jabba_gg Path to jabba.simple.rds file
#' @param events Path to complex.rds file
#' @param signature_counts Path to signature_counts.txt file
#' @param karyograph Optional path to the karyograph.rds file
#' @param oncokb_snv Path to oncokb MAF file
#' @param oncokb_cna Path to oncokb CNA file
#' @param gencode_gr GRanges object with gencode annotations
#' @param amp.thresh SCNA amplification threshold to call an amp as a function of ploidy (4)
#' @param del.thresh SCNA deletion threshold for (het) del as a function of ploidy (by default cn = 1 will be called del, but this allows additoinal regions in high ploidy tumors to be considered het dels)
#' @param verbose logical flag
#' @export
oncotable <- function(
    pair,
    somatic_variant_annotations = NULL,
    fusions = NULL,
    jabba_gg = NULL,
    karyograph = NULL,
    events = NULL,
    signature_counts = NULL,
    multiplicity = NULL,
    oncokb_snv = NULL,
    oncokb_cna = NULL,
    oncokb_fusions = NULL,
    gencode,
    cytoband,
    verbose = TRUE,
    amp.thresh = 4,
    filter = "PASS",
    del.thresh = 0.5) {
  out <- data.table()

  if ("type" %in% names(mcols(gencode))) {
    pge <- gencode %Q% (type == "gene" & gene_type == "protein_coding")
  } else {
    pge <- gencode %Q% (gene_type == "protein_coding")
  }

  ## collect gene fusions
  # prefer fusions from oncokb
  if (!is.null(oncokb_fusions) && file.exists(oncokb_fusions)) {
    out <- rbind(
      out,
      collect_oncokb_fusions(oncokb_fusions, pge, cytoband, verbose),
      fill = TRUE,
      use.names = TRUE
    )
  } else {
    out <- rbind(
      out,
      collect_gene_fusions(fusions, pge, verbose),
      fill = TRUE,
      use.names = TRUE
    )
  }

  ## collect complex events
  out <- rbind(
    out,
    collect_complex_events(events, verbose),
    fill = TRUE,
    use.names = TRUE
  )

  ## collect copy number
  # prefer oncokb cna
  if (!is.null(oncokb_cna) && file.exists(oncokb_cna)) {
    out <- rbind(
      out,
      collect_oncokb_cna(oncokb_cna, jabba_gg, pge, amp.thresh, del.thresh, karyograph, verbose),
      fill = TRUE,
      use.names = TRUE
    )
  } else {
    out <- rbind(
      out,
      collect_copy_number_jabba(jabba_gg, pge, amp.thresh, del.thresh, verbose, karyograph),
      fill = TRUE,
      use.names = TRUE
    )
  }

  ## collect gene mutations
  # prefer oncokb snv
  if (!is.null(oncokb_snv) && file.exists(oncokb_snv)) {
    out <- rbind(
      out,
      collect_oncokb(oncokb_snv, multiplicity, verbose),
      fill = TRUE,
      use.names = TRUE
    )
  } else {
    out <- rbind(
      out,
      collect_gene_mutations(somatic_variant_annotations, jabba_gg, filter, verbose),
      fill = TRUE,
      use.names = TRUE
    )
  }

  ## add gene locations
  gene_locations <- readRDS(system.file("extdata", "data", "gene_locations.rds", package = "Skilift"))

  # Merge gene locations and create location string
  if ("gene" %in% names(out)) {
    out[gene_locations,
      gene_location := paste0(i.seqnames, ":", i.start, "-", i.end),
      on = c(gene = "gene_name")
    ]
  }

  out$id <- pair

  # remove all rows for which data was not passed
  out <- out[!is.na(type)]

  # coerce all empty strings to NA without affecting columns of levels
  out <- out[, lapply(.SD, function(x) {
    if (is.factor(x)) x <- as.character(x) # Convert factors to characters
    ifelse(x == "", NA_character_, x) # Replace empty strings with NA
  })]

  if (verbose) message("done processing sample")
  return(out)
}

#' @name create_oncotable
#' @title create_oncotable
#' @description
#'
#' function to create oncotable for use with filtered_events_json
#'
#' @param cohort Cohort object containing sample information
#' @param amp_thresh_multiplier amp.thresh for oncotable is amp_thresh_multiplier*ploidy
#' @param gencode file to gencode annotations (uses v29lift37 by default)
#' @param outdir path to directory in which to write oncotable outputs
#' @param cores number of cores for parallel processing
#' @return None
#' @export
#' @author Shihab Dider, Joel Rosiene
create_oncotable <- function(
    cohort,
    amp_thresh_multiplier = 1.5,
    gencode = "/gpfs/data/imielinskilab/DB/GENCODE/gencode.v19.annotation.gtf.nochr.rds",
    cytoband = system.file("extdata", "data", "cytoband.rds", package = "Skilift"),
    outdir,
    cores = 1) {
  if (!inherits(cohort, "Cohort")) {
    stop("Input must be a Cohort object")
  }

  if (nrow(cohort$inputs) == 0) {
    stop("No samples found in the cohort")
  }

  if (Sys.which("bcftools") == "") {
    stop("bcftools is not available on the system PATH. Try `module load htslib` first or install it.")
  } else {
    message("bcftools is available.")
  }

  if (gencode == "/gpfs/data/imielinskilab/DB/GENCODE/gencode.v19.annotation.gtf.nochr.rds") {
    message("using default gencode: /gpfs/data/imielinskilab/DB/GENCODE/gencode.v19.annotation.gtf.nochr.rds")
  }

  gencode <- process_gencode(gencode)
  cytoband <- process_cytoband(cytoband)

  if (amp_thresh_multiplier == 1.5) {
    message("using default amp_thres_multiplier: 1.5")
  }

  # Create output directory if it doesn't exist
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  # Create a copy of the cohort to modify
  updated_cohort <- cohort

  # Add oncotable column if it doesn't exist
  if (!"oncotable" %in% names(updated_cohort$inputs)) {
    updated_cohort$inputs[, oncotable := NA_character_]
  }
  
  results <- mclapply(seq_len(nrow(cohort$inputs)), function(i) {
  (
      {
        row <- cohort$inputs[i]

        # Create output directory for this pair
        pair_outdir <- file.path(outdir, row$pair)
        if (!dir.exists(pair_outdir)) {
          dir.create(pair_outdir, recursive = TRUE)
        }

        # Get ploidy from jabba output, default to 2 if missing
        ploidy <- 2 # Default ploidy
        if (file.exists(row$jabba_gg)) {
          ploidy_ggraph <- tryCatch(
            {
              readRDS(row$jabba_gg)
            },
            error = function(e) {
              msg <- sprintf("Error reading JaBbA file for %s: %s. Using default ploidy of 2.", row$pair, e$message)
              warning(msg)
            }
          )

          if (!is.null(ploidy_ggraph)) {
            ploidy <- ifelse(
              !is.null(ploidy_ggraph$meta$ploidy),
              ploidy_ggraph$meta$ploidy,
              ploidy_ggraph$ploidy
            )
          }
        }

        amp_thresh <- amp_thresh_multiplier * ploidy
        message(paste("Processing", row$pair, "using amp.thresh of", amp_thresh))

        futile.logger::flog.threshold("ERROR")
        # Run oncotable for this pair
        oncotable_result <- tryCatchLog(
          {
            oncotable(
              pair = row$pair,
              somatic_variant_annotations = row$somatic_variant_annotations,
              fusions = row$fusions,
              jabba_gg = row$jabba_gg,
              karyograph = row$karyograph,
              events = row$events,
              signature_counts = row$signature_counts,
              multiplicity = row$multiplicity,
              oncokb_snv = row$oncokb_snv,
              oncokb_cna = row$oncokb_cna,
              oncokb_fusions = row$oncokb_fusions,
              gencode = gencode,
              cytoband = cytoband,
              verbose = TRUE,
              amp.thresh = amp_thresh,
              filter = "PASS",
              del.thresh = 0.5
            )
          },
          error = function(e) {
            msg <- sprintf("Error in oncotable for %s: %s", row$pair, e$message)
            print(msg)
            return(NULL)
          }
        )

        if (!is.null(oncotable_result)) {
          # Save successful results
          oncotable_path <- file.path(pair_outdir, "oncotable.rds")
          saveRDS(oncotable_result, oncotable_path)
          fwrite(oncotable_result, file.path(pair_outdir, "oncotable.txt"))
          return(list(index = i, path = oncotable_path))
        }
      }
    )
  }, mc.cores = cores, mc.preschedule = TRUE)


  # Update oncotable paths in the cohort
  results <- Filter(Negate(is.null), results)
  for (result in results) {
    updated_cohort$inputs[result$index, oncotable := result$path][]
  }
  message(sprintf("\nProcessing complete - results written to %s", outdir))

  return(updated_cohort)
}


#' @name create_filtered_events
#' @title create_filtered_events
#' @description
#' Create filtered events for a single sample
#'
#' @param pair patient id
#' @param oncotable oncotable task output
#' @param jabba_gg JaBbA output ggraph or complex
#' @param out_file path to write json
#' @param return_table TRUE/FALSE whether to return the data.table
#' @return data.table or NULL
#' @export
create_filtered_events <- function(
    pair,
    oncotable,
    jabba_gg,
    out_file,
    return_table = FALSE,
    cohort_type = "paired") {
  ot <- readRDS(oncotable)

  possible_drivers = empty_oncotable = structure(list(gene = character(0), gene_summary = character(0), 
    role = character(0), value = integer(0), vartype = character(0), 
    type = character(0), tier = integer(0), tier_description = character(0), 
    variant_summary = character(0), therapeutics = character(0), 
    resistances = logical(0), diagnoses = character(0), prognoses = character(0), 
    effect = character(0), effect_description = character(0), 
    fusion_genes = character(0), fusion_gene_coords = character(0), 
    track = character(0), source = character(0), variant.g = character(0), 
    variant.c = character(0), variant.p = character(0), annotation = character(0), 
    distance = logical(0), major_count = numeric(0), minor_count = numeric(0), 
    major_snv_copies = numeric(0), minor_snv_copies = numeric(0), 
    altered_copies = numeric(0), segment_cn = integer(0), ref = integer(0), 
    alt = integer(0), VAF = numeric(0), gene_location = character(0), 
    id = character(0)), row.names = c(NA, 0L), class = c("data.table", 
"data.frame"))

  if (NROW(ot) > 0) {
      # add a fusion_gene_coords column of NAs if no fusions
      if (!"fusion_gene_coords" %in% colnames(ot)) {
        ot[, fusion_genes := NA]
        ot[, fusion_gene_coords := NA]
      }
      # snvs <- ot[grepl("frameshift|missense|stop|disruptive", annotation, perl = TRUE)]
      snvs <- ot[ot$vartype == "SNV" & !is.na(ot$type)]
      if ("tier" %in% colnames(ot)) {
        snvs <- snvs[order(is.na(tier)), ]
      }
      if (NROW(snvs) > 0) {
        snvs[, is_unique_p := !is.na(variant.p) & !duplicated(cbind(gene, variant.p))]
        snvs[, is_unique_g := !duplicated(cbind(gene, variant.g))]
        remove_variant_c = FALSE
        if (is.null(snvs$variant.c)) {
          snvs$variant.c = 1:NROW(snvs)
          remove_variant_c = TRUE
        }
        snvs[, is_unique_c := !duplicated(cbind(gene, variant.c))]
        snvs <-  snvs[is_unique_p | (is_unique_g & is_unique_c)]
        snvs$is_unique_p = NULL
        snvs$is_unique_g = NULL
        snvs$is_unique_c = NULL
        if (remove_variant_c) snvs$variant.c = NULL
      }

      homdels <- ot[ot$type == "homdel",][, vartype := "HOMDEL"][, type := "SCNA"]
      amps <- ot[ot$type == "amp",][, vartype := "AMP"][, type := "SCNA"]
      fusions <- ot[ot$type == "fusion",] ## fusion vartype is either fusion or outframe_fusion
      possible_drivers <- rbind(snvs, homdels, amps, fusions)
  }

  oncotable_col_to_filtered_events_col <- c(
    "id" = "id",
    "gene" = "gene",
    "fusion_genes" = "fusion_genes",
    "fusion_gene_coords" = "fusion_gene_coords",
    "value" = "fusion_cn",
    "vartype" = "vartype",
    "type" = "type",
    "variant.g" = "Variant_g",
    "variant.p" = "Variant",
    "altered_copies" = "estimated_altered_copies",
    "segment_cn" = "segment_cn",
    "ref" = "ref",
    "alt" = "alt",
    "VAF" = "VAF",
    "gene_location" = "Genome_Location",
    "tier" = "Tier",
    "role" = "role",
    "gene_summary" = "gene_summary",
    "variant_summary" = "variant_summary",
    "effect" = "effect",
    "effect_description" = "effect_description",
    "therapeutics" = "therapeutics",
    "resistances" = "resistances",
    "diagnoses" = "diagnoses",
    "prognoses" = "prognoses"
  )
  filtered_events_columns <- names(possible_drivers)[names(possible_drivers) %in% names(oncotable_col_to_filtered_events_col)]

  res <- possible_drivers[, ..filtered_events_columns]
  intersected_columns <- intersect(filtered_events_columns, names(res))
  setnames(res, old = intersected_columns, new = oncotable_col_to_filtered_events_col[intersected_columns])

  res <- res %>% unique(., by = c("gene", "vartype", "Variant"))
  if (nrow(res) > 0) {
    res[, seqnames := tstrsplit(Genome_Location, ":", fixed = TRUE, keep = 1)]
    res[, start := tstrsplit(Genome_Location, "-", fixed = TRUE, keep = 1)]
    res[, start := tstrsplit(start, ":", fixed = TRUE, keep = 2)]
    res[, end := tstrsplit(Genome_Location, "-", fixed = TRUE, keep = 2)]
    res.mut <- res[vartype == "SNV"]
    if (nrow(res.mut) > 0) {
      # res.mut[, Variant := gsub("p.", "", Variant)]
      res.mut[, vartype := "SNV"]
      # TODO:
      # truncating mutations are not always deletions
      # initial logic may be misleading calling all small mutations "SNV"
      # but we should encode this as something more robust
      # res.mut[type=="trunc", vartype := "DEL"]
    }
    res.fus = res[type == "fusion",] ## need to deal each class explicitly
    if (NROW(res.fus) > 0) {
      res.fus$gene = res.fus$fusion_genes
      is_inframe = res.fus$vartype == "fusion"
      is_outframe = res.fus$vartype == "outframe_fusion"
      res.fus$Variant = ifelse(
        is_inframe,
        "In-Frame Fusion",
        ifelse(
          is_outframe,
          "Out-of-Frame Fusion",
          NA_character_
        )
      )

      res.fus$estimated_altered_copies = res.fus$fusion_cn
    }
    res.cn.dt = res.cn <- res[(
      type == "SCNA" ## FIXME: redundant logic for now
      | vartype %in% c("AMP", "HOMDEL")
    )]
    if (nrow(res.cn) > 0) {
      # jab <- readRDS(jabba_gg)
      # res.cn.gr <- GRanges(res.cn)
      # res.cn.dt
      # res.cn.gr <- gr.val(res.cn.gr, jab$nodes$gr, c("cn", "cn.low", "cn.high"))
      # res.cn.dt <- as.data.table(res.cn.gr)
      # FIXME: value is actually min_cn is renamed fusion_cn!!
      # Also, we already picked out the corresponding segment cn to use (min_cn)
      res.cn.dt$min_cn = res.cn.dt$fusion_cn 
      res.cn.dt[, estimated_altered_copies := abs(min_cn - 2)]
      res.cn.dt[, segment_cn := min_cn]
      res.cn.dt[, Variant := vartype]
      res.cn.dt$fusion_cn = NA_real_

      # To fix very small CNAs that can appear due to forcing in reciprocal
      # junctions with very small gaps.
      ## res.cn.dt = res.cn.dt[!res.cn.dt$node_width_jab < 1000,]
      # remove redundant columns since already added to Variant
      
      res.cn.dt[, c("min_cn", "cn", "cn.high", "cn.low", "width", "strand") := NULL]
    }
    res.final <- rbind(res.mut, res.cn.dt, res.fus, fill = TRUE)
    res.final[, sample := pair]
    if (identical(cohort_type, "heme")) {
      res.final <- select_heme_events(res.final)
    }
    res.final$type = tools::toTitleCase(res.final$type)
    ### FIXME: REMOVING Y CHROMOSOME UNTIL WE UPDATE DRYCLEAN
    res.final <- res.final[res.final$seqnames != "Y"]
    ### FIXME ^^^: REMOVING Y CHROMOSOME UNTIL WE UPDATE DRYCLEAN
    write_json(res.final, out_file, pretty = TRUE)
    
    ## FIXME: Decide whether to either propagate return_table to lift_mvp in lift_heme somehow, or change this
    ## conditional logic to be based on cohort_type. Can also be left alone.
    ## Note that return_table for debugging purposes as well.
    if (return_table || identical(cohort_type, "heme")) {
      return(res.final)
    }
  }
  NULL
}


#' @name lift_filtered_events
#' @title lift_filtered_events
#' @description
#' Create filtered events for all samples in a cohort
#'
#' @param cohort Cohort object containing sample information
#' @param output_data_dir Base directory for output files
#' @param cores Number of cores for parallel processing (default: 1)
#' @return None
#' @export
lift_filtered_events <- function(cohort, output_data_dir, cores = 1, return_table = FALSE) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }
    
    if (!dir.exists(output_data_dir)) {
        dir.create(output_data_dir, recursive = TRUE)
    }
    
    # Validate required columns exist
    required_cols <- c("pair", "oncotable", "jabba_gg")
    missing_cols <- required_cols[!required_cols %in% names(cohort$inputs)]
    if (length(missing_cols) > 0) {
        stop("Missing required columns in cohort: ", paste(missing_cols, collapse = ", "))
    }
    
    cohort_type = cohort$type
    # Process each sample in parallel
    lst_outs = mclapply(seq_len(nrow(cohort$inputs)), function(i) {
        row <- cohort$inputs[i,]
        pair_dir <- file.path(output_data_dir, row$pair)
        
        if (!dir.exists(pair_dir)) {
            dir.create(pair_dir, recursive = TRUE)
        }
        
        out_file <- file.path(pair_dir, "filtered.events.json")
        highlights_out_file <- file.path(pair_dir, "highlights.json")

        out = NULL
        futile.logger::flog.threshold("ERROR")
        tryCatchLog({
          out <- create_filtered_events(
              pair = row$pair,
              oncotable = row$oncotable,
              jabba_gg = row$jabba_gg,
              out_file = out_file,
              return_table = return_table,
              cohort_type = cohort_type
          )
          if (identical(cohort_type, "heme")) {
            create_heme_highlights(events_tbl = out, jabba_gg = row$jabba_gg, out_file = highlights_out_file)
          }
            
        }, error = function(e) {
            print(sprintf("Error processing %s: %s", row$pair, e$message))
            NULL
        }
        )

        
        return(out)
    }, mc.cores = cores, mc.preschedule = TRUE)
    
    invisible(lst_outs)
}


#' Global HemeDB path
HEMEDB = "/gpfs/data/imielinskilab/projects/Clinical_NYU/db/master_heme_database.20250128_095937.790322.rds"

#' Select Heme events from Addy's hemedb
#' 
#' Coarse selection of Heme events for first pass filtering
#'
#' @export
select_heme_events <- function(
  filtered_events, 
  hemedb_path = Skilift:::HEMEDB
) {
  hemedb = readRDS(hemedb_path) 

  # wtf = merge(events_tbl, hemedb_genes, by.x = "gene", by.y = "GENE")[5]
  is_small_mutation_heme <- (
    filtered_events$vartype == "SNV" &
      (
        filtered_events$gene %in% hemedb$GENE |
          filtered_events$Tier %in% c(1, 2)
      )
  )
  
  is_other_event = filtered_events$vartype != "SNV"
  # if (any(filtered_events$gene %in% c("FLT3", "DUX4", "KMT2A"))) {
  #   .NotYetImplemented()
  # }
  return(
    filtered_events[is_small_mutation_heme | is_other_event]
  )
}


#' Merge oncokb and multiplicity
#'
#' Finds overlaps between oncokb maf coordinates and multiplicity
#' Oncokb is in MAF coordinates which is based on the altered bases,
#' while multiplicity encodes VCF coordinates based on reference bases.
#'
#' @export
merge_oncokb_multiplicity <- function(oncokb, multiplicity, overwrite = FALSE) {
  is_character_oncokb = is.character(oncokb)
  is_length_one = NROW(oncokb) == 1
  is_exists_oncokb = is_character_oncokb && is_length_one && file.exists(oncokb)
  is_oncokb_rds = is_exists_oncokb && grep("rds$", oncokb)
  is_oncokb_txt = is_exists_oncokb && grep("(c|t)sv$|maf$|txt$", oncokb)

  if (is_oncokb_txt) {
    oncokb <- data.table::fread(oncokb)
  } else if (is_oncokb_rds) {
    oncokb <- readRDS(oncokb)
  }

  is_character_mult = is.character(multiplicity)
  is_length_one = NROW(multiplicity) == 1
  is_exists_mult = is_character_mult && is_length_one && file.exists(multiplicity)
  is_mult_rds = is_exists_mult && grep("rds$", multiplicity)
  is_mult_txt = is_exists_mult && grep("(c|t)sv$|maf$|txt$", multiplicity)

  if (is_mult_rds) {
    multiplicity <- readRDS(multiplicity)
  } else if (is_mult_txt) {
    multiplicity = fread(multiplicity)
  }
  if (inherits(multiplicity, "GRanges")) {
    gr_multiplicity <- multiplicity
  } else {
    multiplicity_tmp = tryCatch(
      {
        dt2gr(multiplicity)
      }, error = function(e) tryCatch(
      {
        as(multiplicity, "GRanges")
      }, error = function(e) NULL)
    )
    if (is.null(multiplicity_tmp)) {
      stop("multiplicity must be coercible to GRanges")
    }
    gr_multiplicity = multiplicity_tmp
  }

  if (inherits(oncokb, "GRanges")) {
    gr_oncokb <- oncokb
    oncokb = gr2dt(oncokb)
  } else {
    oncokb_tmp = tryCatch(
      {
        dt2gr(oncokb)
      }, error = function(e) tryCatch(
      {
        as(oncokb, "GRanges")
      }, error = function(e) NULL)
    )
    if (is.null(oncokb_tmp)) {
      stop("oncokb must be coercible to GRanges")
    }
    gr_oncokb = oncokb_tmp
  }

  gr_oncokb$ALT <- gr_oncokb$Allele

  mc <- S4Vectors::mcols(gr_multiplicity)
  checknormalcols <- c("normal.ref", "normal.alt")
  for (col in checknormalcols) {
    if (!col %in% names(mc)) {
      mc[[col]] <- NA_integer_
    }
  }
  S4Vectors::mcols(gr_multiplicity) <- mc
  ov <- gUtils::gr.findoverlaps(gr_oncokb, gr_multiplicity, by = "ALT", type = "equal")
  ovQuery <- data.table(query.id = integer(0), subject.id = integer(0))
  if (NROW(ov) > 0) {
    ovQuery <- gUtils::gr2dt(ov)[, .(query.id, subject.id)]
  }
  missingIds <- setdiff(1:NROW(gr_oncokb), ov$query.id)

  missingOvQuery <- data.table(query.id = integer(0), subject.id = integer(0))

  if (length(missingIds) > 0) {
    dt_oncokb <- gUtils::gr2dt(gr_oncokb[missingIds])
    invisible(dt_oncokb[, Reference_Allele_Fixed := ifelse(Reference_Allele == "-", "", Reference_Allele)])
    invisible({
      dt_oncokb[Variant_Type == "INS", end := start + nchar(Reference_Allele_Fixed)]
      dt_oncokb[Variant_Type == "DEL", start := end - nchar(Reference_Allele_Fixed)]
    })
    dt_oncokb$oid <- missingIds
    ovMissing <- gUtils::gr.findoverlaps(
      gUtils::dt2gr(dt_oncokb), gr_multiplicity,
      by = "ALT", type = "equal", qcol = c("oid")
    )
    if (NROW(ovMissing) > 0) {
      missingOvQuery <- gr2dt(ovMissing)[, .(query.id = oid, subject.id)]
    }
  }
  ovQuery <- rbind(ovQuery, missingOvQuery)

  # Logic for when there's exact match in coordinates but REF/ALT slight mismatch due to VCF -> maf parsing.
  missingIds = setdiff(1:NROW(gr_oncokb), ovQuery$query.id)
  missingOvQuery = data.table(query.id = integer(0), subject.id = integer(0))
  
  if (length(missingIds) > 0) {
    gr_oncokb_missing = gr_oncokb[missingIds]
    gr_oncokb_missing$oid = missingIds
    ovMissing = gUtils::gr.findoverlaps(
      gr_oncokb_missing, gr_multiplicity
     ,
     type = "equal",
      qcol = c("oid"),
     )
    if (NROW(ovMissing) > 0) {
      missingOvQuery = (
        gr2dt(ovMissing)
        [, .(query.id = oid, subject.id)]
        [!duplicated(query.id)]
      )
    }
  }
  ovQuery = rbind(ovQuery, missingOvQuery)

  # Logic for when there's inexact coordinate match due to VCF -> maf parsing.
  missingIds = setdiff(1:NROW(gr_oncokb), ovQuery$query.id)
  missingOvQuery = data.table(query.id = integer(0), subject.id = integer(0))

  if (length(missingIds) > 0) {
    gr_oncokb_missing = gr_oncokb[missingIds]
    gr_oncokb_missing$oid = missingIds
    ovMissing = gUtils::gr.findoverlaps(
      gr_oncokb_missing, gr_multiplicity
     ,
      qcol = c("oid")
    )
    if (NROW(ovMissing) > 0) {
      missingOvQuery = (
        gr2dt(ovMissing)
        [, .(query.id = oid, subject.id)]
        [!duplicated(query.id)]
      )
    }
  }
  ovQuery = rbind(ovQuery, missingOvQuery)

  cols.keep <- c(
    "ref", "alt", "ref_denoised", "alt_denoised", "normal.ref", "normal.alt",
    "variant.g", "major.count", "minor.count", "major_snv_copies", "minor_snv_copies",
    "total_snv_copies", "total_copies", "VAF", "cn", "altered_copies"
  )
  # if(!any(cols.keep %in% names(mcols(gr_multiplicity)))) {
  #   stop("")
  # }
  cols.keep <- cols.keep[cols.keep %in% names(mcols(gr_multiplicity))]
  subject <- base::subset(gUtils::gr2dt(gr_multiplicity), select = cols.keep)
  skey <- data.table::setkey(ovQuery, query.id)
  skey <- skey[list(1:NROW(oncokb))]

  if (!identical(overwrite, TRUE)) {
    oncokb_multiplicity <- cbind(oncokb, subject[skey$subject.id])
  } else {
    oncokb_multiplicity <- oncokb
    # mc_oncokb_multiplicity = S4Vectors::mcols(oncokb_multiplicity)
    subject_ord <- subject[skey$subject.id]
    for (multiplicity_col in cols.keep) {
      oncokb_multiplicity[[multiplicity_col]] <- NULL
      oncokb_multiplicity[[multiplicity_col]] <- subject_ord[[multiplicity_col]]
    }
    # S4Vectors::mcols(oncokb_multiplicity) = mc_oncokb_multiplicity
  }

  oncokb_multiplicity$multiplicity_id_match <- skey$subject.id
  return(oncokb_multiplicity)
}

