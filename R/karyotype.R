#' SV in GRangesList to bedpe
#' 
#' Convert GRangesList representation of SVs to bedpe-like format
#' 
#' @export
grl2bedpe = function (
	grl, 
	add_breakend_mcol = FALSE, 
	flip = FALSE, 
	as.data.table = TRUE, 
    zerobased = TRUE
) {
    grpiv = gUtils::grl.pivot(grl)
    if (zerobased) {
      GenomicRanges::start(grpiv[[1]]) = GenomicRanges::start(grpiv[[1]]) - 1
      GenomicRanges::start(grpiv[[2]]) = GenomicRanges::start(grpiv[[2]]) - 1
    }
    mcgrl = as.data.frame(mcols(grl))
    df1 = as.data.frame(grpiv[[1]])[, c(1:3, 5), drop = F]
    df2 = as.data.frame(grpiv[[2]])[, c(1:3, 5), drop = F]
    colnames(df1) = c("chrom1", "start1", "end1", "strand1")
    colnames(df2) = c("chrom2", "start2", "end2", "strand2")
    mc1 = data.frame()[seq_len(NROW(df1)), , drop = F]
    mc2 = data.frame()[seq_len(NROW(df2)), , drop = F]
    if (isTRUE(add_breakend_mcol)) {
        mc1 = as.data.frame(grpiv[[1]])[, -c(1:5), drop = F]
        mc2 = as.data.frame(grpiv[[2]])[, -c(1:5), drop = F]
        colnames(mc1) = paste0("first.", colnames(mc1))
        colnames(mc2) = paste0("second.", colnames(mc2))
    }
    out = cbind(df1, df2, name = as.character(seq_len(NROW(df1))), 
        score = rep_len(0, NROW(df1)), mcgrl, mc1, mc2)
    canon_col = c(1, 2, 3, 5, 6, 7, 9, 10, 4, 8)
    nix = seq_len(ncol(out))
    out = out[, c(canon_col, nix[!nix %in% canon_col]), drop = F]
    if (flip) {
        out$strand1 = c(`+` = "-", `-` = "+")[out$strand1]
        out$strand2 = c(`+` = "-", `-` = "+")[out$strand2]
    }
    if (as.data.table) 
        return(data.table::as.data.table(out))
    else return(out)
}

#' Annotate Karyotype
#' 
#' Convert gGraph to karyotype
#' 
#' @export
annotate_karyotype = function(
	gg, 
	segment_size_threshold = 5e6, 
	chr_fraction_threshold = 0.9, 
	arm_fraction_threshold = 0.9, 
	band_segment_ratio_threshold = 0.1
) {
  cyto = Skilift:::process_cytoband()
  cyto$arm = gsub("^(p|q).*", "\\1", cyto$band)
  cyto$band_width = as.integer(width(cyto))

  cyto_by_arm = gGnome::gr_construct_by(cyto, "arm")
  rcyto_by_arm = GenomicRanges::reduce(cyto_by_arm)
  rcyto_arm = gGnome::gr_deconstruct_by(rcyto_by_arm, by = "arm", meta = TRUE)
  rcyto_arm_dt = gUtils::gr2dt(rcyto_arm)

  annotated_tra_inv = ""

  ge = gg$edges
  geAlt = ge[ge$dt$type == "ALT"]
  if (length(geAlt) > 0) {
    sortedGrl = BiocGenerics::sort(GenomeInfoDb::sortSeqlevels(geAlt$grl), ignore.strand = TRUE)

    widths = width(range(sortedGrl, ignore.strand = TRUE))
    widths[base::lengths(widths) == 2] = list(NA_integer_)
    mcols(sortedGrl)$widths = unlist(widths)
    breakpoints = gUtils::grl.unlist(sortedGrl)
    breakpoints = gUtils::gr.end(breakpoints, ignore.strand = FALSE)
    breakpoints = breakpoints %*% cyto

    grl = S4Vectors::split(breakpoints, breakpoints$grl.ix)
	bedpe = grl2bedpe(grl, add_breakend_mcol = TRUE, zerobased = FALSE)

    isTranslocation = !(bedpe$first.chrom_name == bedpe$second.chrom_name)
    k_type = dplyr::case_when(
      bedpe$first.class == "TRA-like" ~ "t",
      bedpe$first.class == "INV-like" ~ "inv",
      TRUE ~ NA_character_
    )
    k_chrom = ifelse(
      bedpe$first.chrom_name == bedpe$second.chrom_name,
      bedpe$first.chrom_name,
      paste(bedpe$first.chrom_name, bedpe$second.chrom_name, sep = ";")
    )
    isBig = (
		(
			isTranslocation | 
			(bedpe$first.widths > segment_size_threshold)
		)
		& !is.na(k_type)
	)
    k_bands = paste(bedpe$first.band, bedpe$second.band, sep = ";")

    if (any(isBig))
      annotated_tra_inv = glue::glue('{k_type}({k_chrom})({k_bands})')[isBig]
  }

  segments = gg$nodes$gr
  segments = BiocGenerics::sort(
	GenomeInfoDb::sortSeqlevels(segments), 
	ignore.strand = TRUE
  )
  segments = GenomeInfoDb::keepStandardChromosomes(
	segments, 
	pruning.mode = "coarse"
  )
  genome_sl = GenomeInfoDb::seqlengths(segments)
  genome_sl = data.table::data.table(
    seqnames = names(genome_sl),
    chr_length = genome_sl
  )

  annotated_cna = ""

  events = segments[
    segments$cn != 2 
    ## FIXME: account properly for Y chromosomes
    & seqnames(segments) != "Y"
  ]

  if (length(events) > 0) {
    events$label = ifelse(events$cn > 2, "dup", "del")

    events_by_label = gGnome::gr_construct_by(events, "label")
    revents_by_label = GenomicRanges::reduce(events_by_label, ignore.strand = TRUE)
    revents = gGnome::gr_deconstruct_by(revents_by_label, by = "label", meta = TRUE)
    revents = revents %*% cyto

    revents_dt = gUtils::gr2dt(revents)

    revents_dt[, segment_arm_id := .GRP, by = .(query.id, arm)]

    revents_dt = data.table::merge.data.table(revents_dt, genome_sl, by = "seqnames")

    chr_stats = revents_dt[, sum(width) / chr_length[1], by = .(seqnames, label)]

    chr_events = data.table::copy(chr_stats[V1 > chr_fraction_threshold])
    chr_events[, label := ifelse(label == "dup", "+", "-")]
    annotated_chr_cna = with(chr_events, glue::glue('{label}{seqnames}'))

    revents_arm_dt = data.table::merge.data.table(revents_dt, rcyto_arm_dt[, .(seqnames, arm, arm_width = width)], by = c("seqnames", "arm"))
    arm_stats = revents_arm_dt[, sum(width) / arm_width[1], by = .(seqnames, label, arm)]
    arm_events = arm_stats[V1 > arm_fraction_threshold & !seqnames %in% chr_events$seqnames]

    annotated_arm_cna = with(arm_events, glue::glue('{label}({seqnames})({arm})'))

    band_events_prefiltered = data.table::copy(revents_dt)
    band_events_prefiltered = band_events_prefiltered[!seqnames %in% chr_events$seqnames]
    band_events_prefiltered = band_events_prefiltered[!paste(seqnames, arm) %in% arm_events[, paste(seqnames, arm)]]
    band_events_prefiltered = band_events_prefiltered[stain != "acen"]

    band_events_prefiltered[, band_segment_ratio := sum(width) / band_width[1], by = .(segment_arm_id)]
    band_events_prefiltered[, total_segment_width := sum(width), by = segment_arm_id]

    band_events = band_events_prefiltered[total_segment_width > segment_size_threshold & band_segment_ratio >= band_segment_ratio_threshold]

    band_events_collapsed = band_events[, .(
      seqnames = seqnames[1],
      min_band = band[1],
      max_band = band[.N],
      arm = arm[1],
      label = label[1]
    ), by = .(segment_arm_id)]

    band_events_collapsed[, first_band := ifelse(arm == "p", max_band, min_band)]
    band_events_collapsed[, last_band := ifelse(arm == "p", min_band, max_band)]

    annotated_band_cna = with(band_events_collapsed, {
      glue::glue('{label}({seqnames})({first_band},{last_band})')
    })

  }
  
  karyotype_string = paste(c(
    annotated_tra_inv,
    annotated_chr_cna,
    annotated_arm_cna,
    annotated_band_cna
  ), collapse = ",")

  karyotype_string = gsub("^,|,$", "", karyotype_string)
  
  return(karyotype_string)
}