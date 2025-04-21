#' create_heme_highlights
#'
#' Create highlighted heme events for a single sample. 
#' This code sits downstream of create_filtered_events
#' and is used within the lift_filtered_events method.
#' 
#' @export
create_heme_highlights = function(
  events_tbl, ## filtered events R output
  jabba_gg,
  out_file,
  hemedb_path = Skilift:::HEMEDB(),
  duncavage_path = Skilift:::DUNCAVAGEDB()
) {
  ## {
  ##   "karotype": <string> | null,
  ##   "gene_mutations": [
  ##     {
  ##       "gene_name": <string>,
  ##       "variant_p": <string> | null,
  ##       "vaf": <float> | null,
  ##       "altered_copies": <float> | null,
  ##       "total_copies": <float> | null,
  ##       "alteration_type": <enum("small", "cna", "rearrangement")> | null,
  ##       "aggregate_label": <string> | null
  ##     }
  ##     ] | null
  ##   }
  
  karyotype_string = ""
  if (NROW(jabba_gg) == 1 && is.character(jabba_gg) && file.exists(jabba_gg))
    karyotype_string = annotate_karyotype(readRDS(jabba_gg))

  emptyDfForJson = structure(list(gene_name = character(0), variant = character(0), 
    vaf = numeric(0), tier = integer(0), altered_copies = numeric(0), 
    total_copies = numeric(0), alteration_type = character(0), 
    aggregate_label = character(0), indication = list()), class = c("data.table", 
  "data.frame"), row.names = c(NA, 0L), sorted = "gene_name")

  hemedb = readRDS(hemedb_path)
  hemedb_guideline = hemedb[, .(GENE, GUIDELINE, DISEASE)][GUIDELINE==TRUE] %>% unique()
  hemedb_guideline = hemedb_guideline[, .(GUIDELINE = GUIDELINE[1], DISEASE = list(DISEASE)), by = GENE]

  ## events_tbl = events_tbl[[1]]

  small_muts = events_tbl[events_tbl$vartype == "SNV",]
  smallForJson = data.table::copy(emptyDfForJson)
  criterias = list(
    is_small_in_guidelines = small_muts$gene %in% hemedb_guideline$GENE, ## accounts for merge step
    is_frequent = small_muts$gene %in% hemedb[hemedb$FREQ >= 5]$GENE,
    is_tier2_or_better = small_muts$Tier <= 2
  )
  is_small_mutation_heme_relevant = base::Reduce("&", criterias)
  if (
    NROW(small_muts) > 0
    && any(
      is_small_mutation_heme_relevant
    )
  ) {
    small_guidelines = small_muts[is_small_mutation_heme_relevant]
    ## Inner join
    small_guidelines = data.table::merge.data.table(small_guidelines, hemedb_guideline, all.x = TRUE, by.x = "gene", by.y = "GENE")

    small_guidelines$variant = ifelse(is.na(small_guidelines$Variant), small_guidelines$type %>% paste(., "Variant"), small_guidelines$Variant)  %>% tools::toTitleCase()
    small_guidelines$alteration_type = "small"

    small_guidelines$aggregate_label = NA_character_

    small_guidelines$aggregate_label = glue::glue(
      '{small_guidelines[, signif(estimated_altered_copies, 2)]} out of {small_guidelines$segment_cn} copies mutated (VAF: {signif(small_guidelines$VAF, 2)})' 
    ) %>% as.character()


    changemap = c(
      "gene" = "gene_name",
      "variant" = "variant",
      "VAF" = "vaf",
      "Tier" = "tier",
      "estimated_altered_copies" = "altered_copies",
      "segment_cn" = "total_copies",
      "alteration_type" = "alteration_type",
      "aggregate_label" = "aggregate_label",
      "DISEASE" = "indication"
    )


    smallForJson = base::subset(
      change_names(small_guidelines, changemap),
      select = changemap
    )
  }



  dunc = readRDS(duncavage_path)
  ## base::dput(dunc[Gene == "BCR::ABL1"][1])
  dunc = rbind(
    dunc,
    (
      structure(list(Indication = "MPN", Type = structure(2L, .Label = c("SNV", 
        "Fusions/SV", "SV"), class = "factor"), Gene = "BCR::ABL1", DISEASE = "MPN"), row.names = c(NA, 
          -1L), class = c("data.table", "data.frame"))
      [, DISEASE := "ALL"]
      [, Indication := "ALL"]
    ),
    fill = TRUE
  )
  dunc[, DISEASE := gsub("^B-|^T-", "", Indication)]

  hemedb_fusions = dunc[Type == "Fusions/SV", .(Type = Type[1], DISEASE = list(DISEASE)),  by = .(Gene)]

  hemedb_fusions_rev = data.table::copy(hemedb_fusions)
  fg_split = data.table::tstrsplit(hemedb_fusions_rev$Gene, "::")
  hemedb_fusions_rev$Gene = paste(fg_split[[2]], fg_split[[1]], sep = "::")

  hemedb_fusions_fr = rbind(hemedb_fusions, hemedb_fusions_rev)

  svs = events_tbl[grepl("fusion", events_tbl$type, ignore.case = TRUE)]
  svsForJson = data.table::copy(emptyDfForJson)
  fg_exon = svs$fusion_genes ## this still works if svs is empty
  fg_exon = gsub("@.*$", "", fg_exon)
  fg = gsub("\\([0-9]+\\)", "", fg_exon)
  any_svs_in_guidelines = length(intersect(fg, hemedb_fusions_fr$Gene)) > 0 ## accounts for merge step later
  if (NROW(svs) > 0 && any_svs_in_guidelines) {
    lst_exons = lapply(strsplit(fg_exon, "::"), function(x) {
      exons = gsub(".*(\\([0-9]+\\)).*", "\\1", x)
      exons = gsub("\\(|\\)", "", exons)
      exons
    })

    svs$lst_exons = lst_exons
    svs$fg = fg


    svs_guidelines = data.table::merge.data.table(svs, hemedb_fusions_fr, by.x = "fg", by.y = "Gene")
    
    exons = do.call(Map, c(f = c, svs_guidelines$lst_exons))

    element_type = dplyr::case_when(
      grepl("UTR", exons[[1]]) & grepl("UTR", exons[[2]]) ~ "UTRs",
      !grepl("UTR", exons[[1]]) & !grepl("UTR", exons[[2]]) ~ "exons",
      TRUE ~ "UTR/exon"
    )

    svs_guidelines$variant = glue::glue('Fusion involving {element_type} {exons[[1]]} <> {exons[[2]]}') %>% as.character()
    svs_guidelines$vaf = NA_real_
    svs_guidelines$alteration_type = "rearrangement"
    svs_guidelines$aggregate_label = glue::glue('{svs_guidelines$estimated_altered_copies} fusion cop{ifelse(svs_guidelines$estimated_altered_copies == 1, "y", "ies")}') %>% as.character()


    changemap = c(
      "fg" = "gene_name",
      "variant" = "variant",
      "VAF" = "vaf",
      "Tier" = "tier",
      "estimated_altered_copies" = "altered_copies",
      "segment_cn" = "total_copies",
      "alteration_type" = "alteration_type",
      "aggregate_label" = "aggregate_label",
      "DISEASE" = "indication"
    )

    svsForJson = base::subset(
      change_names(svs_guidelines, changemap),
      select = changemap
    )
  }
  

  cna = events_tbl[grepl("SCNA", events_tbl$type, ignore.case = TRUE),]
  cnaForJson = data.table::copy(emptyDfForJson)
  # any_cna_in_guidelines = length(intersect(cna$gene, hemedb_guideline$GENE)) > 0 ## accounts for merge stelater
  criterias = list(
    is_cna_in_guidelines = cna$gene %in% hemedb_guideline$GENE,
    is_tier2_or_better = cna$Tier <= 2
  )
  is_cna_heme_relevant = Reduce("|", criterias)
  if (NROW(cna) > 0 && any(is_cna_heme_relevant)) {
    cna$aggregate_label = glue::glue('{cna$estimated_altered_copies} out of {cna$estimated_altered_copies} copies altered') %>% as.character()
    cna$alteration_type = "cna"
    cna_guidelines = cna[is_cna_heme_relevant,]
    cna_guidelines = data.table::merge.data.table(cna_guidelines, hemedb_guideline, by.x = "gene", by.y = "GENE", all.x = TRUE)
    cna_guidelines$vartype = tools::toTitleCase(tolower(cna_guidelines$vartype))
    cna_guidelines$DISEASE = lapply(
        cna_guidelines$DISEASE,
        function(x) {
            if (is.null(x)) {
                "MULTIPLE"
            } else {
                x
            }
        }
    )

    changemap = c(
      "gene" = "gene_name",
      "vartype" = "variant",
      "VAF" = "vaf",
      "Tier" = "tier",
      "estimated_altered_copies" = "altered_copies",
      "segment_cn" = "total_copies",
      "alteration_type" = "alteration_type",
      "aggregate_label" = "aggregate_label",
      "DISEASE" = "indication"
    )

    cnaForJson = base::subset(
      change_names(cna_guidelines, changemap),
      select = changemap
    )
  }

  allOutputsForJson = rbind(
    smallForJson,
    cnaForJson,
    svsForJson
  )

  highlights_output = list(
    karyotype = jsonlite::unbox(karyotype_string),
    gene_mutations = allOutputsForJson
  )

  jsonlite::write_json(
    highlights_output, 
    path = out_file,
    na = "null", 
    null = "list",
    pretty = TRUE
  )
  # jsonlite::toJSON(list(gene_mutations = allOutputsForJson), na = "null", null = "list")s

}



#' Change names of data.frame/matrix
#' 
#' Mean to be an easy wrapper that uses either two
#' equal length character vectors, or a single named
#' character vector to rename the columns of a data frame
#' or matrix.
#' 
#' 
change_names = function(obj, old, new) {
  onames = names(obj)
  if (missing(new) && length(names(old)) > 0 && length(old) > 0) {
    new = old
    old = names(old)
  }
  newnames = onames
  if (length(old) != length(new)) stop("lengths of old and new vectors must match")
  for (i in 1:NROW(old)) {
    if (!nzchar(old[i])) next
    if (old[i] == new[i]) next
    newnames[newnames %in% old[i]] = new[i]
  }
  names(obj) = newnames
  return(obj)
}




#' create_summary
#'
#' Create high level summary of
#' mutations that will get passed
#' to metadata.json.
#' 
#' @export
create_summary = function(
  events_tbl, ## filtered events R output
  altered_copies_threshold = 0.9
) {


  small_muts = events_tbl[events_tbl$vartype == "SNV",]
  criterias = list(
    is_tier2_or_better = small_muts$Tier <= 2
   ,
    is_clonal = small_muts$estimated_altered_copies >= altered_copies_threshold
  )
  is_small_mutation_relevant = base::Reduce("&", criterias)
  ## small_muts_parsed = data.table()
  small_muts_parsed = ""
  if (
    NROW(small_muts) > 0
    && any(
      is_small_mutation_relevant
    )
  ) {
    small_muts_out = small_muts[is_small_mutation_relevant]

    small_muts_tally = (
      base::subset(small_muts_out)[, .(gene, type)]
      [, table(type, gene)]
      %>% reshape2::melt()
      %>% as.data.table()
    )
    small_muts_tally = small_muts_tally[value > 0]
    small_muts_parsed = paste(
      small_muts_tally[, paste(type, ": ", paste(gene, collapse = ","), sep = ""), by = type]$V1,
      collapse = "; "
    )
    
  }




  svs = events_tbl[grepl("fusion", events_tbl$type, ignore.case = TRUE)]
  #   svsForJson = data.table::copy(emptyDfForJson)
  fg_exon = svs$fusion_genes ## this still works if svs is empty
  fg_exon = gsub("@.*$", "", fg_exon)
  fg = gsub("\\([0-9]+\\)", "", fg_exon)
  #   any_svs_in_guidelines = length(intersect(fg, hemedb_fusions_fr$Gene)) > 0 ## accounts for merge step later
  criterias = list(
    is_tier2_or_better = svs$Tier <= 2
  )
  is_svs_relevant = Reduce("&", criterias)

  svs_parsed = ""
  if (
    NROW(svs) > 0
    && any(
      is_svs_relevant
    )
  ) {
    svs_parsed = (
      base::subset(svs, is_svs_relevant)[, .(gene, type)]
      [, paste(type, ": ", paste(gene, collapse = ","), sep = "")]
    )
  }
  

  cna = events_tbl[grepl("SCNA", events_tbl$type, ignore.case = TRUE),]
  cna_parsed = ""
  criterias = list(
    # is_cna_in_guidelines = cna$gene %in% hemedb_guideline$GENE,
    is_tier2_or_better = cna$Tier <= 2
  )
  is_cna_relevant = Reduce("&", criterias)
  if (NROW(cna) > 0 && any(is_cna_relevant)) {
    cna_out = cna[is_cna_relevant,]
    cna_out$vartype = tools::toTitleCase(tolower(cna_out$vartype))


    cna_tally = (
      base::subset(cna_out)[, .(gene, vartype)]
      [, table(vartype, gene)]
      %>% reshape2::melt()
      %>% as.data.table()
    )
    cna_tally = cna_tally[value > 0]
    cna_parsed = paste(
      cna_tally[, paste(vartype, ": ", paste(gene, collapse = ","), sep = ""), by = vartype]$V1,
      collapse = "; "
    )
  }

  summary_string = paste(small_muts_parsed, cna_parsed, svs_parsed, sep = "\n")
  summary_string = trimws(summary_string)
  summary_string = gsub("\n{2,}", "\n", summary_string)

  return(summary_string)
  

}

# results <- Filter(Negate(is.null), results)
# 	for (result in results) {
# 	updated_cohort$inputs[result$index, oncotable := result$path][]
# 	}
