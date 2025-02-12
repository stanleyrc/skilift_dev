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
  hemedb_path = "/gpfs/data/imielinskilab/projects/Clinical_NYU/db/master_heme_database.20250128_095937.790322.rds",
  duncavage_path = "~/projects/Clinical_NYU/db/Duncavage_Blood_22.rds"
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
  
  karyotype_string = annotate_karyotype(readRDS(jabba_gg))

  
  hemedb = readRDS(hemedb_path)
  hemedb_guideline = hemedb[, .(GENE, GUIDELINE, DISEASE)][GUIDELINE==TRUE] %>% unique()
  hemedb_guideline = hemedb_guideline[, .(GUIDELINE = GUIDELINE[1], DISEASE = list(DISEASE)), by = GENE]

  ## events_tbl = events_tbl[[1]]

  events_tbl$ix = seq_len(NROW(events_tbl))

  small_muts = events_tbl[vartype == "SNV"]
  other_muts = events_tbl[!vartype == "SNV"]

  small_guidelines = data.table::merge.data.table(small_muts, hemedb_guideline, all.x = TRUE, by.x = "gene", by.y = "GENE")

  is_frequent = small_guidelines$gene %in% hemedb[hemedb$FREQ >= 5]$GENE
  is_guideline = small_guidelines$gene %in% hemedb_guideline$GENE

  small_guidelines = small_guidelines[is_frequent & is_guideline & Tier <= 2]

  small_guidelines$variant = ifelse(is.na(small_guidelines$Variant), small_guidelines$type %>% paste(., "Variant"), small_guidelines$Variant)  %>% tools::toTitleCase()
  small_guidelines$alteration_type = "small"

  small_guidelines$aggregate_label = NA_character_

  small_guidelines$aggregate_label = glue::glue(
    '{small_guidelines[, signif(estimated_altered_copies, 2)]} out of {small_guidelines$segment_cn} copies mutated (VAF: {signif(small_guidelines$VAF, 2)})' 
  )


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


  svs = events_tbl[type == "fusion"]
  fg = svs$fusion_genes
  fg = gsub("@.*$", "", fg)
  
  lst_exons = lapply(strsplit(fg, "::"), function(x) {
    exons = gsub(".*(\\([0-9]+\\)).*", "\\1", x)
    exons = gsub("\\(|\\)", "", exons)
    exons
  })

  svs$lst_exons = lst_exons

  fg = gsub("\\([0-9]+\\)", "", fg)
  svs$fg = fg


  svs_guidelines = data.table::merge.data.table(svs, hemedb_fusions_fr, by.x = "fg", by.y = "Gene")
  exons = do.call(Map, c(f = c, svs_guidelines$lst_exons))

  element_type = dplyr::case_when(
    grepl("UTR", exons[[1]]) & grepl("UTR", exons[[2]]) ~ "UTRs",
    !grepl("UTR", exons[[1]]) & !grepl("UTR", exons[[2]]) ~ "exons",
    TRUE ~ "UTR/exon"
  )

  svs_guidelines$variant = glue::glue('Fusion involving {element_type} {exons[[1]]} <> {exons[[2]]}')
  svs_guidelines$vaf = NA_real_
  svs_guidelines$alteration_type = "rearrangement"
  svs_guidelines$aggregate_label = glue::glue('{svs_guidelines$estimated_altered_copies} fusion cop{ifelse(svs_guidelines$estimated_altered_copies == 1, "y", "ies")}')


  changemap = c(
    "fg" = "gene_name",
    "type" = "variant",
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

  cna_guidelines = events_tbl[type == "SCNA"]

  cna_guidelines$aggregate_label = glue::glue('{cna_guidelines$estimated_altered_copies} out of {cna_guidelines$estimated_altered_copies} copies altered')
  cna_guidelines$alteration_type = "cna"
  cna_guidelines = data.table::merge.data.table(cna_guidelines, hemedb_guideline, by.x = "gene", by.y = "GENE")

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


  allOutputsForJson = rbind(
    smallForJson,
    cnaForJson,
    svsForJson
  )

  jsonlite::write_json(
    list(gene_mutations = allOutputsForJson), 
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
    newnames[newnames %in% old[i]] = new[i]
  }
  names(obj) = newnames
  return(obj)
}

