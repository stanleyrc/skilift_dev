# #' @name create_highlights
# #' @title create_highlights
# #' @description
# #' Create highlighted events for a single sample
# #'
# #' @param pair patient id
# #' @param oncotable oncotable task output
# #' @param jabba_gg JaBbA output ggraph or complex
# #' @param out_file path to write json
# #' @param temp_fix TRUE/FALSE whether to apply temporary fix
# #' @param return_table TRUE/FALSE whether to return the data.table
# #' @return data.table or NULL
# #' @export
# create_highlights <- function(
#     pair,
#     oncotable,
#     jabba_gg,
#     out_file,
#     temp_fix = FALSE,
#     return_table = FALSE,
#     type = "paired") {

#     ot <- readRDS(oncotable)

#     # add a fusion_gene_coords column of NAs if no fusions
#     if (!"fusion_gene_coords" %in% colnames(ot)) {
#         ot[, fusion_genes := NA]
#         ot[, fusion_gene_coords := NA]
#     }
#     # snvs <- ot[grepl("frameshift|missense|stop|disruptive", annotation, perl = TRUE)]
#     snvs <- ot[vartype == "SNV"][!is.na(type)]
#     if ("tier" %in% colnames(ot)) {
#       snvs <- snvs[order(is.na(tier)), ]
#     }
#     snvs[, is_unique_p := !is.na(variant.p) & !duplicated(cbind(gene, variant.p))]
#     snvs[, is_unique_g := !duplicated(cbind(gene, variant.g))]
#     snvs[, is_unique_c := !duplicated(cbind(gene, variant.c))]
#     snvs <-  snvs[is_unique_p | (is_unique_g & is_unique_c)]
#     snvs$is_unique_p = NULL
#     snvs$is_unique_g = NULL
#     snvs$is_unique_c = NULL
#     homdels <- ot[type == "homdel"][, vartype := "HOMDEL"][, type := "SCNA"]
#   	amps <- ot[type == "amp"][, vartype := "AMP"][, type := "SCNA"]
#     fusions <- ot[type == "fusion"]
#     possible_drivers <- rbind(snvs, homdels, amps, fusions)

#     oncotable_col_to_filtered_events_col <- c(
#         "id" = "id",
#         "gene" = "gene",
#         "fusion_genes" = "fusion_genes",
#         "fusion_gene_coords" = "fusion_gene_coords",
#         "value" = "fusion_cn",
#         "vartype" = "vartype",
#         "type" = "type",
#         "variant.g" = "Variant_g",
#         "variant.p" = "Variant",
#         "total_copies" = "estimated_altered_copies",
#         "segment_cn" = "segment_cn",
#         "ref" = "ref",
#         "alt" = "alt",
#         "VAF" = "VAF",
#         "gene_location" = "Genome_Location",
#         "tier" = "Tier",
#         "role" = "role",
#         "gene_summary" = "gene_summary",
#         "variant_summary" = "variant_summary",
#         "effect" = "effect",
#         "effect_description" = "effect_description",
#         "therapeutics" = "therapeutics",
#         "resistances" = "resistances",
#         "diagnoses" = "diagnoses",
#         "prognoses" = "prognoses"
#     )
#     filtered_events_columns <- names(possible_drivers)[names(possible_drivers) %in% names(oncotable_col_to_filtered_events_col)]
    
#     res <- possible_drivers[, ..filtered_events_columns]
#     intersected_columns <- intersect(filtered_events_columns, names(res))
#     setnames(res, old = intersected_columns, new = oncotable_col_to_filtered_events_col[intersected_columns])

#   res <- res %>% unique(., by = c("gene","vartype", "Variant"))
#   if (nrow(res) > 0) {
#       res[, seqnames := tstrsplit(Genome_Location, ":", fixed = TRUE, keep = 1)]
#       res[, start := tstrsplit(Genome_Location, "-", fixed = TRUE, keep = 1)]
#       res[, start := tstrsplit(start, ":", fixed = TRUE, keep = 2)]
#       res[, end := tstrsplit(Genome_Location, "-", fixed = TRUE, keep = 2)]
#       res.mut <- res[!is.na(Variant), ]
#       if (nrow(res.mut) > 0) {
#           #res.mut[, Variant := gsub("p.", "", Variant)]
#           res.mut[, vartype := "SNV"]
#           # TODO:
#           # truncating mutations are not always deletions
#           # initial logic may be misleading calling all small mutations "SNV"
#           # but we should encode this as something more robust
#           # res.mut[type=="trunc", vartype := "DEL"]
#       }
#       res.cn <- res[is.na(Variant) & !is.na(Genome_Location), ]
#       if (nrow(res.cn) > 0) {
#           jab <- readRDS(jabba_gg)
#           res.cn.gr <- GRanges(res.cn)
#           res.cn.gr <- gr.val(res.cn.gr, jab$nodes$gr, c("cn", "cn.low", "cn.high"))
#           res.cn.dt <- as.data.table(res.cn.gr)
#           res.cn.dt[, estimated_altered_copies := abs(cn - 2)]
#           res.cn.dt[, segment_cn := cn]
#           res.cn.dt[, Variant := vartype]
#           ## res.cn.dt[!is.na(cn) & !is.na(cn.low) & !is.na(cn.high), Variant := paste0("Total CN:", round(cn, digits = 3), "; CN Minor:", round(cn.low, digits = 3), "; CN Major:", round(cn.high, digits = 3))]
#           ## res.cn.dt[!is.na(cn) & is.na(cn.low) & is.na(cn.high), Variant := paste0("Total CN:", round(cn, digits = 3)                                                                                     )]
#           if (temp_fix) {
#               res.cn.dt <- res.cn.dt[!(type == "homdel" & cn != 0), ]
#               res.cn.dt <- res.cn.dt[!(type == "amp" & cn <= 2), ]
#           }
#           res.cn.dt[, c("cn", "cn.high", "cn.low", "width", "strand") := NULL] # make null, already added to Variant
#           res.final <- rbind(res.mut, res.cn.dt, fill = TRUE)
#       } else {
#           res.final <- res.mut
#       }
#       if (type == "heme") {
#         res.final = select_heme_events(res.final)
#       }
#       write_json(res.final, out_file, pretty = TRUE)
#       res.final[, sample := pair]
#       if (return_table) {
#           return(res.final)
#       }
#     }
# }

# #' @name lift_highlights
# #' @title lift_highlights
# #' @description
# #' Create highlights events for all samples in a cohort
# #'
# #' @param cohort Cohort object containing sample information
# #' @param output_data_dir Base directory for output files
# #' @param cores Number of cores for parallel processing (default: 1)
# #' @return None
# #' @export
# lift_highlights <- function(cohort, output_data_dir, cores = 1) {
#     if (!inherits(cohort, "Cohort")) {
#         stop("Input must be a Cohort object")
#     }
    
#     if (!dir.exists(output_data_dir)) {
#         dir.create(output_data_dir, recursive = TRUE)
#     }
    
#     # Validate required columns exist
#     required_cols <- c("pair", "oncotable", "jabba_gg")
#     missing_cols <- required_cols[!required_cols %in% names(cohort$inputs)]
#     if (length(missing_cols) > 0) {
#         stop("Missing required columns in cohort: ", paste(missing_cols, collapse = ", "))
#     }
    
#     cohort_type = cohort$type
#     # Process each sample in parallel
#     mclapply(seq_len(nrow(cohort$inputs)), function(i) {
#         row <- cohort$inputs[i,]
#         pair_dir <- file.path(output_data_dir, row$pair)
        
#         if (!dir.exists(pair_dir)) {
#             dir.create(pair_dir, recursive = TRUE)
#         }
        
#         out_file <- file.path(pair_dir, "filtered.events.json")
        
#         tryCatch({
#             create_filtered_events(
#                 pair = row$pair,
#                 oncotable = row$oncotable,
#                 jabba_gg = row$jabba_gg,
#                 out_file = out_file,
#                 temp_fix = FALSE,
#                 return_table = FALSE,
#                 cohort_type = cohort_type
#             )
            
#         }, error = function(e) {
#             warning(sprintf("Error processing %s: %s", row$pair, e$message))
#         })
#     }, mc.cores = cores, mc.preschedule = FALSE)
    
#     invisible(NULL)
# }
