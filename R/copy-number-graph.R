#' @name lift_copy_number_graph
#' @title lift_copy_number_graph
#' @description
#' Create copy number graph JSON files for all samples in a cohort
#'
#' @param cohort Cohort object containing sample information
#' @param output_data_dir Base directory for output files
#' @param is_allelic Boolean indicating whether to process allelic (TRUE) or total (FALSE) copy number
#' @param cores Number of cores for parallel processing (default: 1)
#' @return None
#' @export
lift_copy_number_graph <- function(
    cohort,
    output_data_dir,
    settings = Skilift:::default_settings_path,
    is_allelic = FALSE,
    cores = 1
) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }
    
    if (!dir.exists(output_data_dir)) {
        dir.create(output_data_dir, recursive = TRUE)
    }
    
    # Determine which column to use based on is_allelic
    if (is_allelic) {
        message("Processing allelic copy number")
        cn_column <- "allelic_jabba_gg"
        out_filename <- "allelic.json"
    } else {
        message("Processing total copy number")
        cn_column <- "jabba_gg"
        out_filename <- "complex.json"
    }
    
    # Validate required columns exist
    required_cols <- c("pair", cn_column)
    missing_cols <- required_cols[!required_cols %in% names(cohort$inputs)]
    if (length(missing_cols) > 0) {
        stop("Missing required columns in cohort: ", paste(missing_cols, collapse = ", "))
    }
    
    # Process each sample in parallel
    mclapply(seq_len(nrow(cohort$inputs)), function(i) {
        row <- cohort$inputs[i,]
        pair_dir <- file.path(output_data_dir, row$pair)
        
        if (!dir.exists(pair_dir)) {
            dir.create(pair_dir, recursive = TRUE)
        }
        
        out_file <- file.path(pair_dir, out_filename)
        
        futile.logger::flog.threshold("ERROR")
        tryCatchLog({
            ggraph_path <- row[[cn_column]]
            
            if (!is.null(ggraph_path) && file.exists(ggraph_path)) {
                message(sprintf("Reading gGraph for %s", row$pair))
                ggraph <- readRDS(ggraph_path)
                
                if (!any(class(ggraph) == "gGraph")) {
                    warning(sprintf("Input for %s is not a gGraph object", row$pair))
                    return(NULL)
                }
                
                # Check sequence names overlap with reference
                seq_lengths <- gGnome::parse.js.seqlengths(
                    settings,
                    js.type = "PGV",
                    ref = cohort$reference_name
                )
                
                # Reduce gGraph to only sequences that overlap with reference
                ggraph.reduced <- ggraph[seqnames %in% names(seq_lengths)]
                if (length(ggraph.reduced) == 0) {
                    warning(sprintf(
                        "No overlap between reference sequences and gGraph sequences for %s",
                        row$pair
                    ))
                    return(NULL)
                }
                
                # Set parameters for json export
                mc = mcols(ggraph.reduced$nodes$gr)
                nfields_arg = NULL
                nfields_present = names(mc)[names(mc) %in% c(unlist(row$copy_number_graph_annotations))]
                ggraph.reduced$nodes$mark(Events = NA_character_)
                if (NROW(nfields_present) > 0) {
                    events_string = gGnome:::.dtstring(
                        dt = mc[,nfields_present,drop=FALSE],
                        field_val_sep = ": ",
                        delimiter = "; "
                    )
                    ggraph.reduced$nodes$mark(Events = events_string)
                    nfields_arg = c("col", "Events")
                }   
                
                # if("col" %in% names(mcols(ggraph$nodes$gr))) "col" else NULL
                params <- list(
                    filename = out_file,
                    verbose = TRUE,
                    maxcn = row$copy_number_graph_max_cn,
                    # nfields = if("col" %in% names(mcols(ggraph$nodes$gr))) "col" else NULL,
                    nfields = nfields_arg,
                    annotations = unlist(row$copy_number_graph_annotations),
                    node_field_val_sep = ":",
                    node_delimiter = ";"
                )
                
                # Generate and write JSON
                message(sprintf("Writing copy number graph JSON for %s", row$pair))
                do.call(gGnome::refresh(ggraph.reduced)$json, params)
            } else {
                warning(sprintf("Copy number graph file missing for %s", row$pair))
            }
        }, error = function(e) {
            print(sprintf("Error processing %s: %s", row$pair, e$message))
            NULL
        })
    }, mc.cores = cores, mc.preschedule = TRUE)
    
    invisible(NULL)
}
