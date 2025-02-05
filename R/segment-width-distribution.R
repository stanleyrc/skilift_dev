#' @name create_ppfit_json
#' @title create_ppfit_json
#' @description
#'
#' function to create segmentation plots in case reports
#'
#' @param balanced_jabba_gg the balanced_gg ggraph or the path to it
#' @param tumor_coverage path to associated coverage, if null, will try to pull path from pairs table
#' @param coverage_field field in coverage file to use for segstats
#' @param settings path to settings.json (from PGV)
#' @param ref reference name
#' @param max_na max.na used with JaBbA
#' @param cores cores for JaBbA:::segstats
#' @return segstats table
#' @export
#' @author Stanley Clarke, Tanubrata Dey, Shihab Dider

get_segstats <- function(
  balanced_jabba_gg,
  tumor_coverage,
  coverage_field = "foreground",
  settings = internal_settings_path,
  ref = "hg19",
  max_na = 0.9,
  cores = 1) {

  #' drcln_cov = readRDS(thisp$decomposed_cov[i])
  ## make this work with complex where the cov file was not an input and with jabba_gg
  ## x = path_obj %>% sniff %>% inputs %>% select(CovFile, maxna) #get coverage that was used for the jabba run
  if (is.null(balanced_jabba_gg)) {
      stop("Please provide a valid path to a non-integer balanced gGraph file.")
  }

  if (!is.null(tumor_coverage)) {
      cov <- readRDS(tumor_coverage)
  } else {
      stop("Please provide a valid path to a coverage file.")
  }

  ## need to replace NaN with NA or JaBbA:::segstats breaks
  mcols(cov)[[coverage_field]] <- ifelse(is.nan(mcols(cov)[[coverage_field]]), 
                                        NA, 
                                        mcols(cov)[[coverage_field]])
  mcols(cov)[[coverage_field]] <- as.numeric(mcols(cov)[[coverage_field]])

  balanced_gg_gr <- readRDS(balanced_jabba_gg)$nodes$gr

  segstats <- JaBbA:::segstats(
    balanced_gg_gr,
    cov,
    field = coverage_field,
    prior_weight = 1,
    max.chunk = 1e8,
    ## subsample = subsample,
    mc.cores = cores,
    verbose = FALSE,
    max.na = max_na,
    lp = FALSE
  )
  segstats_dt <- gr2dt(segstats)
  names(segstats_dt) <- gsub("\\.", "_", names(segstats_dt))
  return(segstats_dt)
}

#' @name lift_segment_width_distribution
#' @title lift_segment_width_distribution
#' @description
#' Create segment width distribution JSON files for all samples in a cohort
#'
#' @param cohort Cohort object containing sample information
#' @param output_data_dir Base directory for output files
#' @param annotations (list) Optional annotations to include in JSON output (default: NULL)
#' @param cores Number of cores for parallel processing (default: 1)
#' @return None
#' @export
lift_segment_width_distribution <- function(cohort, output_data_dir, annotations = NULL, cores = 1) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }
    
    # Validate required columns exist
    required_cols <- c("pair", "balanced_jabba_gg", "tumor_coverage")
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
        
        out_file <- file.path(pair_dir, "ppfit.json")
        
        tryCatch({
            if (!file.exists(row$balanced_jabba_gg)) {
                warning(sprintf("Balanced JaBbA file not found for %s: %s", row$pair, row$balanced_jabba_gg))
                return(NULL)
            }
            
            # Read the gGraph
            ggraph <- readRDS(row$balanced_jabba_gg)
            if (!any(class(ggraph) == "gGraph")) {
                warning(sprintf("File is not a gGraph for %s: %s", row$pair, row$balanced_jabba_gg))
                return(NULL)
            }
            
            # Get sequence lengths from the gGraph
            seq_lengths <- seqlengths(ggraph$nodes$gr)
            
            # Check for required fields
            colnames_check <- c("start_ix", "end_ix", "eslack_in", "eslack_out", 
                              "edges_in", "edges_out", "tile_id", "snode_id", 
                              "loose_left", "loose_right", "loose_cn_left", 
                              "loose_cn_right", "node_id", "raw_mean", "raw_var", 
                              "nbins", "nbins_tot", "nbins_nafrac", "wbins_nafrac", 
                              "wbins_ok", "mean", "bad", "max_na", "loess_var", 
                              "tau_sq_post", "post_var", "var", "sd")
            
            if (all(colnames_check %in% names(ggraph$nodes$dt))) {
                gg_w_segstats <- ggraph
                fields.keep <- c(colnames_check, "seqnames", "start", "end", "strand", "width", "loose", "index")
            } else {
                # Get segstats information
                segstats.dt <- get_segstats(
                    balanced_jabba_gg = row$balanced_jabba_gg,
                    tumor_coverage = row$tumor_coverage
                )
                segstats.gr <- GRanges(segstats.dt, seqlengths = seq_lengths) %>% trim()
                gg_w_segstats <- gG(nodes = segstats.gr, edges = ggraph$edges$dt)
                fields.keep <- names(segstats.dt) %>% grep("cn", ., invert = TRUE, value = TRUE)
            }
            
            gg_w_segstats$set(y.field = "cn")
            
            # Check for sequence name overlap
            ggraph.reduced <- gg_w_segstats[seqnames %in% names(seq_lengths)]
            if (length(ggraph.reduced) == 0) {
                warning(sprintf("No overlap in sequence names for %s", row$pair))
                return(NULL)
            }
            
            # Create JSON
            gGnome::refresh(ggraph.reduced)$json(
                filename = out_file,
                verbose = TRUE,
                annotations = if (!is.null(annotations)) unlist(annotations) else NULL,
                maxcn = 500,
                nfields = fields.keep,
                save = TRUE
            )
            
        }, error = function(e) {
            warning(sprintf("Error processing %s: %s", row$pair, e$message))
        })
    }, mc.cores = cores, mc.preschedule = FALSE)
    
    invisible(NULL)
}
