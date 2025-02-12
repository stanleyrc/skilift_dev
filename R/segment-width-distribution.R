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





#' Purity Ploidy Plot
#' 
#' Lorem Ipsem
#' 
#' @export
#' @author Aditya Deshpande
create_pp_plot <- function(jabba_gg = NA,
                     het_pileups = NULL,
                     seg.fname = NA,
                     field = 'count',
                     purity = NA,
                     ploidy = NA,
                     plot.min = -2,
                     plot.max = 2,
                     bins = 500,
                     height = 800,
                     width = 800,
                     output.fname = './',
                     is.wgs = FALSE) {
    
    suppressWarnings({
        if (is.na(purity) || is.wgs) {
            if (is.na(jabba_gg) || !file.exists(jabba_gg)) {
                stop("jabba_gg does not exist and no purity and ploidy provided")
            } else {
                jab = readRDS(jabba_gg)
                if (is.null(jab$meta$purity) & is.null(jab$purity)) {
                    stop("jabba_gg does not have purity and ploidy values")
                } else {
                    purity = if (!is.null(jab$meta$purity)) { jab$meta$purity } else { jab$purity }
                    ploidy = if (!is.null(jab$meta$ploidy)) { jab$meta$ploidy } else { jab$ploidy }
                }
            }
        }

        if (is.null(het_pileups) || !file.exists(het_pileups)) {
            stop("het_pileups not supplied")
        }
        hets = Skilift:::grab.hets(het_pileups)
        field = "count"
        if (!field %in% names(values(hets))) {
            stop("hets missing required field")
        }
        if (!field %in% names(mcols(hets))) {
            stop("hets missing required field")
        }
        hets$cn = skitools::rel2abs(hets, field = field, purity = purity, ploidy = ploidy, allele = TRUE) # added allele == T since this are hets
        eqn = skitools::rel2abs(hets, field = field, purity = purity, ploidy = ploidy, allele = TRUE, return.params = TRUE) # Trying to keep it as close to Zi's

        seqlevelsStyle(hets) <- "NCBI"

        # if (is.wgs) {
        #     segs = gr.stripstrand(jab$segstats %Q% (strand(jab$segstats)=="+"))[, c()]
        # } else {
        #     segs = readRDS(seg.fname)
        # }

        segs = jab$nodes$gr[, c()]

        major.segs = gr.val(segs, hets %Q% (allele == "major"), val = "cn", mean = TRUE, na.rm = TRUE)
        minor.segs = gr.val(segs, hets %Q% (allele == "minor"), val = "cn", mean = TRUE, na.rm = TRUE)

        if (is.wgs) {
            tiles = gr.tile(gr = segs, width = 1e4)
            major.tiles = gr.val(tiles, major.segs, val = "cn", mean = TRUE, na.rm = TRUE)
            minor.tiles = gr.val(tiles, minor.segs, val = "cn", mean = TRUE, na.rm = TRUE)
            major.tiles.subs = as.data.table(major.tiles)[sample(.N,500)]
            minor.tiles.subs = as.data.table(minor.tiles)[sample(.N,500)]
            dt = cbind(major.tiles.subs[, .(seqnames, start, end, major.cn = cn)],
                       minor.tiles.subs[, .(minor.cn = cn)])
        } else {
            dt = cbind(as.data.table(major.segs)[, .(seqnames, start, end, major.cn = cn)],
                       as.data.table(minor.segs)[, .(minor.cn = cn)])
        }

        maxval = plot.max * ploidy # max dosage
        minval = plot.min ## min dosage

        dt = dt[major.cn < maxval & minor.cn < maxval &
                major.cn > minval & minor.cn > minval &
                grepl("[0-9]", seqnames)==TRUE,]

        dt[, ratio := major.cn / minor.cn]
        return(
            list(
                pp_plot_data = dt,
                maxval = maxval,
                minval = minval,
                purity = purity,
                ploidy = ploidy,
                eqn = eqn
            )
        )
    })
}


#' Lift 2d purity ploidy plot
#'
#' Lift 2d purity ploidy plot
#'
#' @author Aditya Deshpande
#' @export
lift_pp_plot <- function(cohort, output_data_dir, cores = 1) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }
    
    # Validate required columns exist
    required_cols <- c("pair", "jabba_gg", "het_pileups")

    missing_cols <- required_cols[!required_cols %in% names(cohort$inputs)]
    if (length(missing_cols) > 0) {
        stop("Missing required columns in cohort: ", paste(missing_cols, collapse = ", "))
    }
    
    # Process each sample in parallel
    iterate_function = function(i) {
        row <- cohort$inputs[i,]
        pair_dir <- file.path(output_data_dir, row$pair)
        
        if (!dir.exists(pair_dir)) {
            dir.create(pair_dir, recursive = TRUE)
        }
        
        # out_file <- file.path(pair_dir, "ppfit.json")
        png_path = paste0(normalizePath(pair_dir), "/pp_plot.png")
        
        tryCatch({
            pp_plot_list = create_pp_plot(
                jabba_gg = row$jabba_gg,
                het_pileups = row$het_pileups
            )
            pp_plot_data = pp_plot_list$pp_plot_data
            maxval = pp_plot_list$maxval
            minval = pp_plot_list$minval
            purity = pp_plot_list$purity
            ploidy = pp_plot_list$ploidy
            eqn = pp_plot_list$eqn

            pt = ggplot(pp_plot_data, aes(x = major.cn, y = minor.cn)) +
            scale_x_continuous(breaks = 0:floor(maxval),
                                labels = 0:floor(maxval) %>% as.character,
                                sec.axis = sec_axis(trans = ~(. - eqn["intercept"])/eqn["slope"], # Trying to keep it as close to Zi's
                                                    name = "Major count")) +
            scale_y_continuous(breaks = 0:floor(maxval),
                                labels = 0:floor(maxval) %>% as.character,
                                sec.axis = sec_axis(trans = ~(. - eqn["intercept"])/eqn["slope"], # Trying to keep it as close to Zi's
                                                    name = "Minor count")) +
            labs(x = "Major CN", y = "Minor CN") +
            theme_bw() +
            theme(legend.position = "none",
                    legend.title = element_text(size = 10, family = "sans"),
                    legend.text = element_text(size = 10, family = "sans"),
                    axis.title = element_text(size = 10, family = "sans"),
                    axis.text.x = element_text(size = 10, family = "sans"),
                    axis.text.y = element_text(size = 10, family = "sans")) +
            stat_density_2d(geom = "polygon", contour = TRUE, aes(alpha = 0.5, fill = after_stat(level)),
                            bins = 10) +
            scale_fill_distiller(palette = "Blues", direction = 1) +
            geom_point(size = 2, shape=4, alpha = 0.3) +
            ggtitle(paste0("Purity: ", signif(purity, 2), " Ploidy: ", signif(ploidy, 2)))
            grDevices::png(png_path, units = "in", height = 5, width = 5, res = 600)
            print(pt)
            grDevices::dev.off()
            
        }, error = function(e) {
            warning(sprintf("Error processing %s: %s", row$pair, e$message))
        })
    }

    mclapply(
        seq_len(nrow(cohort$inputs)), 
        iterate_function, 
        mc.cores = cores, 
        mc.preschedule = FALSE
    )
    
    invisible(NULL)
}