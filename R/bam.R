#' @name lift_bam
#' @title lift_bam
#' @description
#' Lift bam file to output data directory
#'
#' @param cohort Cohort object containing sample information
#' @param output_data_dir Base directory for output files
#' @param cores Number of cores for parallel processing (default: 1)
#' @return None
#' @export
lift_bam <- function(cohort, output_data_dir, cores = 1) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }
    
    if (!dir.exists(output_data_dir)) {
        dir.create(output_data_dir, recursive = TRUE)
    }
    
	cohort$inputs
    
    # Check for required column
    if (!"pair" %in% names(cohort$inputs)) {
        stop("Missing required column 'pair' in cohort")
    }
    
	# Checks if cohort needs to be paired
    pairify_cohort_inputs(cohort)
    
    
    # Process each sample in parallel
    mclapply(seq_len(nrow(cohort$inputs)), function(i) {
        row <- cohort$inputs[i,]
        pair_dir <- file.path(output_data_dir, row$pair)
        
        if (!dir.exists(pair_dir)) {
            dir.create(pair_dir, recursive = TRUE)
        }
        
        out_file <- file.path(pair_dir, "metadata.json")
        
        futile.logger::flog.threshold("ERROR")
        tryCatchLog({
           .NotYetImplemented()
        }, error = function(e) {
            print(sprintf("Error processing %s: %s", row$pair, e$message))
            NULL
        })
    }, mc.cores = cores, mc.preschedule = TRUE)
    
    invisible(NULL)
}
