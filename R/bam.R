#' @name lift_bam
#' @title lift_bam
#' @description
#' Lift bam file to output data directory.
#' Currently hardlinking.
#'
#' @param cohort Cohort object containing sample information
#' @param output_data_dir Base directory for output files
#' @param cores Number of cores for parallel processing (default: 1)
#' @return None
#' @export
lift_bam <- function(cohort, output_data_dir, cores = 1, overwrite = FALSE) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }
    
    if (!dir.exists(output_data_dir)) {
        dir.create(output_data_dir, recursive = TRUE)
    }
    
    # Check for required column
    if (!"pair" %in% names(cohort$inputs)) {
        stop("Missing required column 'pair' in cohort")
    }
    
    # Checks if cohort needs to be paired
    # Ensures that "tumor_bam" and "normal_bam" columns are present
    lift_inputs = pairify_cohort_inputs(cohort, tumor_normal_columns = c("bam"), sep = "_")

    if (!any(c("tumor_bam", "normal_bam") %in% names(lift_inputs))) {
        message("No tumor_bam or normal_bam columns present in cohort inputs")
        return(NULL)
    }
    if (!is.null(lift_inputs$tumor_bam)) {
      bambai_path = paste(lift_inputs$tumor_bam, ".bai", sep = "")
      bai_path = paste(tools::file_path_sans_ext(lift_inputs$tumor_bam), ".bai", sep = "")
      lift_inputs$tumor_bam_bai = dplyr::case_when(
                                           file.exists(bambai_path) ~ bambai_path,
                                           file.exists(bai_path) ~ bai_path,
                                           TRUE ~ NA_character_
                                         )
    }
    if (!is.null(lift_inputs$normal_bam)) {
      bambai_path = paste(lift_inputs$normal_bam, ".bai", sep = "")
      bai_path = paste(tools::file_path_sans_ext(lift_inputs$normal_bam), ".bai", sep = "")
      lift_inputs$normal_bam_bai = dplyr::case_when(
                                           file.exists(bambai_path) ~ bambai_path,
                                           file.exists(bai_path) ~ bai_path,
                                           TRUE ~ NA_character_
                                         )
    }

    # Process each sample in parallel
    force_flag = if (overwrite) "f" else ""
    upload_bam = function(i) {
        row <- lift_inputs[i,]
        pair_dir <- file.path(output_data_dir, row$pair)
        
        if (!dir.exists(pair_dir)) {
            dir.create(pair_dir, recursive = TRUE)
        }
        
        tumor_bam_file <- file.path(pair_dir, "tumor.bam")
        tumor_bai_file <- file.path(pair_dir, "tumor.bam.bai")

        normal_bam_file <- file.path(pair_dir, "normal.bam")
        normal_bai_file <- file.path(pair_dir, "normal.bam.bai")
        
        futile.logger::flog.threshold("ERROR")
        tryCatchLog({
          exit_tbam = 0
          exit_tbai = 0
          exit_nbam = 0
          exit_nbai = 0
          if (!is.null(row$tumor_bam)) {
            cmd_bam = glue::glue('ln -n{force_flag} {row$tumor_bam} {tumor_bam_file}')
            cmd_bai = glue::glue('ln -n{force_flag} {row$tumor_bam_bai} {tumor_bai_file}')
            exit_tbam = system(cmd_bam)
            exit_tbai = system(cmd_bai)
          }
          stopifnot(exit_tbam == 0)
          stopifnot(exit_tbai == 0)
          if (!is.null(row$normal_bam)) {
            cmd_bam = glue::glue('ln -n{force_flag} {row$normal_bam} {normal_bam_file}')
            cmd_bai = glue::glue('ln -n{force_flag} {row$normal_bam_bai} {normal_bai_file}')
            exit_nbam = system(cmd_bam)
            exit_nbai = system(cmd_bai)
          }
          stopifnot(exit_nbam == 0)
          stopifnot(exit_nbai == 0)
 
        }, error = function(e) {
            print(sprintf("Error processing %s: %s", row$pair, e$message))
            NULL
        })
    }
    mclapply(seq_len(nrow(lift_inputs)), upload_bam, mc.cores = cores, mc.preschedule = TRUE)
    
    invisible(NULL)
}
