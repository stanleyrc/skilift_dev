#' lift wrapper helper function
#'
#' Helper function to check if required columns exist
has_required_columns <- function(cohort, columns, any = FALSE) {
  if (any) {
    return(any(columns %in% names(cohort$inputs)))
  }
  all(columns %in% names(cohort$inputs))
}


  #' Define required columns for each lifter
  required_columns <- list(
    allelic_copy_number_graph = c("allelic_jabba_gg"),
    total_copy_number_graph = c("events"),
    denoised_coverage = c("tumor_coverage"),
    oncotable = c(
      "somatic_variant_annotations",
      "fusions",
      "jabba_gg",
      "karyograph",
      "events",
      "signature_counts",
      "oncokb_snv",
      "oncokb_cna",
      "oncokb_fusions"
    ),
    filtered_events = c("oncotable", "jabba_gg"),
    aggregated_events = "aggregated_events",
    highlighted_events = "highlighted_events",
    hetsnps = c("het_pileups"),
    germline_multiplicity = c("germline_multiplicity"),
    multiplicity = c("multiplicity"),
    segment_width = c("balanced_jabba_gg", "tumor_coverage"),
    allelic_pp_fit = c("jabba_gg", "het_pileups"),
    multiplicity_fits = c("multiplicity"),
    coverage_jabba_cn = c("jabba_gg", "tumor_coverage"),
    purple_sunrise_plot = c("purple_pp_range", "purple_pp_bestFit"),
    pp_plot = c("jabba_gg", "het_pileups"),
    signatures = c(
      "matrix_sbs_signatures",
      "decomposed_sbs_signatures",
      "matrix_indel_signatures",
      "decomposed_indel_signatures"
    ),
    variant_qc = c("somatic_snvs"),
    metadata = c(
      "tumor_type",
      "disease",
      "primary_site",
      "inferred_sex",
      "jabba_gg",
      "events",
      "somatic_snvs",
      "germline_snvs",
      "tumor_coverage",
      "estimate_library_complexity",
      "alignment_summary_metrics",
      "insert_size_metrics",
      "wgs_metrics",
      "het_pileups",
      "activities_sbs_signatures",
      "activities_indel_signatures",
      "hrdetect",
      "onenesstwoness"
    ),
    bam = c(
      "tumor_bam",
      "normal_bam",
      "bam"
    )
  )


#' Run all lift methods on a cohort
#'
#' This function runs all available lift methods on a cohort to generate the complete
#' set of output files needed for visualization.
#'
#' @param cohort Cohort object containing sample information
#' @param output_data_dir Base directory for output files
#' @param cores Number of cores for parallel processing (default: 1)
#' @param ... Additional arguments
#'
#' @return None
#' @export
lift_all <- function(
    cohort,
    output_data_dir,
    cores = 1,
    ...) {
  # make sure pair is in cohort
  if (!"pair" %in% names(cohort$inputs)) {
    stop("Missing required column in cohort: pair")
  }

  # make sure output_data_dir is passed and exists, create if not
  if (missing(output_data_dir)) {
    stop("output_data_dir is required")
  } else if (!dir.exists(output_data_dir)) {
    dir.create(output_data_dir, recursive = TRUE)
  }

  if (!is.null(cohort$nextflow_results_path)) {
    oncotable_dir <- file.path(cohort$nextflow_results_path, "oncotable")
  } else {
    oncotable_dir <- file.path(output_data_dir, "oncotable")
    message("Oncotable outputs will be placed in: ", oncotable_dir)
  }


  message("Uploading in ", cohort$type, " mode")

  if (cohort$type == "paired") {
    lift_paired(
      cohort = cohort,
      output_data_dir = output_data_dir,
      oncotable_dir = oncotable_dir,
      cores = cores,
      ... = ...
    )
  } else if (cohort$type == "heme") {
    lift_heme(
      cohort = cohort,
      output_data_dir = output_data_dir,
      oncotable_dir = oncotable_dir,
      cores = cores,
      ... = ...
    )
  } else if (cohort$type == "tumor_only") {
    lift_tumor_only(
      cohort = cohort,
      output_data_dir = output_data_dir,
      oncotable_dir = oncotable_dir,
      cores = cores,
      ... = ...
    )
  }

  datafiles_json_path <- file.path(dirname(output_data_dir), "datafiles.json")
  if (!file.exists(datafiles_json_path)) {
    warning("Creating datafiles.json directory")
  }
}

#' Run MVP (Minimum Viable Product) lift methods
#'
#' Helper function containing the common lift methods used across different modes
#'
#' @inheritParams lift_all
#' @param oncotable_dir Directory for oncotable outputs
#'
#' @return Modified cohort object if oncotable is created, otherwise NULL
lift_mvp <- function(
  cohort,
  output_data_dir,
  oncotable_dir,
  cores,
  genome_length = c(1:22, "X", "Y"),
  ...
) {
  list2env(list(...), envir = environment())

  if (has_required_columns(cohort, required_columns$total_copy_number_graph)) {
    lift_copy_number_graph(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      is_allelic = FALSE
    )
  }

  if (has_required_columns(cohort, required_columns$denoised_coverage)) {
    lift_denoised_coverage(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores
    )
  }

  if (has_required_columns(cohort, required_columns$hetsnps)) {
    lift_hetsnps(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores
    )
  }

  if (has_required_columns(cohort, required_columns$filtered_events)) {
    lift_filtered_events(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores
    )
  } else if (has_required_columns(cohort, required_columns$oncotable, any = TRUE) && 'jabba_gg' %in% names(cohort$inputs)) {
    cohort <- create_oncotable(
      cohort = cohort,
      outdir = oncotable_dir,
      cores = cores
    )

    lift_filtered_events(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores
    )
  }

  if (has_required_columns(cohort, required_columns$signatures)) {
    lift_signatures(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores
    )
  }

  if (has_required_columns(cohort, required_columns$multiplicity)) {
    lift_multiplicity(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores
    )
  }

  if (has_required_columns(cohort, required_columns$segment_width)) {
    lift_segment_width_distribution(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores
    )
  }


  if (has_required_columns(cohort, required_columns$multiplicity_fits)) {
    lift_multiplicity_fits(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores
    )
  }

  if (has_required_columns(cohort, required_columns$coverage_jabba_cn)) {
    lift_coverage_jabba_cn(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores
    )
  }

  if (has_required_columns(cohort, required_columns$purple_sunrise_plot)) {
    lift_purple_sunrise_plot(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores
    )
  }

  if (has_required_columns(cohort, required_columns$variant_qc)) {
    lift_variant_qc(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores
    )
  }

  if (has_required_columns(cohort, Skilift:::required_columns$metadata, any = TRUE)) {
    lift_metadata(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      genome_length = genome_length
    )

    lift_datafiles_json(
      output_data_dir = output_data_dir
    )
  }

  if (has_required_columns(cohort, required_columns$pp_plot)) {
    lift_pp_plot(
      cohort = cohort,
      output_data_dir = output_data_dir
    )
  }

  if (has_required_columns(cohort, Skilift:::required_columns$bam, any = TRUE)) {
    lift_bam(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores
    )
  }
  
}

lift_tumor_only <- function(cohort, output_data_dir, oncotable_dir, cores, ...) {
  cohort <- lift_mvp(cohort, output_data_dir, oncotable_dir, cores, ...)
}

lift_heme <- function(cohort, output_data_dir, oncotable_dir, cores, ...) {
  cohort <- lift_mvp(cohort, output_data_dir, oncotable_dir, cores, ...)

  if (has_required_columns(cohort, required_columns$allelic_copy_number_graph)) {
    lift_copy_number_graph(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      is_allelic = TRUE
    )
  }


  if (has_required_columns(cohort, Skilift:::required_columns$multiplicity)) {
    lift_multiplicity(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      is_germline = TRUE
    )
  }

  if (has_required_columns(cohort, required_columns$aggregated_events)) {
    warning("not implemented yet")
  }

  if (has_required_columns(cohort, required_columns$highlighted_events)) {
    warning("not implemented yet")
  }
}

lift_paired <- function(cohort, output_data_dir, oncotable_dir, cores, ...) {
  cohort <- lift_mvp(cohort, output_data_dir, oncotable_dir, cores, ...)

  if (has_required_columns(cohort, required_columns$allelic_copy_number_graph)) {
    lift_copy_number_graph(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      is_allelic = TRUE
    )
  }

  if (has_required_columns(cohort, required_columns$allelic_pp_fit)) {
    lift_allelic_pp_fit(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores
    )
  }

  if (has_required_columns(cohort, required_columns$germline_multiplicity)) {
    lift_multiplicity(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      is_germline = TRUE
    )
  }

}

