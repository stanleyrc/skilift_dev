  #' lift wrapper helper function
  #' 
  #' Helper function to check if required columns exist
  has_required_columns <- function(cohort, columns, any=FALSE) {
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
#' @param settings Settings path for copy number graph
#' @param is_allelic Boolean for copy number graph processing (default: FALSE)
#' @param max_cn Maximum copy number for graph (default: 100)
#' @param annotations List of annotations for copy number graph and segment width distribution
#' @param cohort_column Field in cohort inputs table for coverage track (default: "tumor_coverage")
#' @param coverage_field Field specifying values to write to arrow (default: "foreground")
#' @param color_field Field specifying color values
#' @param bin.width Bin width for coverage track
#' @param is_germline Logical for multiplicity (default: FALSE)
#' @param node_metadata Additional columns for multiplicity
#' @param field Field for multiplicity (default: "total_copies")
#' @param genome_length Chromosome names for metadata
#' @param ... Additional arguments
#'
#' @return None
#' @export
lift_all <- function(
  cohort,
  output_data_dir,
  cores = 1,
  settings = Skilift:::internal_settings_path,
  max_cn = 100,
  annotations = NULL,
  coverage_field = "foreground",
  color_field = NULL,
  bin.width = 1e4,
  node_metadata = c("gene", "feature_type", "annotation", "REF", "ALT", "variant.c", "variant.p", "vaf", "transcript_type", "impact", "rank"),
  multiplicity_field = "altered_copies",
  genome_length = c(1:22, "X", "Y"),
  path_to_nodejs = "/gpfs/share/apps/nodejs/22.9.0/bin/node",
  ...
) {
  
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
    oncotable_dir = file.path(cohort$nextflow_results_path, "oncotable")
  } else {
    oncotable_dir = file.path(output_data_dir, "oncotable")
    message("Oncotable outputs will be placed in: ", oncotable_dir)
  }
  

  message("Uploading in ", cohort$cohort_type, " mode")

  if (cohort$cohort_type == "paired") {
    lift_paired(
      cohort = cohort, 
      output_data_dir = output_data_dir,  
      oncotable_dir = oncotable_dir,
      cores = cores,
      settings = settings,
      max_cn = max_cn,
      annotations = annotations,
      coverage_field = coverage_field,
      color_field = color_field,
      bin.width = bin.width,
      node_metadata = node_metadata,
      multiplicity_field = multiplicity_field,
      genome_length = genome_length,
      ... = ...
    )
  } else if (cohort$cohort_type == "heme") {
    lift_heme(
      cohort = cohort, 
      output_data_dir = output_data_dir,  
      oncotable_dir = oncotable_dir,
      cores = cores,
      settings = settings,
      max_cn = max_cn,
      annotations = annotations,
      coverage_field = coverage_field,
      color_field = color_field,
      bin.width = bin.width,
      node_metadata = node_metadata,
      multiplicity_field = multiplicity_field,
      genome_length = genome_length,
      ... = ...
    )
  } else if (cohort$cohort_type == "tumor_only") {
    lift_tumor_only(
      cohort = cohort, 
      output_data_dir = output_data_dir,  
      oncotable_dir = oncotable_dir,
      cores = cores,
      settings = settings,
      max_cn = max_cn,
      annotations = annotations,
      coverage_field = coverage_field,
      color_field = color_field,
      bin.width = bin.width,
      node_metadata = node_metadata,
      multiplicity_field = multiplicity_field,
      genome_length = genome_length,
      ... = ...
    )
  }

  datafiles_json_path = file.path(dirname(output_data_dir), "datafiles.json")
  if (!file.exists(datafiles_json_path)) {
    warning("Creating datafiles.json directory")
  }
}



lift_tumor_only = function(cohort, output_data_dir, ...) {

  if (has_required_columns(cohort, Skilift:::required_columns$allelic_copy_number_graph)) {
    lift_copy_number_graph(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      settings = settings,
      is_allelic = TRUE,
      max_cn = max_cn,
      annotations = annotations
    )
  }
  
  if (has_required_columns(cohort, Skilift:::required_columns$total_copy_number_graph)) {
    lift_copy_number_graph(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      settings = settings,
      is_allelic = FALSE,
      max_cn = max_cn,
      annotations = annotations
    )
  }

  if (has_required_columns(cohort, Skilift:::required_columns$denoised_coverage)) {
    lift_denoised_coverage(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      coverage_field = coverage_field,
      color_field = color_field,
      bin.width = bin.width
    )
  }

  if (has_required_columns(cohort, Skilift:::required_columns$hetsnps)) {
    lift_hetsnps(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores
    )
  }



  if (has_required_columns(cohort, Skilift:::required_columns$filtered_events)) {
    lift_filtered_events(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores
    )
  # oncotable doesn't need every single column, just any one of them
  } else if (has_required_columns(cohort, Skilift:::required_columns$oncotable, any = TRUE)) {
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

  if (has_required_columns(cohort, Skilift:::required_columns$multiplicity)) {
    lift_multiplicity(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      is_germline = FALSE,
      node_metadata = node_metadata,
      field = multiplicity_field
    )
  }

  if (has_required_columns(cohort, Skilift:::required_columns$germline_multiplicity)) {
    lift_multiplicity(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      is_germline = TRUE,
      node_metadata = node_metadata,
      field = multiplicity_field
    )
  }

  if (has_required_columns(cohort, Skilift:::required_columns$segment_width)) {
    lift_segment_width_distribution(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      annotations = annotations
    )
  }


  if (has_required_columns(cohort, Skilift:::required_columns$signatures)) {
    lift_signatures(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores
    )
  }

  if (has_required_columns(cohort, Skilift:::required_columns$variant_qc)) {
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
  }

  if (has_required_columns(cohort, Skilift:::required_columns$pp_plot)) {
    lift_pp_plot(
      cohort = cohort, 
      output_data_dir = output_data_dir
    )
  }
  
}

lift_heme = function(cohort, output_data_dir, ...) {
  list2env(list(...), envir = environment())

  if (has_required_columns(cohort, Skilift:::required_columns$allelic_copy_number_graph)) {
    lift_copy_number_graph(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      settings = settings,
      is_allelic = TRUE,
      max_cn = max_cn,
      annotations = annotations
    )
  }
  
  if (has_required_columns(cohort, Skilift:::required_columns$total_copy_number_graph)) {
    lift_copy_number_graph(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      settings = settings,
      is_allelic = FALSE,
      max_cn = max_cn,
      annotations = annotations
    )
  }

  if (has_required_columns(cohort, Skilift:::required_columns$denoised_coverage)) {
    lift_denoised_coverage(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      coverage_field = coverage_field,
      color_field = color_field,
      bin.width = bin.width
    )
  }

  if (has_required_columns(cohort, Skilift:::required_columns$hetsnps)) {
    lift_hetsnps(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores
    )
  }

  # oncotable doesn't need every single column, just any one of them
  if (has_required_columns(cohort, Skilift:::required_columns$oncotable, any = TRUE)) {

    cohort <- create_oncotable(
      cohort = cohort,
      outdir = oncotable_dir,
      cores = cores
    )

    events_tbl = lift_filtered_events(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      return_table = TRUE
    )
  
  } else if (has_required_columns(cohort, Skilift:::required_columns$filtered_events)) {
    events_tbl = lift_filtered_events(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      return_table = TRUE
    )
  }

  if (has_required_columns(cohort, Skilift:::required_columns$multiplicity)) {
    lift_multiplicity(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      is_germline = FALSE,
      node_metadata = node_metadata,
      field = multiplicity_field
    )
  }

  if (has_required_columns(cohort, Skilift:::required_columns$germline_multiplicity)) {
    lift_multiplicity(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      is_germline = TRUE,
      node_metadata = node_metadata,
      field = multiplicity_field
    )
  }

  if (has_required_columns(cohort, Skilift:::required_columns$segment_width)) {
    lift_segment_width_distribution(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      annotations = annotations
    )
  }

  if (has_required_columns(cohort, Skilift:::required_columns$variant_qc)) {
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
  }

  if (has_required_columns(cohort, Skilift:::required_columns$pp_plot)) {
    lift_pp_plot(
      cohort = cohort, 
      output_data_dir = output_data_dir
    )
  }
  

  if (has_required_columns(cohort, Skilift:::required_columns$karyotype)) {
    warning("not implemented yet")
  }

  if (has_required_columns(cohort, Skilift:::required_columns$aggregated_events)) {
    warning("not implemented yet")
  }

  if (has_required_columns(cohort, Skilift:::required_columns$highlighted_events)) {
    warning("not implemented yet")
  }

  
}

lift_paired = function(cohort, output_data_dir, ...) {

  if (has_required_columns(cohort, Skilift:::required_columns$allelic_copy_number_graph)) {
    lift_copy_number_graph(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      settings = settings,
      is_allelic = TRUE,
      max_cn = max_cn,
      annotations = annotations
    )
  }
  
  if (has_required_columns(cohort, Skilift:::required_columns$total_copy_number_graph)) {
    lift_copy_number_graph(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      settings = settings,
      is_allelic = FALSE,
      max_cn = max_cn,
      annotations = annotations
    )
  }

  if (has_required_columns(cohort, Skilift:::required_columns$denoised_coverage)) {
    lift_denoised_coverage(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      coverage_field = coverage_field,
      color_field = color_field,
      bin.width = bin.width
    )
  }

  if (has_required_columns(cohort, Skilift:::required_columns$hetsnps)) {
    lift_hetsnps(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores
    )
  }


  if (has_required_columns(cohort, Skilift:::required_columns$filtered_events)) {
    lift_filtered_events(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores
    )
  # oncotable doesn't need every single column, just any one of them
  } else if (has_required_columns(cohort, Skilift:::required_columns$oncotable, any = TRUE)) {
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

  if (has_required_columns(cohort, Skilift:::required_columns$multiplicity)) {
    lift_multiplicity(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      is_germline = FALSE,
      node_metadata = node_metadata,
      field = multiplicity_field
    )
  }

  if (has_required_columns(cohort, Skilift:::required_columns$germline_multiplicity)) {
    lift_multiplicity(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      is_germline = TRUE,
      node_metadata = node_metadata,
      field = multiplicity_field
    )
  }

  if (has_required_columns(cohort, Skilift:::required_columns$pp_plot)) {
    lift_pp_plot(
      cohort = cohort, 
      output_data_dir = output_data_dir
    )
  }

  if (has_required_columns(cohort, Skilift:::required_columns$segment_width)) {
    lift_segment_width_distribution(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores,
      annotations = annotations
    )
  }

  if (has_required_columns(cohort, Skilift:::required_columns$signatures)) {
    lift_signatures(
      cohort = cohort,
      output_data_dir = output_data_dir,
      cores = cores
    )
  }

  if (has_required_columns(cohort, Skilift:::required_columns$variant_qc)) {
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
  }
}

