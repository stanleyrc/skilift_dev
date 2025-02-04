#' @description
#' Creates a multiplicity data table from SNV copy number data, handling both somatic and germline cases.
#'
#' @param snv_cn Path to SNV copy number RDS file or data.table object
#' @param is_germline Logical indicating if mutations are germline (default: FALSE)
#' @param field Column name to use for copy number values (default: "total_copies")
#'
#' @return data.table containing processed mutation data
#' @export
create_multiplicity <- function(snv_cn, is_germline = FALSE, field = "total_copies") {
    if (is.character(snv_cn)) {
        if (!grepl("\\.rds$", snv_cn)) {
            message("Expected .rds ending for mutations. Attempting to read anyway: ", snv_cn)
        }
        mutations.dt <- as.data.table(readRDS(snv_cn))
    } else {
        stop("Input must be a file path to an RDS file.")
    }

    if (!any(class(mutations.dt) == "data.table")) {
        stop("Input must be a data.table.")
    }

    # Process mutations
    setnames(mutations.dt, old = "VAF", new = "vaf", skip_absent = TRUE)
    mutations.dt <- mutations.dt[!is.na(get(field)), ]
    mutations.dt[start == end, end := end + 1]
    
    if ("strand" %in% colnames(mutations.dt)) {
        mutations.dt[, strand := NULL]
    }

    # Create annotation string
    mut_ann <- ""
    annotation_fields <- list(
        annotation = "Type",
        gene = "Gene",
        variant.c = "Variant",
        variant.p = "Protein_variant",
        variant.g = "Genomic_variant",
        vaf = "VAF",
        alt = "Alt_count",
        ref = "Ref_count",
        normal.alt = "Normal_alt_count",
        normal.ref = "Normal_ref_count",
        FILTER = "Filter"
    )

    for (col in names(annotation_fields)) {
        if (col %in% colnames(mutations.dt)) {
            mut_ann <- paste0(mut_ann, annotation_fields[[col]], ": ", mutations.dt[[col]], "; ")
        }
    }
    
    mutations.dt[, annotation := mut_ann]
    
    return(mutations.dt)
}

#' @title Convert Multiplicity Data to Intervals
#' @description
#' Converts multiplicity data into a list of intervals and settings for visualization
#'
#' @param multiplicity data.table containing mutation data with seqnames, start, and end columns
#' @param field Column name to use for y-axis values (default: "total_copies")
#' @param settings Path to settings JSON file
#' @param node_metadata Additional columns to include in node data (default: NULL)
#' @param reference_name Reference genome name (default: "hg19")
#'
#' @return List containing settings and intervals for visualization
#' @export
multiplicity_to_intervals <- function(
    multiplicity,
    field = "total_copies",
    settings = internal_settings_path,
    node_metadata = NULL,
    reference_name = "hg19"
) {
    # Load chromosome lengths from settings
    settings_data <- jsonlite::fromJSON(settings)
    chrom_lengths <- as.data.table(
        settings_data$coordinates$sets[[reference_name]])[
        , .(chromosome, startPoint, endPoint)
    ]
    setnames(chrom_lengths, c("seqnames", "start", "end"))

    # Ensure chromosome naming consistency
    if (nrow(chrom_lengths[grepl("chr", seqnames), ]) > 0) {
        chrom_lengths[!grepl("chr", seqnames), 
                     seqnames := paste0("chr", seqnames)]
    }

    # Set y-values from specified field
    multiplicity[, y_value := get(field)]

    # Convert to GRanges
    gr <- if (nrow(chrom_lengths[grepl("chr", seqnames), ]) > 0) {
        dt2gr(multiplicity[order(seqnames, start), ]) %>%
            sortSeqlevels() %>%
            gr.chr()
    } else {
        dt2gr(multiplicity[order(seqnames, start), ]) %>%
            sortSeqlevels() %>%
            gr.nochr()
    }

    # Validate ranges
    if (any(gr@seqinfo@seqlengths > 
            chrom_lengths[seqnames %in% names(seqlengths(gr))]$end)) {
        stop(paste("Ranges exceed chromosome lengths in", reference_name))
    }

    # Create graph object and convert to data.table
    jab <- gG(nodes = gr)
    node_cols <- c("snode.id", "y_value", node_metadata)
    node_dt <- gr2dt(jab$nodes$gr[, node_cols])

    # Create final node data
    intervals <- node_dt[, .(
        chromosome = seqnames,
        startPoint = start,
        endPoint = end,
        iid = snode.id,
        title = snode.id,
        type = "interval",
        y = y_value,
        annotation = node_dt$annotation
    )]

    # Construct return list
    list(
        settings = list(
            y_axis = list(
                title = "copy number",
                visible = TRUE
            )
        ),
        intervals = intervals,
        connections = data.table()
    )
}


#' @name lift_multiplicity
#' @title lift_multiplicity
#' @description
#' Create multiplicity JSON files for all samples in a cohort
#'
#' @param cohort Cohort object containing sample information
#' @param is_germline Logical indicating if mutations are germline (default: FALSE)
#' @param node_metadata Additional columns to include in node data (default: c("gene", "feature_type", "annotation", "REF", "ALT", "variant.c", "variant.p", "vaf", "transcript_type", "impact", "rank"))
#' @param output_data_dir Base directory for output files
#' @param cores Number of cores for parallel processing (default: 1)
#' @return None
#' @export
lift_multiplicity <- function(
    cohort,
    is_germline = FALSE,
    node_metadata = c("gene", "feature_type", "annotation", "REF", "ALT", "variant.c", "variant.p", "vaf", "transcript_type", "impact", "rank"),
    output_data_dir,
    field = "total_copies",
    cores = 1
) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }
    
    if (!dir.exists(output_data_dir)) {
        dir.create(output_data_dir, recursive = TRUE)
    }
    
    # Determine which column to use based on is_germline
    snv_cn_col <- if(is_germline) "germline_multiplicity" else "multiplicity"
    
    # Validate required column exists
    if (!snv_cn_col %in% names(cohort$inputs)) {
        stop(sprintf("Missing required column in cohort: %s", snv_cn_col))
    }
    
    # Get reference name from cohort
    reference_name <- cohort$reference_name
    if (is.null(reference_name)) {
        stop("Reference name not found in cohort object")
    }
    
    # Process each sample in parallel
    mclapply(seq_len(nrow(cohort$inputs)), function(i) {
        row <- cohort$inputs[i,]
        pair_dir <- file.path(output_data_dir, row$pair)
        
        if (!dir.exists(pair_dir)) {
            dir.create(pair_dir, recursive = TRUE)
        }
        
        # Determine output filename based on is_germline
        out_file <- file.path(
            pair_dir,
            if(is_germline) "germline_mutations.json" else "mutations.json"
        )
        
        tryCatch({
            # Create multiplicity data.table
            mult_dt <- create_multiplicity(
                snv_cn = row[[snv_cn_col]],
                is_germline = is_germline,
                field = field
            )
            
            # Convert to intervals
            intervals_list <- multiplicity_to_intervals(
                multiplicity = mult_dt,
                reference_name = reference_name,
                node_metadata = node_metadata,
                field = field
            )
            
            # Write to JSON
            jsonlite::write_json(
                intervals_list,
                out_file,
                pretty = TRUE,
                auto_unbox = TRUE,
                digits = 4
            )
            
        }, error = function(e) {
            warning(sprintf("Error processing %s: %s", row$pair, e$message))
        })
    }, mc.cores = cores, mc.preschedule = FALSE)
    
    invisible(NULL)
}
