internal_settings_path = system.file("extdata", "test_data", "settings.json", package = "Skilift")

#' Converts a color from hexadecimal to numeric format
#'
#' @param color A character vector of colors in hexadecimal format (e.g., "#FF0000").
#' @return A numeric vector representing the colors.
color2numeric <- function(color) {
    # Remove the '#' character if present
    color <- gsub("#", "", color)
    # Convert hexadecimal to numeric
    return(as.numeric(strtoi(color, base = 16L)))
}


#'
#' @param ref_seqinfo_json Path to the JSON file containing metadata.
#' @param ref Optional reference name. If not provided, the default reference is used.
#' @return A data.table containing the sequence information for the specified reference.
get_ref_metadata <- function(ref_seqinfo_json, ref = NULL) {
    meta <- jsonlite::read_json(ref_seqinfo_json)
    
    if (!("coordinates" %in% names(meta))) {
        stop("Input meta file is not a proper settings.json format.")
    }
    
    coord <- meta$coordinates
    
    if (is.null(ref)) {
        if (!("default" %in% names(coord))) {
            stop(ref_seqinfo_json, " is missing a default coordinates attribute")
        }
        ref <- coord$default
        message('No reference name provided so using default: "', ref, '".')
    }
    
    if (!("sets" %in% names(coord))) {
        stop("Could not find reference sets in ", ref_seqinfo_json, "in attribute 'sets'.")
    }
    
    sets <- coord$sets
    
    if (!(ref %in% names(sets))) {
        stop(ref, " does not appear in set of references in ", ref_seqinfo_json)
    }
    
    seq.info <- data.table::rbindlist(sets[[ref]])
    return(seq.info)
}

#' @name granges_to_arrow_scatterplot
#' @description
#'
#' Converts GRanges into an arrow file
#'
#' @param gr_path input GRanges
#' @param field name of GRanges column to use as arrow Y axis
#' @param ref the name of the reference to use.
#' 
#' @importFrom arrow Table write_feather schema float32
#' @importFrom gGnome cov2cov.js
#' @param cov.color.field name of GRanges column containing color of each point
#' @param overwrite (logical) by default, if the output path already exists, it will not be overwritten.
#' @param ref_seqinfo_json path to JSON file with ref seqlengths
#' @param bin.width (integer) bin width for rebinning the coverage (default: 1e4)
#' @return arrow table
#' @author Alon Shaiber, Shihab Dider
granges_to_arrow_scatterplot = function(
    gr_path,
    field = "foreground.X",
    ref = 'hg19',
    cov.color.field = NULL,
    ref_seqinfo_json = internal_settings_path,
    bin.width = 1e4,
    ...){

    if (!requireNamespace("arrow", quietly = TRUE)) {
        stop('You must have the package "arrow" installed in order for this function to work. Please install it.')
    }

    message('Preparing GRanges for conversion to arrow format')
    dat = gGnome::cov2cov.js(
        gr_path,
        meta.js = ref_seqinfo_json,
        js.type = 'PGV',
        field = field,
        ref = ref,
        cov.color.field = cov.color.field,
        ...
    )

    if (!is.null(cov.color.field)){
        dat[, color := color2numeric(get(cov.color.field))]
    } else {
        if (!is.null(ref_seqinfo_json)){
            ref_meta = get_ref_metadata(ref_seqinfo_json, ref)
            data.table::setkey(ref_meta, 'chromosome')
            dat$color = color2numeric(ref_meta[dat$seqnames]$color)
        } else {
            dat$color = 0
        }
    }

    outdt = dat[, .(x = new.start, y = get(field), color)]

    # if there are any NAs for colors then set those to black
    outdt[is.na(color), color := 0]

    # remove NAs
    outdt = outdt[!is.na(y)]

    # sort according to x values (that is what front-end expects)
    outdt = outdt[order(x)]

    arrow_table = arrow::Table$create(
        outdt,
        schema = arrow::schema(
            x = arrow::float32(),
            y = arrow::float32(),
            color = arrow::float32()
        )
    )
    return(arrow_table)
}

#' @description
#' Create arrow scatterplot 
#'
#' @param plot_metadata data.table with Plot metadata for a single plot, columns = (patient.id, source, x (contains path to data file or reference to data file object itself), ref, overwrite).
#' @param datadir Path to data directory with patient directories.
#' @param settings Path to settings.json file.
#'
#' @return NULL.
create_scatterplot_arrow = function(plot_metadata, datadir, settings = internal_settings_path, color_field = NULL) {
    cov_json_path <- file.path(
        datadir,
        plot_metadata$patient.id,
        plot_metadata$source
    )

    if (!("field" %in% names(plot_metadata))) {
        stop(warning("Please include a 'field' column which indicates the column name that contains the coverage data."))
    }

    # check if the input is a nested list
    if (is(plot_metadata$x, "list")) {
        x <- plot_metadata$x[[1]]
    } else {
        x <- plot_metadata$x
    }

    temp_file <- NULL
    if (is(x, "GRanges")) {
        # save granges to a temp file and assign path to plot_metadata$x
        message("Saving GRanges to temp file")
        temp_file <- tempfile(fileext = ".rds")
        saveRDS(x, temp_file)
        x <- temp_file
    } else if (is(x, "character")) {
        message("Using path to GRanges")
    } else {
        stop(warning("Please provide a GRanges object or a path to a GRanges object."))
    }

    if (!file.exists(cov_json_path) || plot_metadata$overwrite) {
        if (file.exists(x)) {
            arrow_table <- granges_to_arrow_scatterplot(x,
                field = plot_metadata$field,
                ref_seqinfo_json = settings,
                ref = plot_metadata$ref,
                cov.color.field = color_field
            )
            message('Writing arrow file...')

            if (!dir.exists(dirname(cov_json_path))) {
                dir.create(dirname(cov_json_path), recursive = TRUE)
            }
            arrow::write_feather(arrow_table, cov_json_path, compression = "uncompressed")
        } else {
            warning(
                paste0(
                    "Input coverage file does not exist for name: ",
                    plot_metadata$sample,
                    " so no coverage will be generated."
                )
            )
        }
    } else {
        message(cov_json_path, " already exists! Set overwrite = TRUE if you want to overwrite it.")
    }

    if (!is.null(temp_file) && file.exists(temp_file)) {
        file.remove(temp_file)
    }
}

#' @name lift_denoised_coverage
#' @title lift_denoised_coverage
#' @description
#' Create denoised coverage arrow files for all samples in a cohort
#'
#' @param cohort Cohort object containing sample information
#' @param output_data_dir Base directory for output files
#' @param cores Number of cores for parallel processing (default: 1)
#' @return None
#' @export
lift_denoised_coverage <- function(cohort, output_data_dir, cores = 1) {
    if (!inherits(cohort, "Cohort")) {
        stop("Input must be a Cohort object")
    }
    
    if (!dir.exists(output_data_dir)) {
        dir.create(output_data_dir, recursive = TRUE)
    }
    
    # Validate required columns exist
    required_cols <- c("pair", "tumor_coverage")
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
        
        out_file <- file.path(pair_dir, "coverage.arrow")
        
        tryCatch({
            if (!is.null(row$tumor_coverage) && file.exists(row$tumor_coverage)) {
                # Create arrow table
                arrow_table <- granges_to_arrow_scatterplot(
                    gr_path = row$tumor_coverage,
                    field = "foreground.X",
                    ref = cohort$reference_name,
                )
                
                # Write arrow table
                message(sprintf("Writing coverage arrow file for %s", row$pair))
                arrow::write_feather(arrow_table, out_file, compression = "uncompressed")
            } else {
                warning(sprintf("Tumor coverage file missing for %s", row$pair))
            }
        }, error = function(e) {
            warning(sprintf("Error processing %s: %s", row$pair, e$message))
        })
    }, mc.cores = cores)
    
    invisible(NULL)
}
