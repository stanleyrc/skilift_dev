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

#' @name gr2arrow
#' @description
#'
#' Converts GRanges into an arrow file
#'
#' @param gr input GRanges
#' @param field name of GRanges column to use as arrow Y axis
#' @param output_file output file path.
#' @param ref the name of the reference to use.
#' 
#' @importFrom arrow Table write_feather schema float32
#' @importFrom gGnome cov2cov.js
#' @param cov.color.field name of GRanges column containing color of each point
#' @param overwrite (logical) by default, if the output path already exists, it will not be overwritten.
#' @param ref_seqinfo_json path to JSON file with ref seqlengths
#' @param bin.width (integer) bin width for rebinning the coverage (default: 1e4)
#' @author Alon Shaiber, Shihab Dider
gr2arrow = function(
    gr_path,
    field = "ratio",
    output_file = 'coverage.arrow',
    ref = 'hg19',
    cov.color.field = NULL,
    overwrite = TRUE,
    ref_seqinfo_json = internal_settings_path,
    bin.width = 1e4,
    ...){

    if (!requireNamespace("arrow", quietly = TRUE)) {
        stop('You must have the package "arrow" installed in order for this function to work. Please install it.')
    }

    outdir = dirname(output_file)
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    if (file.exists(output_file) && !overwrite){
        message('arrow file, "', output_file, '" already exists. Set overwrite = TRUE to overwrite it.')
    } else {
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

        message('Writing arrow file...')
        arrow_table = arrow::Table$create(
            outdt,
            schema = arrow::schema(
                x = arrow::float32(),
                y = arrow::float32(),
                color = arrow::float32()
            )
        )
        arrow::write_feather(arrow_table, output_file, compression = "uncompressed")
    }
    return(output_file)
}

