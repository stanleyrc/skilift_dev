library(R6)
library(jsonlite)
library(data.table)

#' @name cov2arrowPGV
#' @description
#'
#' Prepares an scatter plot arrow file with coverage info for PGV (https://github.com/mskilab/pgv)
#'
#' @param cov input coverage data (GRanges)
#' @param field which field of the input data to use for the Y axis
#' @param output_file output file path.
#' @param ref the name of the reference to use. If not provided, then the default reference that is defined in the meta.js file will be loaded.
#' @param cov.color.field a field in the input GRanges object to use to determine the color of each point
#' @param overwrite (logical) by default, if the output path already exists, it will not be overwritten.
#' @param meta.js path to JSON file with metadata for PGV (should be located in "public/settings.json" inside the repository)
#' @param bin.width (integer) bin width for rebinning the coverage (default: 1e4)
#' @author Alon Shaiber, Max Chao
#' @export
cov2arrowPGV <- function(cov,
                         field = "ratio",
                         output_file = "coverage.arrow",
                         ref = "hg19",
                         meta.js = NULL,
                         ...) {
  if (!file.exists(output_file)) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      stop('You must have the package "arrow" installed in order for converting a
           coverage file to arrow file to work. Please install it.')
    }
    message("Converting coverage format")
    dat <- cov2cov.js(cov,
      meta.js = meta.js,
      js.type = "PGV", field = field,
      ref = ref, ...
    )
    message("Done converting coverage format")
    if (!is.null(meta.js)) {
      ref_meta <- get_ref_metadata_from_PGV_json(meta.js, ref)
      setkey(ref_meta, "chromosome")
      # create a map
      # 3.981s
      map_cols <- data.table(
        color = unique(ref_meta$color),
        numcolor = color2numeric(unique(ref_meta$color))
      )
      dat$color <- merge(ref_meta[dat$seqnames], map_cols,
        by = "color", sort = FALSE
      )$numcolor
    } else {
      # no cov.color.field and no meta.js so set all colors to black
      dat$color <- 0
    }
    outdt <- dat[, .(x = new.start, y = get(field), color)]
    # if there are any NAs for colors then set those to black
    outdt[is.na(color), color := 0]
    # remove NAs
    outdt <- outdt[!is.na(y)]

    # sort according to x values (that is what PGV expects)
    outdt <- outdt[order(x)]

    message("Writing arrow file (using write_feather)")
    arrow_table <- arrow::Table$create(outdt,
      schema = arrow::schema(
        x = arrow::float32(),
        y = arrow::float32(),
        color = arrow::float32()
      )
    )
    arrow::write_feather(arrow_table, output_file)
  } else {
    message('arrow file, "', output_file, '" already exists.')
  }
  return(output_file)
}


PGVdb <- R6Class("PGVdb",
  private = list(
    json_file = NULL
  ),
  public = list(
    metadata = NULL,
    plots = NULL,
    datadir = NULL,
    publicdir = NULL,
    settings = NULL,
    initialize = function(json_file, datadir, publicdir, settings) {
      private$json_file <- json_file
      self$load_json(json_file)
      self$datadir <- datadir
      self$publicdir <- publicdir
      self$settings <- settings
    },

    load_json = function(json_file) {
      json_data <- jsonlite::fromJSON(json_file)

      # Extract metadata
      metadata_list <- lapply(names(json_data), function(patient_id) {
        patient_data <- json_data[[patient_id]]
        ref <- patient_data$reference
        tags <- list(patient_data$description) # Store tags as a list

        data.table(
          patient.id = patient_id,
          ref = ref,
          tags = tags
        )
      })

      self$metadata <- data.table::rbindlist(metadata_list)
      setnames(self$metadata, old = names(self$metadata), new = c("patient.id", "ref", "tags")) # Rename columns

      plots_list <- lapply(names(json_data), function(patient_id) {
        patient_data <- json_data[[patient_id]]
        patient_plots <- patient_data$plots

        plot_data <- data.table(
          sample = patient_plots$sample,
          type = patient_plots$type,
          source = patient_plots$source,
          visible = patient_plots$visible,
          title = patient_plots$title,
          figure = patient_plots$figure,
          server = patient_plots$server,
          uuid = patient_plots$uuid
        )

        data.table(
          patient.id = rep(patient_id, nrow(plot_data)),
          plot_data
        )
      })

      self$plots <- data.table::rbindlist(plots_list, fill = TRUE)
    },

    update_datafiles_json = function() {
      json_file <- private$json_file
      # Create a backup file with timestamp
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      backup_file <- paste0(json_file, ".", timestamp)
      file.copy(json_file, backup_file)

      # Build the JSON structure
      json_data <- list()

      # Add metadata section
      metadata <- list(
        "description" = self$metadata$tags[[1]],
        "reference" = self$metadata$ref
      )
      json_data$DEMO <- metadata

      # Add plots section
      plots <- list()
      for (i in seq_len(nrow(self$plots))) {
        plot_entry <- list(
          "type" = self$plots$type[i],
          "visible" = self$plots$visible[i],
          "source" = self$plots$source[i],
          "title" = self$plots$title[i]
        )
        if (!is.na(self$plots$sample[i])) {
          plot_entry$sample <- self$plots$sample[i]
        }
        if (!is.na(self$plots$figure[i])) {
          plot_entry$figure <- self$plots$figure[i]
        }
        if (!is.na(self$plots$server[i])) {
          plot_entry$server <- self$plots$server[i]
        }
        if (!is.na(self$plots$uuid[i])) {
          plot_entry$uuid <- self$plots$uuid[i]
        }
        plots[[i]] <- plot_entry
      }
      json_data$DEMO$plots <- plots

      # Convert JSON to string and write to file
      json_string <- jsonlite::toJSON(json_data, auto_unbox = TRUE, pretty = TRUE)
      write(json_string, file = json_file)
    },

    to_datatable = function(filter = NULL) {
      filtered_datatable <- data.table()

      if (!is.null(filter)) {
        column_name <- filter[[1]]
        filter_string <- filter[[2]]

        if (column_name %in% names(self$metadata)) {
          subset_metadata <- self$metadata[grepl(filter_string, self$metadata[[column_name]])]
          filtered_datatable <- merge(subset_metadata, self$plots, by = "patient.id")
        } else if (column_name %in% names(self$plots)) {
          subset_plots <- self$plots[grepl(filter_string, self$plots[[column_name]])]
          filtered_datatable <- merge(self$metadata, subset_plots, by = "patient.id")
        } else {
          stop(paste("Column", column_name, "not found in PGVdb."))
        }
      } else {
        filtered_datatable <- merge(self$metadata, self$plots, by = "patient.id")
      }

      return(filtered_datatable)
    },

    add_plots = function(plots_to_add, overwrite_files = FALSE) {
      new_plots <- data.table::setDT(plots_to_add) # Convert to data.table if required

      # Check if required columns exist
      required_columns <- c("patient.id", "sample", "type", "source", "path", "visible")
      if (!all(required_columns %in% names(new_plots))) {
        print("Required columns not found, creating an empty data.table with required columns instead...")
        return(data.table::data.table(patient.id = character(),
          sample = character(),
          type = character(),
          source = character(),
          path = character(),
          visible = logical(),
          stringsAsFactors = FALSE))
      }

      # Loop through each row of the plots_to_add table
      for (i in seq_len(nrow(new_plots))) {
        plot <- new_plots[i, ]

        # Check if patient exists, if not, create patient directory
        patient_dir <- file.path(self$datadir, plot$patient.id)
        if (!dir.exists(patient_dir)) {
          dir.create(patient_dir)
        }

        # Check if plot source does not exist for that patient or if overwrite = TRUE
        plot_file <- file.path(patient_dir, plot$source)
        if (file.exists(plot_file) && !overwrite_files) {
          message(paste0("Plot source file already exists for patient ", plot$patient.id, ". Specify overwrite_files = TRUE to overwrite."))
        } else {
          # Use the plot type to determine the conversion function to use
          if (plot$type == "genome") {
            self$create_ggenome_json(plot)
          } else if (plot$type == "scatterplot") {
            self$create_cov_arrow(plot)
          } else if (plot$type == "walk") {
            self$create_gwalk_json(plot)
          }

          # Save plot to patient directory
          if (file.exists(plot$path)) {
            file.copy(plot$path, plot_file, overwrite = TRUE)
          }
        }
      }

      # Update metadata and plots tables
      self$update_datafiles_json()
    },

    remove_plots = function(plots_to_remove, delete_files = FALSE) {
      remove_plots <- data.table::setDT(plots_to_remove) # Convert to data.table if required

      # Determine if only patient IDs are provided
      only_patient_ids <- "patient.id" %in% names(remove_plots) && !(c("sample", "source") %in% names(remove_plots))

      # Check if required columns exist
      required_columns <- if (only_patient_ids) "patient.id" else c("patient.id", "sample", "source")
      if (!all(required_columns %in% names(remove_plots))) {
        return("Error: Required columns (patient.id OR patient.id, source, sample) not found in plots_to_remove data.table.")
      }

      # Remove plots from PGVdb$plots
      if (only_patient_ids) {
        self$plots <- self$plots[!patient.id %in% remove_plots$patient.id]
      } else {
        self$plots <- self$plots[!(patient.id %in% remove_plots$patient.id & source %in% remove_plots$source)]
      }

      # Update datafiles.json
      self$update_datafiles_json()

      # Call validate to remove any patients that have no plots
      self$validate()

      # If delete_files is TRUE, remove plot files and patient directories if empty
      if (delete_files) {
        # Loop through each row of the plots_to_remove table
        for (i in seq_len(nrow(remove_plots))) {
          plot <- remove_plots[i, ]

          # Remove plot files
          if (only_patient_ids) {
            patient_plots <- self$plots[patient.id == plot$patient.id, .(patient.id, source)]
          } else {
            patient_plots <- self$plots[patient.id == plot$patient.id & source %in% plot$source, .(patient.id, source)]
          }

          # Loop through each plot file and remove it
          for (j in seq_len(nrow(patient_plots))) {
            plot_file <- file.path(self$datadir, patient_plots$patient.id[j], patient_plots$source[j])
            if (file.exists(plot_file)) {
              file.remove(plot_file)
            }
          }

          # Check if patient directories are empty and remove if so
          unique_patients <- unique(patient_plots$patient.id)
          for (patient in unique_patients) {
            patient_dir <- file.path(self$datadir, patient)
            if (dir.exists(patient_dir) && length(dir(paste0(patient_dir, "/*"))) == 0) {
              dir.remove(patient_dir)
            }
          }
        }
      }
    },

    validate = function() {
      # Check if all patients have at least one plot, otherwise remove patient from metadata
      patients_without_plots <- self$metadata[!patient.id %in% unique(self$plots$patient.id), patient.id]
      if (length(patients_without_plots) > 0) {
        self$metadata <- self$metadata[!patient.id %in% patients_without_plots]
      }

      # Check if all source files exist in proper directory
      missing_files <- self$plots[!file.exists(file.path(self$datadir, patient.id, source)), .(patient.id, source)]

      # Check if all plots have required column values != NA or NULL
      missing_values <- self$plots[sapply(.SD, function(x) any(is.na(x) | is.null(x))), .SDcols = c("patient.id", "source")]

      # Construct error message with missing files and missing values
      error_message <- ""
      if (nrow(missing_files) > 0) {
        error_message <- paste(error_message, "Missing Files:\n")
        error_message <- paste(error_message, paste(missing_files$patient.id, missing_files$source, sep = " - "), collapse = "\n")
      }
      if (nrow(missing_values) > 0) {
        error_message <- paste(error_message, "Missing Values:\n")
        error_message <- paste(error_message, paste(missing_values$patient.id, missing_values$source, sep = " - "), collapse = "\n")
      }

      # Return error message if there are any missing files or values
      if (error_message != "") {
        return(error_message)
      }
    },

    create_cov_arrow = function(plot_metadata, overwrite = FALSE) {
      cov_json_path <- file.path(self$datadir, plot_metadata$source)
      if (!file.exists(cov_json_path) || overwite == TRUE) {
        if (file.exists(plot_metadata$path)) {
          cov2arrowPGV(plot_metadata$path,
            field = plot_metadata$path,
            meta.js = self$settings,
            ref = plot_metadata$path,
            output_file = cov_json_path
          )
        } else {
          warning(paste0(
            "Input coverage file does not exist for name: ",
            plot_metadata$sample,
            " so no coverage will be generated."
          ))
        }
      } else {
        message(plot_metadata$path, " already exists! Set overwrite = TRUE if you want to overwrite it.")
      }
    },

    create_ggraph_json = function(plot_metadata, overwrite = FALSE) {
      ggraph_json_path <- file.path(self$datadir, plot_metadata$source)
      if (!file.exists(ggraph_json_path) || overwrite) {
        print(paste0("reading in ", plot_metadata$path))
        # TODO: at some point we need to do a sanity check to see that a valid rds of gGraph was provided
        if (grepl(plot_metadata$path, pattern = ".rds")) {
          ggraph <- readRDS(plot_metadata$path)
        } else {
          message("Expected .rds ending for gGraph. Attempting to read anyway: ", plot_metadata$path)
          ggraph <- readRDS(plot_metadata$path)
        }
        if (any(class(ggraph) == "gGraph")) {
          seq_lengths <- parse.js.seqlengths(
            self$settings,
            js.type = "PGV",
            ref = plot_metadata$ref
          )
          # check for overlap in sequence names
          ggraph.reduced <- gg[seqnames %in% names(seq_lengths)]
          if (length(ggraph.reduced) == 0) {
            stop(sprintf(
              'There is no overlap between the sequence names in the reference
              used by PGV and the sequences in your gGraph. Here is an
              example sequence from your gGraph: "%s". And here is an
              example sequence from the reference used by gGnome.js: "%s"',
              seqlevels(ggraph$nodes$gr)[1], names(seq_lengths)[1]
            ))
          }
          # sedge.id or other field
          if (plot_metadata$annotation) {
            # probably check for other cid.field names?
            # field = 'sedge.id'
            refresh(ggraph[seqnames %in% names(seq_lengths)])$json(
              filename = ggraph_json_path,
              verbose = TRUE,
              annotation = plot_metadata$annotation
            ) # ,
            # cid.field = field)
          } else {
            refresh(ggraph[seqnames %in% names(seq_lengths)])$json(
              filename = ggraph_json_path,
              verbose = TRUE
            )
          }
        } else {
          warning(plot_metadata$path, " rds read was not a gGraph")
        }
      } else {
        warning("file ", ggraph_json_path, "already exists. Set overwrite = TRUE if you want to overwrite it.")
      }
    },

    create_gwalk_json = function(overwrite = FALSE) {
      gwalk_json_path <- file.path(self$datadir, plot_metadata$source)
      if (!file.exists(gwalk_json_path) || overwrite == TRUE) {
        print(paste0("reading in ", plot_metadata$path))
        # TODO: at some point we need to do a sanity check to see that a valid rds of gWalk was provided
        gwalk <- readRDS(plot_metadata$path) %>%
          refresh()
        if (gwalk$length == 0) {
          warning(sprintf("Zero walks in gWalk .rds file provided for sample %s!
                          No walks json will be produced!", plot_metadata$sample))
          return(NA)
        }
        gwalk$json(
          filename = gwalk_json_path, verbose = TRUE,
          annotation = plot_metadat$annotation,
          include.graph = FALSE
        )
      } else {
        message(gwalk_json_path, "already exists! Set overwrite = TRUE if you want to overwite it.")
      }
    }
  )
)

datafiles.json <- "~/projects/pgv/public/datafiles.json"
datadir <- "~/projects/pgv/public/data"
publicdir <- "~/projects/pgv/public"
settings <- "~/projects/pgv/public/settings.json"

db <- PGVdb$new(datafiles.json, datadir, publicdir, settings)
# db$drop_patient("TEST")
db$update_datafiles_json()
# writeLines(db$to_json(), "test.datafiles.json")
# filtered_patients <- db$filter_by_patient_id("E")
