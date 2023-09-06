source("src/higlass_and_ggraph.R")  # import auxiliary functions for import datatables and granges bigwigs

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
#' @author Alon Shaiber
cov2arrowPGV = function(cov,
        field = "ratio",
        output_file = 'coverage.arrow',
        ref = 'hg19',
        cov.color.field = NULL,
        overwrite = FALSE,
        meta.js = NULL,
        ...){

    outdir = dirname(output_file)
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    if (!file.exists(output_file) | overwrite){
        if (!requireNamespace("arrow", quietly = TRUE)) {
            stop('You must have the package "arrow" installed in order for this function to work. Please install it.')
        }

        message('Converting coverage format')
        dat = cov2cov.js(cov, meta.js = meta.js, js.type = 'PGV', field = field,
                         ref = ref, cov.color.field = cov.color.field, ...)
        message('Done converting coverage format')

        if (!is.null(cov.color.field)){
            dat[, color := color2numeric(get(cov.color.field))]
        } else {
            if (!is.null(meta.js)){
                ref_meta = get_ref_metadata_from_PGV_json(meta.js, ref)
                setkey(ref_meta, 'chromosome')
                dat$color = color2numeric(ref_meta[dat$seqnames]$color)
            } else {
                # no cov.color.field and no meta.js so set all colors to black
                dat$color = 0
            }
        }

        outdt = dat[, .(x = new.start, y = get(field), color)]

        # if there are any NAs for colors then set those to black
        outdt[is.na(color), color := 0]

        # remove NAs
        outdt = outdt[!is.na(y)]

        # sort according to x values (that is what PGV expects)
        outdt = outdt[order(x)]

        message('Writing arrow file (using write_feather)')
        arrow_table = arrow::Table$create(outdt, schema = arrow::schema(x = arrow::float32(), y = arrow::float32(), color = arrow::float32()))
        arrow::write_feather(arrow_table, output_file, compression="uncompressed")
    } else {
        message('arrow file, "', output_file, '" already exists.')
    }
    return(output_file)
}

#' @title PGVdb Object
#'
#' @description
#' Class representing a PGV database. Contains metadata, plots, and methods
#' for interacting with the database and converting to and from JSON.
#'
#' @section Methods:
#' \code{initialize()} Initialize the PGVdb object
#' \code{load_json()} Load data from JSON file into metadata and plots
#' \code{update_datafiles_json()} Update the JSON data files
#' \code{to_datatable()} Convert metadata and plots to a single data table
#' \code{add_plots()} Add new plots to the database
#' \code{remove_plots()} Remove plots from the database
#' \code{validate()} Validate metadata and plot data
#' \code{create_cov_arrow()} Create coverage arrow plot JSON
#' \code{create_ggraph_json()} Create gGraph JSON
#' \code{create_gwalk_json()} Create gWalk JSON
#' \code{upload_to_higlass()} Upload file to higlass server
#' \code{init_pgv()} Download and launch PGV instance
#'
#' @export
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom gGnome parse.js.seqlengths refresh
#' @import parallel
#' @import R6
#' @import httr
#' @import data.table
PGVdb <- R6Class( "PGVdb",
  private = list(
    #' @field datafiles_json_path (`character(1)`)
    #' Path to the datafiles.json
    datafiles_json_path = NULL
  ),
  public = list(
    #' @field metadata (`data.table`)\cr
    #' Data table containing patient metadata
    metadata = NULL,

    #' @field plots (`data.table`)\cr
    #' Data table containing plot metadata
    plots = NULL,

    #' @field datadir (`character(1)`).
    datadir = NULL,

    #' @field settings (`character(1)`).
    settings = NULL,

    #' @field higlass_metadata (`list`)\cr
    #' Attribute containing HiGlass metadata
    higlass_metadata = NULL,

    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param datafiles_json_path (`character(1)`)\cr
    #'   JSON file path.
    #' @param datadir (`character(1)`)\cr
    #'   Data directory path.
    #' @param settings (`character(1)`)\cr
    #'   Settings object path.
    #' @param higlass_metadata (`list`, optional)\cr
    #'   HiGlass metadata containing endpoint, username, and password.
    initialize = function(datafiles_json_path, datadir, settings, higlass_metadata = NULL) {
      private$datafiles_json_path <- datafiles_json_path
      self$load_json(datafiles_json_path)
      self$datadir <- datadir
      self$settings <- settings
      
      if (!is.null(higlass_metadata)) {
        self$higlass_metadata <- list(
          endpoint = ifelse(is.null(higlass_metadata$endpoint), "http://10.1.29.225:41800/", higlass_metadata$endpoint),
          username = ifelse(is.null(higlass_metadata$username), "admin", higlass_metadata$username),
          password = ifelse(is.null(higlass_metadata$password), "higlass_test", higlass_metadata$password)
        )
      } else {
        self$higlass_metadata <- list(
          endpoint = "http://10.1.29.225:41800/",
          username = "admin",
          password = "higlass_test"
        )
      }
    },

    #' @description
    #' Loads patient metadata and plot data from a JSON file.
    #'
    #' @param datafiles_json_path (`character(1)`)\cr
    #'   Datafiles.json file path.
    #' @return NULL.
    load_json = function(datafiles_json_path) {
      json_data <- jsonlite::fromJSON(datafiles_json_path)

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

      self$metadata <- rbindlist(metadata_list)
      setnames(self$metadata, old = names(self$metadata), new = c("patient.id", "ref", "tags")) # Rename columns

      plots_list <- lapply(names(json_data), function(patient_id) {
        patient_data <- json_data[[patient_id]]
        patient_plots <- patient_data$plots

        plot_data  <- copy(patient_plots)

        data.table(
          patient.id = rep(patient_id, nrow(plot_data)),
          plot_data
        )
      })

      self$plots <- rbindlist(plots_list, fill = TRUE)
      self$validate()
    },

    #' @description
    #' Update the JSON data files with current metadata and plot data.
    #' @return NULL.
    update_datafiles_json = function() {
      missing_data <- self$validate()
      if (!is.null(missing_data)) {
        stop(warning("Did not update datafiles because of malformed pgvdb object."))
      }

      datafiles_json_path <- private$datafiles_json_path
      # Create a backup file with timestamp
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      backup_file <- paste0(datafiles_json_path, ".", timestamp)
      file.copy(datafiles_json_path, backup_file)

      # Build the JSON structure
      json_data <- list()

      # Iterate over patient.ids
      for (patient_id in self$metadata$patient.id) {
        # Subset metadata for the current patient
        metadata <- self$metadata[self$metadata$patient.id == patient_id, ]
        description <- metadata$tags[[1]]
        reference <- metadata$ref

        # Subset plots for the current patient
        patient_plots <- self$plots[self$plots$patient.id == patient_id, ]
        plots <- list()
        for (i in seq_len(nrow(patient_plots))) {
          plot_entry <- list(
            "type" = patient_plots$type[i],
            "visible" = patient_plots$visible[i]
          )

          # Get the column names of patient_plots
          cols <- names(patient_plots)

          # Loop over the column names
          for (col in cols) {
            # Use [[ ]] to access the column by name and check if it's NA
            if (!is.na(patient_plots[[col]][i])) {
              # If it's not NA, add it to plot_entry
              plot_entry[[col]] <- patient_plots[[col]][i]
            }
          }

          plots[[i]] <- plot_entry
        }

        # Create patient data entry
        patient_data <- list(
          "description" = description,
          "reference" = reference,
          "plots" = plots
        )

        # Add the patient data to the main json_data list
        json_data[[patient_id]] <- patient_data
      }

      # Convert JSON to string and write to file
      json_string <- jsonlite::toJSON(json_data, auto_unbox = TRUE, pretty = TRUE)
      write(json_string, file = datafiles_json_path)
    },

    #' @description
    #' Convert metadata and plots to a single data table.
    #'
    #' @param filter (`list()`)\cr
    #'  Filter to apply. Consists of two elements: the column to filter on and the string to filter with.
    #'
    #' @return (`data.table`).
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

    #' @description
    #' Add new plots to the PGVdb.
    #'
    #' @param plots_to_add (`data.table`)\cr
    #'  Data table of plots to add. Must have a minimal set of columns: patient.id, x. 
    #'  x is either a list of the [server, uuid] or an object (GRanges, gWalk, gGraph, JSON), or a filepath
    #'  If adding a new patient, table must include a ref column.
    #' @param cores (`number(1)`)\cr
    #'  Number of cores to use for parallel execution
    #' @return NULL.
    add_plots = function(plots_to_add, cores=2) {
      new_plots <- data.table::setDT(plots_to_add) # Convert to data.table if required

      if (!("overwrite" %in% tolower(names(new_plots)))) {
        new_plots[, overwrite := FALSE]
      }
      # Check if required columns exist
      required_columns <- c("patient.id", "x")
      empty <- data.table::data.table(
        patient.id = character(),
        x = list(), # either a list of the [server, uuid] or an object (GRanges, gWalk, gGraph, JSON), or a filepath
        stringsAsFactors = FALSE
      )

      if (!all(required_columns %in% names(new_plots))) {
        warning("Required columns not found, creating an empty data.table (called 'empty') with required columns instead...")
        message("(All plots must have an 'x' column containing either list(list(server=, uuid=)) or an object (GRanges, gWalk, gGraph, JSON), or a filepath to a valid RDS data file)")

        return(empty)
      }

      if (any(!new_plots$patient.id %in% self$metadata$patient.id)) {
        if (!("ref" %in% names(new_plots)) | any(is.na(new_plots$ref)) | any(is.null(new_plots$ref))) {
          warning("You are trying to add a new patient without specifying its reference. Make sure all references are non-null")
          message("Creating an empty data.table (called 'empty') with required columns instead...")
          empty$ref <- character()
          return(empty)
        }
      }

      # Loop through each row of the plots_to_add table
      for (i in seq_len(nrow(new_plots))) {
        plot <- new_plots[i, ]

        # Check if patient dir exists, if not, create patient directory
        patient_dir <- file.path(self$datadir, plot$patient.id)
        if (!dir.exists(patient_dir)) {
          warning(paste("Patient directory does not exist. Creating a directory for new patient: ", plot$patient.id))
          dir.create(patient_dir)
        }

        if (is(plot$x[[1]], "GRanges")) {
          if (!"type" %in% names(plot) || is.na(plot$type) || !(plot$type %in% c("scatterplot", "bigwig"))) {
            warning("Plot type must be specified by the user for GRanges objects and should be either type='scatterplot' or type='bigwig'. This plot will be skipped.")
            next
          }
          if (plot$type == "scatterplot") {
            if (!"source" %in% names(plot) || is.na(plot$source)) {
              plot$source <- 'coverage.arrow'
            }
          } else if (plot$type == "bigwig") {
            unique_filename <- paste("file_", gsub(" ", "_", gsub(":|-", "_", Sys.time())), ".bw", sep = "")
            plot$source  <- unique_filename
          }
        } else if (is(plot$x[[1]], "gGraph")) {
          plot$type <- 'genome'
          if (!"source" %in% names(plot) || is.na(plot$source)) {
            plot$source <- 'genome.json'
          }
        } else if (is(plot$x[[1]], "gWalk")) {
          plot$type <- 'walk'
          if (!"source" %in% names(plot) || is.na(plot$source)) {
            plot$source <- 'walks.json'
          }
        } else if (is.list(plot$x) && length(plot$x[[1]]) == 2) {
          plot$type <- 'bigwig'
          plot$server <- plot$x[[1]][[1]]
          plot$uuid  <- plot$x[[1]][[2]]
          new_plots[i, "server"] <- plot$server
          new_plots[i, "uuid"] <- plot$uuid
        } else if (file.exists(plot$x) && tools::file_ext(plot$x) == 'rds') {
          rds_object <- readRDS(plot$x)
          if (is(rds_object, "GRanges")) {
            if (!"type" %in% names(plot) || is.na(plot$type) || !(plot$type %in% c("scatterplot", "bigwig"))) {
              warning("Plot type must be specified by the user for GRanges objects and should be either type='scatterplot' or type='bigwig'. This plot will be skipped.")
              next
            }
            if (plot$type == "scatterplot") {
              if (!"source" %in% names(plot) || is.na(plot$source)) {
                plot$source <- 'coverage.arrow'
              }
            } else if (plot$type == "bigwig") {
              unique_filename <- paste("file_", Sys.time(), ".bw", sep = "")
              plot$source  <- unique_filename
            }
          } else if (is(rds_object, "gGraph")) {
            plot$type <- 'genome'
            if (!"source" %in% names(plot) || is.na(plot$source)) {
              plot$source <- 'genome.json'
            }
          } else if (is(rds_object, "gWalk")) {
            plot$type <- 'walk'
            if (!"source" %in% names(plot) || is.na(plot$source)) {
              plot$source <- 'walks.json'
            }
          }
        } else if (file.exists(plot$x) && tools::file_ext(plot$x) == 'json') {
          if (!"type" %in% names(plot) || is.na(plot$type)) {
            warning("Type must be specified by the user for JSON files. This plot will be skipped.")
            next
          } else {
            if (!"source" %in% names(plot) || is.na(plot$source)) {
              plot$source  <- basename(plot$x)
            }
          }
        }

        # Check if source already exists, if so, increment
        if (!is.null(plot$source)) {
          source_full_path <- file.path(self$datadir, plot$patient.id, plot$source)
          if (plot$overwrite && file.exists(source_full_path)) {
            warning("Existing source for plot will be overwritten.")
          } else if (!plot$overwrite && file.exists(source_full_path)) {
            base_name <- tools::file_path_sans_ext(plot$source)
            ext <- tools::file_ext(plot$source)
            counter <- 2
            while (file.exists(file.path(self$datadir, plot$patient.id, paste0(base_name, counter, ".", ext)))) {
              counter <- counter + 1
            }
            plot$source <- paste0(base_name, counter, ".", ext)
          }
          new_plots[i, "source"] <- plot$source
        }
        new_plots[i, "type"] <- plot$type
      } # Break the loop into three pieces to parallize the plot creation

      # Define the function to create plot files
      create_plot_file <- function(plot) {
        tryCatch({
          plot_file <- file.path(self$datadir, plot$patient.id, plot$source)

          if (plot$type == "bigwig" && !is.null(plot$source)) {

            # field contains the column name in GRanges that corresponds to bigwig scores
            if (is.null(plot$field)) {
              warning("The column name in GRanges containing scores for bigwig was not specified, using 'foreground' as default")
              score_col_name <- "foreground" # Use default value
            } else {
              score_col_name <- plot$field
            }

            # Find ref chrom lengths with matching patient.id
            settings_data <- jsonlite::fromJSON(self$settings)
            chrom_lengths <- as.data.table(settings_data$coordinates$sets[[plot$ref]])[,.(chromosome,startPoint,endPoint)]
            colnames(chrom_lengths) = c("seqnames","start","end")
            chrom_lengths[!grepl("chr",seqnames), seqnames := paste0("chr",seqnames)] # weird fix because hg38_chr does not have chr on Y and M

            if (endsWith(plot$x, ".rds")) {
              bigwig_grange <- readRDS(plot$x)
            } else if (inherits(plot$x, "GRanges")) {
              bigwig_grange <- plot$x
            }

            gr2bw(gr = bigwig_grange,
                  output_filepath = plot$source,
                  score_col_name = score_col_name,
                  chrom_lengths = chrom_lengths)

            # Call upload_to_higlass function with parameters
            uuid <- self$upload_to_higlass(patient.id = plot$patient.id,
                                      datafile = plot$source,
                                      filetype = "bigwig",
                                      datatype = "vector",
                                      name = basename(plot$source),
                                      coordSystem = plot$ref)

            # Bigwig should be deleted after uploading
            if (plot$overwrite && !is.null(uuid)) {
              warning("Removing bigwig file after successful upload to higlass")
              file.remove(plot$source)
            }
          } else {
            if (grepl("\\.json$", plot$x)) {
              file.copy(plot$x, plot_file)
            } else {
              # Use the plot type to determine the conversion function to use
              if (plot$type == "genome") {
                self$create_ggraph_json(plot)
              } else if (plot$type == "scatterplot") {
                self$create_cov_arrow(plot)
              } else if (plot$type == "walk") {
                self$create_gwalk_json(plot)
              }
            }
          }
        }, error = function(e) {
          message("Error in creating plot file: ", e$message)
          print(e)
          traceback()
        })
      }

      if (!any(is.null(plot$source))) {
            # Use mclapply to create the plot files in parallel
            parallel::mclapply(seq_len(nrow(new_plots)), function(i) {
                                 plot <- new_plots[i, ]
                                 create_plot_file(plot)
        }, mc.cores = cores)
      }

      # Last piece of the loop
      for (i in seq_len(nrow(new_plots))) {
        plot <- new_plots[i, ]

        if (!is.null(plot$source)) {
          if (plot$type == "bigwig") {
            next
          }
          source_full_path  <- file.path(self$datadir, plot$patient.id, plot$source)
          if (!file.exists(source_full_path)) {
            warning("File does not exist: ", source_full_path, " skipping adding to pgvdb...")
            next
          }
        }

        common_columns <- intersect(names(plot), names(self$plots))
        extended_data <- plot[, ..common_columns]

        is_duplicate <- duplicated(rbind(self$plots, extended_data, fill = TRUE))

        if (any(is_duplicate)) {
          warning("Row already in plots, skipping (plot file will still be overwritten if overwrite flag was set): ")
          print(plot)
        } else {
          self$plots <- rbind(self$plots, extended_data, fill = TRUE)
        }

        if (plot$patient.id %in% self$metadata$patient.id) {
          # If patient.id exists, update tags and ref
          if ("tags" %in% colnames(plot)) {
            self$metadata[self$metadata$patient.id == plot$patient.id, "tags"] <- plot$tags
          }
          if ("ref" %in% colnames(plot) & !is.null(plot$ref) & !is.na(plot$ref) ) {
            self$metadata[self$metadata$patient.id == plot$patient.id, "ref"] <- plot$ref
          }
        } else {
          if ("tags" %in% colnames(plot)) {
            self$metadata <- rbind(self$metadata, data.frame(patient.id = plot$patient.id, tags = plot$tags, ref = plot$ref), fill = TRUE)
          } else {
            self$metadata <- rbind(self$metadata, data.frame(patient.id = plot$patient.id, ref = plot$ref), fill = TRUE)
          }
        }
      }

      self$update_datafiles_json()
    },

    #' @description
    #' Remove plots from the PGVdb.
    #'
    #' @param plots_to_remove (`data.table`)\cr 
    #'   Data table of plots to remove
    #' @param delete (`logical(1)`)\cr 
    #'   Delete plot files if TRUE.
    #'
    #' @return NULL
    remove_plots = function(plots_to_remove, delete = FALSE) {
      plots_to_remove <- data.table::setDT(plots_to_remove) # Convert to data.table if required

      # Determine if only patient IDs are provided
      is_remove_patients <- length(names(plots_to_remove)) == 1 && names(plots_to_remove) == "patient.id"
      is_remove_server <- "server" %in% names(plots_to_remove) && any(!is.na(plots_to_remove$server)) && any(!is.na(plots_to_remove$uuid))

      # Check if required columns exist
      required_columns <- if (is_remove_patients) "patient.id" else c("patient.id", "source")
      server_required_columns <- c("patient.id", "server", "uuid")
      if (!(all(required_columns %in% names(plots_to_remove)) || all(server_required_columns %in% names(plots_to_remove)))) {
        return("Error: Required columns (patient.id OR patient.id, source OR patient.id, server, uuid) not found in plots_to_remove data.table.")
      }

      # If delete is TRUE, remove plot files and patient directories if empty
      if (delete) {
        # Loop through each row of the plots_to_remove table
        print("Deleting plots from data directory...")
        for (i in seq_len(nrow(plots_to_remove[!is.na("source")]))) {
          plot <- plots_to_remove[i, ]

          # Select plots/patients to remove
          if (is_remove_patients) {
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
              message(patient_dir, "has no plots, deleting...")
              unlink(patient_dir)
              self$metadata <- self$metadata[patient.id != patient]
            }
          }
        }
      }

      # Remove plots from PGVdb$plots
      print("Removing plots from PGVdb...")
      initial_rows <- nrow(self$plots)

      if (is_remove_patients) {
        self$plots <- self$plots[!self$plots$patient.id %in% plots_to_remove$patient.id | is.na(self$plots$patient.id) | is.na(plots_to_remove$patient.id), ]
      } else {
        if (!is.null(plots_to_remove$source)) {
          self$plots <- self$plots[
                                   !(self$plots$patient.id %in% plots_to_remove$patient.id & !is.na(self$plots$patient.id) & !is.na(plots_to_remove$patient.id) &
                                     self$plots$source %in% plots_to_remove$source & !is.na(self$plots$source) & !is.na(plots_to_remove$source))
                                   ]
        }
        if (is_remove_server) {
          self$plots <- self$plots[
                                   !(self$plots$patient.id %in% plots_to_remove$patient.id & !is.na(self$plots$patient.id) & !is.na(plots_to_remove$patient.id) &
                                     self$plots$server %in% plots_to_remove$server & !is.na(self$plots$server) & !is.na(plots_to_remove$server) &
                                     self$plots$uuid %in% plots_to_remove$uuid & !is.na(self$plots$uuid) & !is.na(plots_to_remove$uuid))
                                   ]
        }
      }

      final_rows <- nrow(self$plots)
      if (initial_rows == final_rows) {
        warning("No rows were removed. Are you sure the row exists in plots?")
      }
      self$update_datafiles_json()

    },

    #' @description
    #' Validate metadata and plot data.
    #'
    #' @return NULL.
    validate = function() {
      # Check if there are any duplicate columns in self$plots
      if (any(duplicated(colnames(self$plots)))) {
        warning("Duplicate columns found in plots table. Removing duplicate...")
        print(colnames(self$plots))
        self$plots[, which(duplicated(names(self$plots))) := NULL]
      }
      
      # Check if there are any duplicate columns in self$metadata
      if (any(duplicated(colnames(self$metadata)))) {
        warning("Duplicate columns found in plots table. Removing duplicate...")
        print(colnames(self$metadata))
        self$metadata[, which(duplicated(names(self$metadata))) := NULL]
      }

      # Check if all patients have at least one plot, otherwise remove patient from metadata
      patients_without_plots <- self$metadata[!patient.id %in% unique(self$plots$patient.id), patient.id]
      if (length(patients_without_plots) > 0) {
        self$metadata <- self$metadata[!patient.id %in% patients_without_plots]
      }

      # Check if all source files exist in proper directory and server is not null
      missing_files <- self$plots[(is.na(server) | is.null(server)) & !file.exists(file.path(self$datadir, patient.id, source)), .(patient.id, source)]

      missing_servers <- self$plots[(is.na(source) | is.null(source)) & (is.na(server) | is.null(server) & (is.na(uuid) | is.null(uuid)))]
      # All plots must have a patient.id and either a source or server and uuid
      missing_values <- self$plots[
        (is.na(patient.id) | is.null(patient.id)) |
        (
          (is.na(source) | is.null(source)) &
          (
            (is.na(server) | is.null(server)) &
            (is.na(uuid) | is.null(uuid))
          )
        ),
        .SDcols = c("patient.id", "source", "server", "uuid"),
      ]

      # Construct error message and table with missing source and missing values
      error_message <- ""
      missing_data <- data.table()
      if (nrow(missing_files) > 0) {
        error_message <- paste(error_message, "Missing Files:\n")
        error_message <- paste(error_message, paste(missing_files$patient.id, missing_files$source, sep = " - "), collapse = "\n")
        missing_data <- rbind(missing_data, missing_files, fill=TRUE)
      }
      if (nrow(missing_servers) > 0) {
        error_message <- paste(error_message, "Missing Servers (or uuids for the servers):\n")
        error_message <- paste(error_message, paste(missing_servers$patient.id, missing_servers$server, missing_servers$uuid, sep = " - "), collapse = "\n")
        missing_data <- rbind(missing_data, missing_servers, fill=TRUE)
      }
      if (nrow(missing_values) > 0) {
        error_message <- paste(error_message, "Missing Values:\n")
        error_message <- paste(error_message, paste(missing_values$patient.id, missing_values$source, sep = " - "), collapse = "\n")
        missing_data <- rbind(missing_data, missing_values, fill=TRUE)
      }

      # Coerce visible to be boolean
      self$plots$visible = as.logical(self$plots$visible)

      # Return error message if there are any missing files or values
      if (error_message != "") {
        warning(error_message)
        print("Returning data.table with the invalid rows...")
        return(missing_data)
      }
    },

    #' @description
    #' Create coverage arrow plot JSON file.
    #'
    #' @param plot_metadata (`data.table`)\cr 
    #'   Plot metadata.
    #'
    #' @return NULL.
    create_cov_arrow = function(plot_metadata) {
      cov_json_path <- file.path(
        self$datadir,
        plot_metadata$patient.id,
        plot_metadata$source
      )

      if (!("field" %in% names(plot_metadata))) {
          stop(warning("Please include a 'field' column which indicates the column name that contains the coverage data."))
      }

      if (!file.exists(cov_json_path) || plot_metadata$overwrite) {
        if (is(plot_metadata$x[[1]], "GRanges")) {
          cov2arrowPGV(plot_metadata$x[[1]],
            field = plot_metadata$field,
            meta.js = self$settings,
            ref = plot_metadata$ref,
            output_file = cov_json_path
          )
        } else if (file.exists(plot_metadata$x)) {
          cov2arrowPGV(plot_metadata$x,
            field = plot_metadata$field,
            meta.js = self$settings,
            ref = plot_metadata$ref,
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
        message(cov_json_path, " already exists! Set overwrite = TRUE if you want to overwrite it.")
      }
    },


    #' @description
    #' Create gGraph JSON file.
    #'
    #' @param plot_metadata (`data.table`)\cr 
    #'   Plot metadata.
    #'
    #' @return NULL.
    create_ggraph_json = function(plot_metadata) {
      ggraph_json_path <- file.path(
        self$datadir,
        plot_metadata$patient.id,
        plot_metadata$source
      )
      if (!file.exists(ggraph_json_path) || plot_metadata$overwrite) {
        if (is(plot_metadata$x[[1]], "gGraph")) {
          ggraph <- plot_metadata$x[[1]]
        } else {
          message(paste0("reading in ", plot_metadata$x))
          if (grepl(plot_metadata$x, pattern = ".rds")) {
            ggraph <- readRDS(plot_metadata$x)
          } else {
            message("Expected .rds ending for gGraph. Attempting to read anyway: ", plot_metadata$x)
            ggraph <- readRDS(plot_metadata$x)
          }
        }
        if (any(class(ggraph) == "gGraph")) {
          seq_lengths <- gGnome::parse.js.seqlengths(
            self$settings,
            js.type = "PGV",
            ref = plot_metadata$ref
          )
          # check for overlap in sequence names
          ggraph.reduced <- ggraph[seqnames %in% names(seq_lengths)]
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
          if ("annotation" %in% colnames(plot_metadata)) {
            # probably check for other cid.field names?
            # field = 'sedge.id'
            gGnome::refresh(ggraph[seqnames %in% names(seq_lengths)])$json(
              filename = ggraph_json_path,
              verbose = TRUE,
              annotation = plot_metadata$annotation
            # cid.field = field
            )
          } else {
            gGnome::refresh(ggraph[seqnames %in% names(seq_lengths)])$json(
              filename = ggraph_json_path,
              verbose = TRUE
            )
          }
        } else {
          warning(plot_metadata$x, " rds read was not a gGraph")
        }
      } else {
        warning("file ", ggraph_json_path, "already exists. Set overwrite = TRUE if you want to overwrite it.")
      }
    },

    #' @description
    #' Create gWalk JSON file
    #'
    #' @param plot_metadata (`data.table`)\cr 
    #'   Plot metadata.
    #'
    #' @return NULL.
    create_gwalk_json = function(plot_metadata) {
      gwalk_json_path <- file.path(self$datadir, plot_metadata$patient.id, plot_metadata$source)
      if (!file.exists(gwalk_json_path) || plot_metadata$overwrite == TRUE) {
        if (is(plot_metadata$x[[1]], "gWalk")) {
          gwalk <- plot_metadata$x[[1]] %>% gGnome::refresh()
        } else {
          message(paste0("reading in ", plot_metadata$x))
          if (grepl(plot_metadata$x, pattern = ".rds")) {
            gwalk <- readRDS(plot_metadata$x) %>% gGnome::refresh()
          } else {
            message("Expected .rds ending for gWalk Attempting to read anyway: ", plot_metadata$x)
            gwalk <- readRDS(plot_metadata$x) %>% gGnome::refresh()
          }
        }

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
        message(gwalk_json_path, "already exists! Set overwrite = TRUE if you want to overwrite it.")
      }
    },

    #' @description
    #' Upload file to higlass server
    #'
    #' @return uuid
    #'  The uuid returned by the response from higlass
    upload_to_higlass = function(endpoint = self$higlass_metadata$endpoint, 
                                 patient.id = "TEST_HIGLASS",
                                 datafile, 
                                 filetype, 
                                 datatype, 
                                 coordSystem, 
                                 name, 
                                 username = self$higlass_metadata$username, 
                                 password = self$higlass_metadata$password) {
      # Define the API endpoint
      print(paste("Uploading datafile to higlass:", datafile))

      # Convert datafile argument to a multipart object
      datafile <- httr::upload_file(datafile)

      # Package other data as key-value pairs in a list
      body <- list(
                   'datafile' = datafile,
                   'filetype' = filetype,
                   'datatype' = datatype,
                   'coordSystem' = coordSystem,
                   'name' = name
      )

      url <- paste0(endpoint, "api/v1/tilesets/")

      # Create the response object
      response <- httr::POST(
                             url,
                             authenticate(username, password, "basic"),
                             body = body,
                             encode = "multipart"
      )

      # Parse the response
      response_content <- httr::content(response, "parsed")

      # Store the UUID
      uuid <- response_content$uuid
      filetype  <- response_content$filetype

      print(paste("UUID:", uuid))
      print(paste("filetype:", filetype))
      if (filetype == "bigwig") {
        new_higlass <-  data.table(
                                   patient.id = patient.id, 
                                   ref=coordSystem, 
                                   x = list(list(server = endpoint, uuid = uuid)), 
                                   visible=TRUE
        )
        self$add_plots(new_higlass)
      }

      return(uuid)
    },

    #' @description
    #' Remove file in higlass server
    #'
    #' @return httr:response
    delete_from_higlass = function(endpoint = "http://10.1.29.225:41800/", 
                                   patient.id = "TEST_HIGLASS",
                                   uuid,
                                   username = "admin", 
                                   password = "higlass_test") {
      # Define the API endpoint
      url <- paste0(endpoint, "api/v1/tilesets/", uuid, "/")

      # Create the response object
      response <- DELETE(
                         url,
                         authenticate(username, password, "basic")
      )
      response_content <- httr::content(response, "parsed")
      print(response_content)

      remove_higlass <-  data.table(patient.id = patient.id, server = endpoint, uuid = uuid)
      self$remove_plots(remove_higlass)
    },


    #' @description
    #' Download and instantiate a PGV instance with symlinked data
    #'
    #' @param pgv_dir (`character(1)`)\cr 
    #'   Directory where the pgv instance will be installed
    #' @param build (`logical(1)`)\cr 
    #'   Flag to indicate whether to build pgv or launch local instance
    #'
    #' @return NULL
    init_pgv = function(pgv_dir, build = FALSE) {
      init_script_path  <- system.file("src", "init_pgv.sh", package="PGVdb")
      # Check if node is installed
      cmd <- paste("which", "node")
      is_installed <- length(system(cmd, intern = TRUE)) != 0
      if (!is_installed) {
        stop(paste(program_name, "is not installed. Please install it before proceeding."))
      } else {
        system(paste("bash", init_script_path, private$datafiles_json_path, self$datadir, self$settings, pgv_dir, build))
      }
    }
  )
)

