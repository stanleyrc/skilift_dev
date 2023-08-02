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
#' \code{init_pgv()} Download and launch PGV instance
#'
#' @export
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom gGnome parse.js.seqlengths refresh
#' @import R6
#' @import data.table
PGVdb <- R6Class("PGVdb",
  private = list(
    #' @field json_file (`character(1)`)
    #' Path to the datafiles.json
    json_file = NULL
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

    #' @field publicdir (`character(1)`).
    publicdir = NULL,

    #' @field settings (`charater(1)`).
    settings = NULL,

    #' @description
    #' Creates a new instance of this [R6][R6::R6Class] class.
    #'
    #' @param json_file (`character(1)`)\cr
    #'   JSON file path.
    #' @param datadir (`character(1)`)\cr
    #'   Data directory path.
    #' @param publicdir (`character(1)`)\cr
    #'   Public directory path.
    #' @param settings (`character(1)`)\cr
    #'   Settings object path.
    initialize = function(json_file, datadir, publicdir, settings) {
      private$json_file <- json_file
      self$load_json(json_file)
      self$datadir <- datadir
      self$publicdir <- publicdir
      self$settings <- settings
    },

    #' @description
    #' Loads patient metadata and plot data from a JSON file.
    #'
    #' @param json_file (`character(1)`)\cr
    #'   Datafiles.json file path.
    #' @return NULL.
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

      self$metadata <- rbindlist(metadata_list)
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

      self$plots <- rbindlist(plots_list, fill = TRUE)
    },

    #' @description
    #' Update the JSON data files with current metadata and plot data.
    #' @return NULL.
    update_datafiles_json = function() {
      json_file <- private$json_file
      # Create a backup file with timestamp
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      backup_file <- paste0(json_file, ".", timestamp)
      file.copy(json_file, backup_file)

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

          if (!is.na(patient_plots$title[i])) {
            plot_entry$title <- patient_plots$title[i]
          }
          if (!is.na(patient_plots$source[i])) {
            plot_entry$source <- patient_plots$source[i]
          }
          if (!is.na(patient_plots$sample[i])) {
            plot_entry$sample <- patient_plots$sample[i]
          }
          if (!is.na(patient_plots$figure[i])) {
            plot_entry$figure <- patient_plots$figure[i]
          }
          if (!is.na(patient_plots$server[i])) {
            plot_entry$server <- patient_plots$server[i]
          }
          if (!is.na(patient_plots$uuid[i])) {
            plot_entry$uuid <- patient_plots$uuid[i]
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
      write(json_string, file = json_file)
    },

    #' @description
    #' Convert metadata and plots to a single data table.
    #'
    #' @param filter (`character(1)`)\cr
    #'  Filter to apply.
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
    #'   Data table of plots to add.
    #' @param overwrite (`logical(1)`)\cr
    #'   Overwrite existing files if TRUE.
    #' @return NULL.
    add_plots = function(plots_to_add, overwrite = FALSE) {
      new_plots <- data.table::setDT(plots_to_add) # Convert to data.table if required

      # Check if required columns exist
      required_columns <- c("patient.id", "type", "visible")
      empty <- data.table::data.table(
        patient.id = character(),
        type = character(),
        source = character(), # filename in pgv, e.g genome.json
        server = character(),
        uuid = character(),
        path = character(), # path to datafile, e.g HCC1954.gg.rds
        visible = logical(),
        stringsAsFactors = FALSE
      )
      if (!all(required_columns %in% names(new_plots)) ||
      !(("source" %in% names(new_plots) && "path" %in% names(new_plots)) || ("server" %in% names(new_plots) && "uuid" %in% names(new_plots)))) {
        warning("Required columns not found, creating an empty data.table with required columns instead...")
        message("(All plots must have either a non-null source field or non-null server and uuid fields)")

        return(empty)
      }

      if (any(!new_plots$patient.id %in% self$metadata$patient.id)) {
        if (!("ref" %in% names(new_plots)) | any(is.na(new_plots$ref)) | any(is.null(new_plots$ref))) {
          warning("You are trying to add a new patient without specifying its reference. Make sure all references are non-null")
          message("Creating an empty data.table with required columns instead...")
          empty$ref <- character()
          return(empty)
        }
      }

      # Loop through each row of the plots_to_add table
      for (i in seq_len(nrow(new_plots))) {
        plot <- new_plots[i, ]

        # Check if patient exists, if not, make sure ref exists, then create patient directory
        patient_dir <- file.path(self$datadir, plot$patient.id)
        if (!dir.exists(patient_dir)) {
          warning(paste("Patient directory does not exist. Creating a directory for new patient: ", plot$patient.id))
          dir.create(patient_dir)
        }

        # Check if plot source does not exist for that patient or if overwrite = TRUE
        plot_file <- file.path(patient_dir, plot$source)
        if (file.exists(plot_file) && is.na(plot$server) && !overwrite) {
          warning(paste0("Plot source file already exists for patient ", plot$patient.id, ". Set overwrite = TRUE to overwrite."))
        } else {
          # Use the plot type to determine the conversion function to use
          if (plot$type == "genome") {
            self$create_ggraph_json(plot, overwrite)
          } else if (plot$type == "scatterplot") {
            self$create_cov_arrow(plot, overwrite)
          } else if (plot$type == "walk") {
            self$create_gwalk_json(plot, overwrite)
          }
        }

        common_columns <- intersect(names(plot), names(self$plots))
        extended_data <- plot[, ..common_columns]

        is_duplicate <- duplicated(rbind(self$plots, extended_data, fill=TRUE))

        if (any(is_duplicate)) {
          warning("Skipping duplicate row: ")
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

      # Validate and update metadata and plots tables
      self$validate()
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
      remove_plots <- data.table::setDT(plots_to_remove) # Convert to data.table if required

      # Determine if only patient IDs are provided
      is_remove_patients <- length(names(remove_plots)) == 1 && names(remove_plots) == "patient.id"
      is_remove_server <- "server" %in% names(remove_plots)

      # Check if required columns exist
      required_columns <- if (is_remove_patients) "patient.id" else c("patient.id", "source")
      server_required_columns <- c("patient.id", "server", "uuid")
      if (!(all(required_columns %in% names(remove_plots)) || all(server_required_columns %in% names(remove_plots)))) {
        return("Error: Required columns (patient.id OR patient.id, source OR patient.id, server, uuid) not found in plots_to_remove data.table.")
      }

      # Call validate to remove any patients that have no plots
      self$validate()

      # If delete is TRUE, remove plot files and patient directories if empty
      if (delete) {
        # Loop through each row of the plots_to_remove table
        print("Deleting plots from data directory...")
        for (i in seq_len(nrow(remove_plots[!is.na("source")]))) {
          plot <- remove_plots[i, ]

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
      if (is_remove_patients) {
        self$plots <- self$plots[!patient.id %in% remove_plots$patient.id]
      } else {
        self$plots <- self$plots[!(patient.id %in% remove_plots$patient.id & source %in% remove_plots$source)]
        if (is_remove_server) {
          self$plots <- self$plots[!(patient.id %in% remove_plots$patient.id & server %in% remove_plots$server)]
        }
      }
    },

    #' @description
    #' Validate metadata and plot data.
    #'
    #' @return NULL.
    validate = function() {
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
        missing_data <- rbind(missing_data, missing_files)
      }
      if (nrow(missing_servers) > 0) {
        error_message <- paste(error_message, "Missing Servers (or uuids for the servers):\n")
        error_message <- paste(error_message, paste(missing_files$patient.id, missing_files$server, missing_files$uuid, sep = " - "), collapse = "\n")
        missing_data <- rbind(missing_data, missing_servers)
      }
      if (nrow(missing_values) > 0) {
        error_message <- paste(error_message, "Missing Values:\n")
        error_message <- paste(error_message, paste(missing_values$patient.id, missing_values$source, sep = " - "), collapse = "\n")
        missing_data <- rbind(missing_data, missing_values)
      }

      # Return error message if there are any missing files or values
      if (error_message != "") {
        warning(error_message)
        print("Returning data.table with invalid rows...")
        return(missing_data)
      } else {
        self$update_datafiles_json()
      }
    },

    #' @description
    #' Create coverage arrow plot JSON file.
    #'
    #' @param plot_metadata (`data.table`)\cr 
    #'   Plot metadata.
    #' @param overwrite (`logical(1)`)\cr 
    #'   Overwrite if file exists.
    #'
    #' @return NULL.
    create_cov_arrow = function(plot_metadata, overwrite = FALSE) {
      cov_json_path <- file.path(
        self$datadir,
        plot_metadata$patient.id,
        plot_metadata$source
      )

      if (!("field" %in% names(plot_metadata))) {
          stop(warning("Please include a 'field' column in the plot_metadata which indicates the column name that contains the coverage data."))
      }

      if (!file.exists(cov_json_path) || overwrite) {
        if (file.exists(plot_metadata$path)) {
          cov2arrowPGV(plot_metadata$path,
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
        message(plot_metadata$path, " already exists! Set overwrite = TRUE if you want to overwrite it.")
      }
    },


    #' @description
    #' Create gGraph JSON file.
    #'
    #' @param plot_metadata (`data.table`)\cr 
    #'   Plot metadata.
    #' @param overwrite (`logical(1)`)\cr 
    #'   Overwrite if file exists.
    #'
    #' @return NULL.
    create_ggraph_json = function(plot_metadata, overwrite = FALSE) {
      ggraph_json_path <- file.path(
        self$datadir,
        plot_metadata$patient.id,
        plot_metadata$source
      )
      if (!file.exists(ggraph_json_path) || overwrite) {
        message(paste0("reading in ", plot_metadata$path))
        # TODO: at some point we need to do a sanity check to see that a valid rds of gGraph was provided
        if (grepl(plot_metadata$path, pattern = ".rds")) {
          ggraph <- readRDS(plot_metadata$path)
        } else {
          message("Expected .rds ending for gGraph. Attempting to read anyway: ", plot_metadata$path)
          ggraph <- readRDS(plot_metadata$path)
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
            ) # ,
            # cid.field = field)
          } else {
            gGnome::refresh(ggraph[seqnames %in% names(seq_lengths)])$json(
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

    #' @description
    #' Create gWalk JSON file
    #'
    #' @param plot_metadata (`data.table`)\cr 
    #'   Plot metadata.
    #' @param overwrite (`logical(1)`)\cr 
    #'   Overwrite if file exists.
    #'
    #' @return NULL.
    create_gwalk_json = function(plot_metadata, overwrite = FALSE) {
      gwalk_json_path <- file.path(self$datadir, plot_metadata$patient.id, plot_metadata$source)
      if (!file.exists(gwalk_json_path) || overwrite == TRUE) {
        message(paste0("reading in ", plot_metadata$path))
        # TODO: at some point we need to do a sanity check to see that a valid rds of gWalk was provided
        gwalk <- readRDS(plot_metadata$path) %>% gGnome::refresh()
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
    #' Download and instantiate a PGV instance with symlinked data
    #'
    #' @param pgv_dir (`character(1)`)\cr 
    #'   Directory where the pgv instance will be installed
    #'
    #' @return NULL
    init_pgv = function(pgv_dir) {
      init_script_path  <- system.file("src", "init_pgv.sh", package="PGVdb")
      # Check if node is installed
      cmd <- paste("which", "node")
      is_installed <- length(system(cmd, intern = TRUE)) != 0
      if (!is_installed) {
        stop(paste(program_name, "is not installed. Please install it before proceeding."))
      } else {
        system(paste("bash", init_script_path, private$json_file, self$datadir, self$publicdir, self$settings, pgv_dir))
      }
    }
  )
)

