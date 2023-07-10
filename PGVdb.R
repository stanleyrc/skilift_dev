# Import necessary libraries
library(R6)
library(jsonlite)

# Define Plot R6 class
Plot <- R6Class(
  "Plot",
  public = list(
    plot_id = NULL,
    sample = NULL,
    type = NULL,
    source = NULL,
    title = NULL,
    visible = NULL,
    figure = NULL,
    server = NULL,
    uuid = NULL,
    initialize = function(type, source, title, visible, sample = NULL, figure = NULL, server = NULL, uuid = NULL) {
      self$plot_id <- as.character(length(self) + 1)
      self$sample <- sample
      self$type <- type
      self$source <- source
      self$title <- title
      self$visible <- visible
      self$figure <- figure
      self$server <- server
      self$uuid <- uuid
    },
    toList = function() {
      return(list(
        plot_id = self$plot_id,
        sample = self$sample,
        type = self$type,
        source = self$source,
        title = self$title,
        visible = self$visible,
        figure = self$figure,
        server = self$server,
        uuid = self$uuid
      ))
    }
  )
)

# Define Patient R6 class
Patient <- R6Class(
  "Patient",
  public = list(
    patient_id = NULL,
    path = NULL,
    ref = NULL,
    tags = NULL,
    plots = list(),
    initialize = function(patient_id, ref, tags, plots, path) {
      self$patient_id <- patient_id
      self$path <- path
      self$ref <- ref
      self$tags <- tags
      self$plots <- apply(plots, 1, function(plot) {
        args <- list(
          type = plot[["type"]],
          source = plot[["source"]],
          title = plot[["title"]],
          visible = plot[["visible"]],
          sample = plot[["sample"]]
        )

        if ("figure" %in% names(plot)) args$figure <- plot[["figure"]]
        if ("server" %in% names(plot)) args$server <- plot[["server"]]
        if ("uuid" %in% names(plot)) args$uuid <- plot[["uuid"]]

        args <- args[!sapply(args, is.null)]
        return(do.call(Plot$new, args))
      })
    },
    toList = function() {
      return(list(
        patient_id = self$patient_id,
        path = self$path,
        ref = self$ref,
        tags = self$tags,
        plots = lapply(self$plots, function(plot) plot$toList()),
        path = self$path
      ))
    }
  )
)

# Define PatientsData R6 class
PGVdb <- R6Class(
  "PGVdb",
  private = list(
    json_file = NULL
  ),
  public = list(
    patients = list(),
    initialize = function(json_file, data_dir) {
      data <- fromJSON(json_file, flatten = TRUE)
      self$patients <- lapply(names(data), function(x) {
        tags <- if (is.null(data[[x]]$description)) NULL else data[[x]]$description
        Patient$new(
          patient_id = x,
          ref = data[[x]]$reference,
          tags = tags,
          plots = data[[x]]$plots,
          path = data_dir
        )
      })
      private$json_file <- json_file # Store the file path
    },
    filter_by_patient_id = function(id_term) {
      return(self$patients[sapply(self$patients, function(patient) {
        grepl(id_term, patient$patient_id)
      })])
    },
    add_patient = function(new_patient) {
      self$patients <- c(self$patients, list(new_patient))
    },
    drop_patient = function(patient_id) {
      self$patients <- self$patients[sapply(self$patients, function(patient) {
        patient$patient_id != patient_id
      })]
    },
    update_datafile = function() {
      # Validate the contents of self$patients
      are_all_patients <- sapply(self$patients, function(patient) {
        inherits(patient, "Patient")
      })

      if (!all(are_all_patients)) {
        stop("The patients list in the PGVdb object contains non-Patient objects.")
      }

      data <- lapply(self$patients, function(patient) {
        patient_data <- patient$toList()
        patient_data <- patient_data[!sapply(patient_data, is.null)]
        return(patient_data)
      })

      # Formulate timestamp and backup filename
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      backup_filename <- paste0(private$json_file, ".", timestamp)

      # Rename the current json file to the backup file
      file.rename(private$json_file, backup_filename)

      # Write the new, updated data to the original filename.
      write(toJSON(data, pretty = TRUE), private$json_file)
    }
  )
)

datafiles.json <- "~/projects/pgv/public/datafiles.json"
datadir <- "~/projects/pgv/public/data"
pgvdir <- "~/projects/pgv/public"

db <- PGVdb$new(datafiles.json, datadir)
db$drop_patient("TEST")
db$update_datafile()
filtered_patients <- db$filter_by_patient_id("E")
