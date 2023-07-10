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
    to_list = function() {
      return(list(
        plot_id = as.character(self$plot_id),
        sample = ifelse(is.null(self$sample), NULL, as.character(self$sample)),
        type = as.character(self$type),
        source = as.character(self$source),
        title = as.character(self$title),
        visible = as.logical(self$visible),
        figure = ifelse(is.null(self$figure), NULL, as.character(self$figure)),
        server = ifelse(is.null(self$server), NULL, as.character(self$server)),
        uuid = ifelse(is.null(self$uuid), NULL, as.character(self$uuid))
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
      self$plots <- lapply(1:nrow(plots), function(i) {
        plot <- plots[i, ]
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
        return(args)
      })
    },
    to_list = function() {
      return(list(
        patient_id = as.character(self$patient_id),
        path = as.character(self$path),
        ref = as.character(self$ref),
        tags = as.character(self$tags),
        plots = lapply(self$plots, function(plot) plot$to_list()),
        path = as.character(self$path)
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
        print(data[[x]]$plots)
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
    to_json = function() {
      data <- vector("list", length(self$patients))
      names(data) <- sapply(self$patients, `[[`, "patient_id")

      for (i in seq_along(self$patients)) {
        patient <- self$patients[[i]]

        plots <- lapply(patient$plots, function(plot) {
          plot_data <- list(
            type = plot$type,
            source = plot$source,
            title = plot$title,
            visible = plot$visible
          )

          if (!is.null(plot$sample)) plot_data$sample <- plot$sample
          if (!is.null(plot$figure)) plot_data$figure <- plot$figure
          if (!is.null(plot$server) || !is.null(plot$uuid)) {
            plot_data$visible <- TRUE
          }
          if (!is.null(plot$server)) plot_data$server <- plot$server
          if (!is.null(plot$uuid)) plot_data$uuid <- plot$uuid

          return(plot_data)
        })
        
        data[[i]] <- list(
          description = patient$tags,
          reference = patient$ref,
          plots = plots
        )
      }

      return(toJSON(data, auto_unbox = TRUE, pretty = TRUE))
    }
  )
)

datafiles.json <- "~/projects/pgv/public/datafiles.json"
datadir <- "~/projects/pgv/public/data"
pgvdir <- "~/projects/pgv/public"

db <- PGVdb$new(datafiles.json, datadir)
db$drop_patient("TEST")
db$to_json()
filtered_patients <- db$filter_by_patient_id("E")
