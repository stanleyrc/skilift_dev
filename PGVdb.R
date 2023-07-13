library(R6)
library(jsonlite)
library(data.table)

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
    # Function to load JSON and assign to the metadata and plots attributes
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
    write_json = function() {
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
      for (i in 1:nrow(self$plots)) {
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
      print(json_string)
      write(json_string, file = json_file)
    }
  )
)

datafiles.json <- "~/projects/pgv/public/datafiles.json"
datadir <- "~/projects/pgv/public/data"
publicdir <- "~/projects/pgv/public"
settings <- "~/projects/pgv/public/settings.json"

db <- PGVdb$new(datafiles.json, datadir, publicdir, settings)
# db$drop_patient("TEST")
db$write_json()
# writeLines(db$to_json(), "test.datafiles.json")
# filtered_patients <- db$filter_by_patient_id("E")
