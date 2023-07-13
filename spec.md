I need you to write some R code to create a class that can represent a JSON object as a flat datatable. I will give you an example JSON so you can understand the structure as well as columns I want the resulting datatable to have. The class will also have some methods for CRUD operations which we will implement afterwards. 

Here is an example JSON:

```
{
  "DEMO": {
    "description": ["tumor_type=ESAD", "dataset=ESAD-UK", "organ=Esophagus", "lineage=Esophageal columnar epithelium", "gender=female", "cluster=INVDUP", "PLAG1 mut", "KRAS amp", "MUTYH germline pathogenic", "CDKN2A homdel"],
    "reference": "hg19",
    "plots": [
      {
        "type": "phylogeny",
        "source": "phylogeny.newick",
        "title": "Phylogenetic Information for DEMO",
        "visible": true
      },
      {
        "type": "anatomy",
        "source": "anatomy.json",
        "title": "Anatomy for DEMO",
        "visible": true,
        "figure": "figure.svg"
      },
      {
        "type": "genome",
        "source": "genome1.json",
        "title": "Tumor sample 7K",
        "visible": true,
        "sample": "7K"
      },
      {
        "type": "walk",
        "source": "walks.json",
        "title": "Walks for sample 7K",
        "visible": true,
        "sample": "7K"
      },
      {
        "type": "genome",
        "source": "genome.json",
        "title": "Tumor sample 7E",
        "visible": true,
        "sample": "7E"
      },
      {
        "type": "barplot",
        "source": "rpkm.arrow",
        "title": "RPKM Probability Distribution 7E",
        "visible": false,
        "sample": "7E"
      },
      {
        "type": "scatterplot",
        "source": "coverage.arrow",
        "title": "Coverage Distribution 7E",
        "visible": false,
        "sample": "7E"
      },
      {
        "type": "genome",
        "source": "genome2.json",
        "title": "Tumor sample 7G",
        "visible": true,
        "sample": "7G"
      },
      {
        "type": "scatterplot",
        "source": "coverage2.arrow",
        "title": "Coverage Distribution 7G",
        "visible": false,
        "sample": "7G"
      },
      {
        "type": "bigwig",
        "visible": true,
        "sample": "7L",
        "server": "http://higlass.io",
        "uuid": "AOR9BgKaS4WPX7esuBM4sQ"
      }
    ]
  },
  "TEST": {
    "description": "patient.id=TEST",
    "reference": "hg19",
    "plots": [
      {
        "type": "genome",
        "source": "TEST.json",
        "title": "TEST",
        "visible": true,
        "sample": "TEST"
      }
    ]
  }
}
```

Here is the schema of the PGVdb class object:

- plots: data.table with columns, where optional indicates that NA or NULL values are allowed 
	- plot.id: Int, an incremental id
	- patient.id: String, maps to root keys of the JSON (e.g "DEMO", "TEST")
	- ref: String, maps to "reference"
    - sample: String
    - type: String
    - source: Filepath
    - visible: Boolean
    - title: String(optional)
    - figure: String(optional)
    - server: URL(optional)
    - uuid: String(optional)
    - tags: List[String](optional), maps to "description"
- datadir: Filepath
- publicdir: Filepath
- settings: Filepath


I want you fill in the write_json method for the following class. This method should write the tables contained in the metadata and plots attribute to a single json that has a structure I will specify shortly.

Here is the class definition:

```
library(R6)
library(jsonlite)
library(data.table)

PGVdb <- R6Class("PGVdb",
  public = list(
    metadata = NULL,
    plots = NULL,
    datadir = NULL,
    publicdir = NULL,
    settings = NULL,

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

      # Extract plots
      plots_list <- lapply(names(json_data), function(patient_id) {
        patient_data <- json_data[[patient_id]]
        patient_plots <- patient_data$plots

        plot_ids <- seq_along(patient_plots)

        sample <- sapply(patient_plots, function(plot) ifelse("sample" %in% names(plot), plot$sample, NA))
        type <- sapply(patient_plots, function(plot) ifelse("type" %in% names(plot), plot$type, NA))
        source <- sapply(patient_plots, function(plot) ifelse("source" %in% names(plot), plot$source, NA))
        visible <- sapply(patient_plots, function(plot) ifelse("visible" %in% names(plot), plot$visible, NA))
        title <- sapply(patient_plots, function(plot) ifelse("title" %in% names(plot), plot$title, NA))
        figure <- sapply(patient_plots, function(plot) ifelse("figure" %in% names(plot), plot$figure, NA))
        server <- sapply(patient_plots, function(plot) ifelse("server" %in% names(plot), plot$server, NA))
        uuid <- sapply(patient_plots, function(plot) ifelse("uuid" %in% names(plot), plot$uuid, NA))

        data.table(
          plot.id = plot_ids,
          patient.id = rep(patient_id, length(patient_plots)),
          sample = sample,
          type = type,
          source = source,
          visible = visible,
          title = title,
          figure = figure,
          server = server,
          uuid = uuid
        )
      })

      self$plots <- data.table::rbindlist(plots_list)
    },

    # Constructor
    initialize = function(json_file, datadir, publicdir, settings) {
      self$load_json(json_file)
      self$datadir <- datadir
      self$publicdir <- publicdir
      self$settings <- settings
    },
    write_json = function() {
      # Create a backup file with timestamp
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      backup_file <- paste0(json_file, ".", timestamp)
      file.copy(json_file, backup_file)

      # write metadata and plots to a single json
    }
  )
)
```

Here is the JSON structure mapping:

```
{
  [patient.id: String]: {
    "description": [tags: List[String],
    "reference": [ref: String],
    "plots": [plots: List[object]
  },
}
```
