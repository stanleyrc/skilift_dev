# Import necessary libraries
library(R6)
library(jsonlite)

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
    tag = NULL,
    defaultChartType = NULL
    initialize = function(patient_id, type, source, visible, title = NULL, sample = NULL, figure = NULL, server = NULL, uuid = NULL, tag = NULL, defaultChartType = NULL) {
      self$plot_id <- as.character(length(self) + 1)
      self$patient_id <- patient_id
      self$sample <- sample
      self$type <- type
      self$source <- source
      self$title <- title
      self$visible <- visible
      self$figure <- figure
      self$server <- server
      self$uuid <- uuid
      self$tag <- tag
      self$defaultChartType <- defaultChartType
    },
    create_cov_arrow = function(datadir, cov_path, field, settings, ref, overwrite = FALSE) {
      cov_json_path <- file.path(datadir, self$source)
      if (!file.exists(cov_json_path) || overwite == TRUE) {
        if (file.exists(cov_path)) {
          # load gGraph
          cov2arrowPGV(cov_path,
            field = field,
            meta.js = settings,
            ref = ref,
            output_file = cov_json_path
          )
        } else {
          warning(paste0(
            "Input coverage file does not exist for name: ",
            self$sample,
            " so no coverage will be generated."
          ))
        }
      } else {
        message(cov_path, " already exists! Set overwrite = TRUE if you want to overwrite it.")
      }
    },
    create_ggraph_json = function(datadir, ggraph_path, settings, ref, overwrite = FALSE, annotation = NULL) {
      ggraph_json_path <- file.path(datadir, self$source)
      if (!file.exists(ggraph_json_path) || overwrite) {
        print(paste0("reading in ", ggraph_path))
        # TODO: at some point we need to do a sanity check to see that a valid rds of gGraph was provided
        if (grepl(ggraph_path, pattern = ".rds")) {
          ggraph <- readRDS(ggraph_path)
        } else {
          message("Expected .rds ending for gGraph. Attempting to read anyway: ", ggraph_path)
          ggraph <- readRDS(ggraph_path)
        }
        if (any(class(ggraph) == "gGraph")) {
          seq_lengths <- parse.js.seqlengths(
            settings,
            js.type = "PGV",
            ref = ref
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
          if (annotation) {
            # probably check for other cid.field names?
            # field = 'sedge.id'
            refresh(ggraph[seqnames %in% names(seq_lengths)])$json(
              filename = ggraph_json_path,
              verbose = TRUE,
              annotation = annotation
            ) # ,
            # cid.field = field)
          } else {
            refresh(ggraph[seqnames %in% names(seq_lengths)])$json(
              filename = ggraph_json_path,
              verbose = TRUE
            )
          }
        } else {
          warning(ggraph_path, " rds read was not a gGraph")
        }
      } else {
        warning("file ", ggraph_json_path, "already exists. Set overwrite = TRUE if you want to overwrite it.")
      }
    },
    create_gwalk_json = function(datadir, gwalk_path, annotation, overwrite) {
      gwalk_json_path <- file.path( datadir, self$source)
      if (!file.exists(gwalk_json_path) || overwrite == TRUE) {
        print(paste0("reading in ", gwalk_path))
        # TODO: at some point we need to do a sanity check to see that a valid rds of gWalk was provided
        gwalk <- readRDS(gwalk_path) %>%
          refresh()
        if (gwalk$length == 0) {
          warning(sprintf("Zero walks in gWalk .rds file provided for sample %s!
                          No walks json will be produced!", self$sample))
          return(NA)
        }
        gwalk$json(
          filename = gwalk_json_path, verbose = TRUE,
          annotation = annotation,
          include.graph = FALSE
        )
      } else {
        message(gwalk_json_path, "already exists! Set overwrite = TRUE if you want to overwite it.")
      }
    },
    delete_plot = function(datadir) {
      plot_json_path <- file.path(datadir, paste0(self$sample, self$source))
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
      self$plots <- lapply(seq_len(nrow(plots)), function(i) {
        plot <- plots[i, ]
        args <- list(
          patient_id = patient_id,
          type = plot[["type"]],
          source = plot[["source"]],
          visible = plot[["visible"]],
          sample = plot[["sample"]]
        )

        if ("title" %in% names(plot)) args$title <- plot[["title"]]
        if ("figure" %in% names(plot)) args$figure <- plot[["figure"]]
        if ("server" %in% names(plot)) args$server <- plot[["server"]]
        if ("uuid" %in% names(plot)) args$uuid <- plot[["uuid"]]
        if ("tag" %in% names(plot)) args$tag <- plot[["tag"]]
        if ("defaultChartType" %in% names(plot)) args$defaultChartType <- plot[["defaultChartType"]]

        args <- args[!sapply(args, is.null)]
        return(Plot$new(args))
      })
    },
    add_plots = function(plots_table, datadir, settings, overwrite=FALSE) {
      for (i in seq_len(nrow(plots_table))) {
        plot <- plots_table[i, ]
        args <- list(
          patient_id = self$patient_id,
          type = plot[["type"]],
          source = plot[["source"]],
          visible = plot[["visible"]],
          sample = plot[["sample"]]
        )

        if ("title" %in% names(plot) && !is.na(plot[["title"]])) args$title <- plot[["title"]]
        if ("figure" %in% names(plot) && !is.na(plot[["figure"]])) args$figure <- plot[["figure"]]
        if ("server" %in% names(plot) && !is.na(plot[["server"]])) args$server <- plot[["server"]]
        if ("uuid" %in% names(plot) && !is.na(plot[["uuid"]])) args$uuid <- plot[["uuid"]]
        if ("tag" %in% names(plot) && !is.na(plot[["tag"]])) args$tag <- plot[["tag"]]
        if ("defaultChartType" %in% names(plot) && !is.na(plot[["defaultChartType"]])) args$defaultChartType <- plot[["defaultChartType"]]

        args <- args[!sapply(args, is.null)]

        # Check if a Plot object with the same source already exists
        existing_plot <- self$plots[sapply(self$plots, function(p) identical(p$source, args$source))]

        if (is.null(existing_plot)) {
          new_plot <- Plot$new(args)

          if (new_plot$type == "scatterplot") {
            new_plot$create_cov_arrow(
              datadir = datadir
              cov_path = plot$path
              field = plot$field,
              settings = settings,
              ref = self$ref
              overwrite = overwrite,
            )
          } else if (new_plot$type == "genome") {
            new_plot$create_ggenome_json(
              datadir = datadir
              ggraph_path = plot$path
              settings = settings,
              ref = self$ref,
              overwrite = overwrite
              annotation = plot$annotation
            )
          } else if (new_plot$type == "walk") {
            new_plot$create_gwalk_json(
              datadir = datadir
              gwalk_path = plot$path
              annotation = plot$annotation
            )
          }

          self$plots <- c(self$plots, new_plot)
        } else {
          print("Plot with the same source already exists. Skipping addition.")
        }
      }
    },
    drop_plots = function(plot_ids, delete = FALSE) {
      for (plot_id in plot_ids) {
        # Retrieve the plot by its id
        matching_plots <- self$plots[which(sapply(self$plots, function(p) p$plot_id == plot_id))]

        if (delete) {
          for (plot in matching_plots) {
            plot$delete_plot()
          }
          self$plots <- self$plots[self$plots != plot]
        }
      }
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
    update = function() {
      filepath <- private$json_file
      backup_filepath <- paste0(
        filepath, ".",
        format(Sys.time(), "%Y%m%d_%H%M%S")
      )

      data <- vector("list", length(self$patients))
      names(data) <- sapply(self$patients, `[[`, "patient_id")

      for (i in seq_along(self$patients)) {
        patient <- self$patients[[i]]

        plots <- lapply(patient$plots, function(plot) {
          plot_data <- list(
            type = plot$type,
            visible = plot$visible
          )

          if (!is.na(plot$source)) plot_data$source <- plot$source
          if (!is.na(plot$title)) plot_data$title <- plot$title
          if (!is.na(plot$sample)) plot_data$sample <- plot$sample
          if (!is.na(plot$figure)) plot_data$figure <- plot$figure
          if (!is.na(plot$server)) plot_data$server <- plot$server
          if (!is.na(plot$uuid)) plot_data$uuid <- plot$uuid

          # Remove attributes with null values
          plot_data <- plot_data[!sapply(plot_data, is.null)]
          return(plot_data)
        })

        data[[i]] <- list(
          description = patient$tags,
          reference = patient$ref,
          plots = plots
        )
      }

      json_string <- toJSON(data, auto_unbox = TRUE, pretty = TRUE)

      file.copy(filepath, backup_filepath, overwrite = TRUE)
      write(json_string, file = filepath)
    }
  )
)

datafiles.json <- "~/projects/pgv/public/datafiles.json"
datadir <- "~/projects/pgv/public/data"
pgvdir <- "~/projects/pgv/public"

db <- PGVdb$new(datafiles.json, datadir)
db$drop_patient("TEST")
db$to_json()
writeLines(db$to_json(), "test.datafiles.json")
filtered_patients <- db$filter_by_patient_id("E")
