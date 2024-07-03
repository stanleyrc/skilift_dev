# Skilift
Faciliates easier loading and manipulation of genomic data for use with [PGV](https://github.com/mskilab-org/pgv) and[Case Report](https://github.com/mskilab-org/case-report).

# Installation
- To install directly from github: `devtools::install_github("mskilab-org/skilift")`
- To install from source:
    1. `git clone https://github.com/mskilab-org/skilift.git`
    2. `devtools::load_all("path/to/clone")`

# Usage
The Skilift class manages a database of patient metadata and genomic plot data. It provides an interface to load, update, query, and manipulate the database programmatically.

The key methods allow converting to/from JSON for data storage, validating the data, adding/removing plots, and generating plot JSON for visualization. The metadata and plots are stored as data tables for easy manipulation.

The conversion methods can also be used in isolation without having to instantiate a Skilift object. 
This is useful for generating JSON files for use with Case Reports.

## Methods
The following methods are for the Skilift class:
- `initialize(datafiles_json_path, datadir, settings)`: Initialize a new Skilift object

- `update_datafiles_json()`: Update JSON data files on disk 

- `to_datatable(filter)`: Convert to a single data table, apply optional filter

- `add_plots(new_plots_dt)`: Add new plots to PGV

- `remove_plots(plots_to_remove_dt)`: Remove plots

- `list_higlass_tilesets(endpoint, username, password)`: List bigwig tilesets on higlass

- `validate()`: Validate metadata and plots

Methods for converting to JSON:
plot_metadata is a data table containing plot metadata and should have the following columns:
    - `patient.id`: Patient identifier (e.g '0124')
    - `source`: output datafile json filename (e.g complex.json)
    - `x`: filepath to the RDS object with the raw data (e.g ~/gg.rds)
    - `overwrite`: Whether to overwrite existing plot files
    - `ref`: Reference genome (e.g 'hg38')

datadir is the parent directory of the sample directory where the raw data file is stored 
e.g if the data file is stored in `data/0124/gg.rds` then the datadir is `data/`

settings is the path to the settings file (you can use the pgv one). It's required for parsing seqlengths

- create_ggraph_json(plot_metadata, datadir, settings)

- create_allelic_json(plot_metadata, datadir, settings)

- create_gwalk_json(plot_metadata, datadir, settings)

- create_somatic_json(plot_metadata, datadir, settings)

- create_ppfit_genome_json(plot_metadata, datadir, settings)

- create_distributions(case_reports_data_folder, common_folder, filter_patients = NULL, write_jsons = TRUE)

- create_ppfit_json(jabba_gg, path_obj, out_file = NULL, write_json = TRUE, overwrite = FALSE, return_table = FALSE, cores = 1)

## Fields

- `metadata`: Data table containing patient metadata 

- `plots`: Data table containing plot metadata

- `datadir`: Path to data directory

- `settings`: Path to settings file

- `higlass_metadata`: List containing endpoint, username, and password to Higlass server

### Constructor
```r
pgvdb <- Skilift$new(datafiles_json_path, datadir, settings)
```

Create a new Skilift object by passing the path to the datafiles.json, data directory, and settings file.

### update_datafiles_json
```r
pgvdb$update_datafiles_json()  
```

Update the datafiles JSON on disk with the current metadata and plots.

### to_datatable
```r
dt <- pgvdb$to_datatable(filter = c("patient.id", "HCC"))
```

Convert metadata and plots to a single data table, applying an optional filter.

### add_plots

```r
pgvdb$add_plots(plots_to_add, cores=2)  
```

Add new plots by passing a `data.table` containing required columns:

- `patient.id`: Patient identifier
- `x`: can be
    - list(list(server="", uuid="")): a list of the [server, uuid] 
    - list(GRanges), list(gWalk), list(gGraph): a list containing an object
    - Filepath to an RDS object
- `visible`: Whether plot is visible

If required columns are missing, an empty table will be returned with the
required column names. The `type` column (indiciating the type of plot:
scatterplot, genome, walk, bigwig, etc) is derived from `x` if not supplied by
the user. For `GRanges`, the `type` and `field` (the name of the score column)
must be specified by the user (either bigwig or scatterplot). The `source`
column (indicating the name of the plot file inside the pgv data directory) is
derived from the `type`, if not supplied by the user. 

Unless an `overwrite` column is set, it will not overwrite existing plot files,
instead it will just increment the filename of the new plot file (i.e
coverage.json -> coverage2.json). For GRanges containing bigwig data, the
GRange will be converted into a bigwig (with a unique name) and automaticaly
uploaded to the default mskilab higlass server (see above for changing to a
different server). If overwrite is set to TRUE, the converted bigwig file will
be deleted after uploading.

The `cores` parameter determines how many cores to use for parallel execution.
By default, it will not run in parallel (i.e use 1 core).

### remove_plots

```r
pgvdb$remove_plots(plots_to_remove, delete = FALSE)
```

Remove plots by passing a `data.table` with either:

- `patient.id`: Remove all plots for a patient 
- `patient.id` and `source`: Remove specific plots

Or alternatively with `patient.id`, `server`, and `uuid`.

The `delete` flag determines whether to also delete the source data files.

Passing just `patient.id` will remove all plots for that patient.

### validate

```r
pgvdb$validate()
```

Validate metadata and plot data, removing invalid entries. Is automatically
called when adding/removing plots. Useful if you make manual changes to the
pgvdb plots or metadata.


### init_pgv

```r 
pgvdb$init_pgv(pgv_dir, build=FALSE)
```

Initialize a pgv instance loaded with the data in pgvdb at `pgv_dir`. This
method will clone/pull the ![pgv repo](https://github.com/mskilab-org/pgv) and
create symlinks in the pgv data directory that point to your pgvdb data. 

If the `build` flag is set to `TRUE` it will build pgv instead of launching a
local instance (useful when running on a remote server or hpc). In that case,
you should set `pgv_dir` to be a directory inside whichever directory is served
by your remote server (e.g `public_html`).

### list_higlass_tilesets

```r 
pgvdb$list_higlass_tilesets(endpoint, username, password)
```

Return all bigwig tileset info on the higlass server as a data.table.

```r 
endpoint <- "http://10.1.29.225:8000/api/v1/tilesets/" # dev endpoint
tilesets <- pgvdb$list_higlass_tilesets(
    endpoint,
    username = "username_here",
    password = "password_here"
)
```



### upload_to_higlass

```r 
pgvdb$upload_to_higlass(endpoint, datafile, filetype, datatype, coordSystem, name, username, password)
```

Upload a file to the higlass server. Will also add the file to the current
pgvdb instance. Note that you will need to upload a `chromSizes.tsv` file
first, before uploading other files. 

```r 
endpoint <- "http://10.1.29.225:8000/api/v1/tilesets/" # dev endpoint
pgvdb$upload_to_higlass(
    endpoint,
    datafile = system.file("extdata", "test_data", "chromSizes.tsv", package = "Skilift"),
    filetype = "chromsizes-tsv",
    datatype = "chromsizes",
    coordSystem = "hg38",
    name = "hg38",
    username = "username_here",
    password = "password_here"
)


pgvdb$upload_to_higlass(
    endpoint,
    datafile = system.file("extdata", "test_data", "higlass_test_bigwig.bw", package = "Skilift"),
    name = "test_bigwig",
    filetype = "bigwig",
    datatype = "vector",
    coordSystem = "hg38",
    username = "username_here",
    password = "password_here"
)
```

### delete_from_higlass

```r 
 pgvdb$delete_from_higlass(endpoint, uuid, username, password)
```

Delete a file from the higlass server. Will also remove the plot from the pgvdb
instance. Tilesets can only be deleted by their uuid.
