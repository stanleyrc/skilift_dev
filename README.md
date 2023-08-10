# PGVdb
Faciliates easier loading and manipulation of genomic data for use with [PGV](https://github.com/mskilab-org/pgv).

# Installation
- To install directly from github: `devtools::install_github("mskilab-org/pgvdb")`
- To install from source:
    1. `git clone https://github.com/mskilab-org/pgvdb.git`
    2. `devtools::load_all("path/to/clone")`

# Usage
The PGVdb class manages a database of patient metadata and genomic plot data. It provides an interface to load, update, query, and manipulate the database programmatically.

The key methods allow converting to/from JSON for data storage, validating the data, adding/removing plots, and generating plot JSON for visualization. The metadata and plots are stored as data tables for easy manipulation.

## Methods
- `initialize(datafiles_json_path, datadir, publicdir, settings)`: Initialize a new PGVdb object

- `update_datafiles_json()`: Update JSON data files on disk 

- `to_datatable(filter)`: Convert to a single data table, apply optional filter

- `add_plots(new_plots_dt)`: Add new plots 

- `remove_plots(plots_to_remove_dt)`: Remove plots

- `validate()`: Validate metadata and plots

## Fields

- `metadata`: Data table containing patient metadata 

- `plots`: Data table containing plot metadata

- `datadir`: Path to data directory

- `publicdir`: Path to public directory 

- `settings`: Path to settings file

### Constructor
```r
pgv <- PGVdb$new(datafiles_json_path, datadir, publicdir, settings)
```

Create a new PGVdb object by passing the path to the datafiles.json, data directory, public directory, and settings file.

### update_datafiles_json
```r
pgv$update_datafiles_json()  
```

Update the datafiles JSON on disk with the current metadata and plots.

### to_datatable
```r
dt <- pgv$to_datatable(filter = c("patient.id", "HCC"))
```

Convert metadata and plots to a single data table, applying an optional filter.

### add_plots

```r
pgv$add_plots(plots_to_add, overwrite = FALSE, cores=1)  
```

Add new plots by passing a `data.table` containing required columns:

- `patient.id`: Patient identifier
- `x`: can be
    - list(list(server="", uuid="")): a list of the [server, uuid] 
    - list(GRanges), list(gWalk), list(gGraph): a list containing an object
    - Filepath to an RDS object
- `visible`: Whether plot is visible

The `overwrite` flag determines whether to overwrite existing source files. 

If required columns are missing, an empty table will be returned with the
required column names. The `type` column (indiciating the type of plot:
scatterplot, genome, walk, bigwig, etc) is derived from `x` if not supplied by
the user. The `source` column (indicating the name of the plot file inside the
pgv data directory) is derived from the `type`, if not supplied by the user.
Unless the overwrite flag is set, it will not overwrite existing plot files,
instead it will just increment the filename of the new plot file (i.e
coverage.json -> coverage2.json)

The `cores` parameter determines how many cores to use for parallel execution.
By default, it will not run in parallel (i.e use 1 core).

### remove_plots

```r
pgv$remove_plots(plots_to_remove, delete = FALSE)
```

Remove plots by passing a `data.table` with either:

- `patient.id`: Remove all plots for a patient 
- `patient.id` and `source`: Remove specific plots

Or alternatively with `patient.id`, `server`, and `uuid`.

The `delete` flag determines whether to also delete the source data files.

Passing just `patient.id` will remove all plots for that patient.

### validate

```r
pgv$validate()
```

Validate metadata and plot data, removing invalid entries. Is automatically
called when adding/removing plots. Useful if you make manual changes to the
pgvdb plots or metadata.


### init_pgv

```r 
pgv$init_pgv(pgv_dir, build=FALSE)
```

Initialize a pgv instance loaded with the data in pgvdb at `pgv_dir`. This
method will clone the ![pgv repo](https://github.com/mskilab-org/pgv) and
create symlinks in the pgv data directory that point to your pgvdb data. 

If the `build` flag is set to `TRUE` it will build pgv instead of launching a
local instance (useful when running on a remote server or hpc). In that case,
you should set `pgv_dir` to be a directory inside whichever directory is served
by your remote server (e.g `public_html`).
