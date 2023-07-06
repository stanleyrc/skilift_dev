library(gGnome)
library(dplyr)
source("./PGVdb.R")
source("./functions.R")

db <- return_PGV_db(
  datafiles.json = "~/projects/pgv/public/datafiles.json",
  data_folder = "~/projects/pgv/public/data",
  PGV_public_dir = "~/projects/pgv/public"
)

reset_db <- function() {
  datafiles <- "~/projects/pgv/public/datafiles.json"
  reset_datafiles <- "~/projects/pgv/public/datafiles0.json"
  file.copy(reset_datafiles, datafiles)
}

test_add_tags <- function() {
  add_tag <- data.table(patient.id = "DEMO", tags = "random tag")
  db$descriptions <- rbind(db$descriptions, add_tag)
  db$descriptions[patient.id = "DEMO"]
  db$descriptions <- db$descriptions[-c(11), ]
  push_PGV_db(db)
}

# Create a new patient called TEST which is an exact duplicate of DEMO
test_create_patient <- function() {
  id <- "TEST"
  new_graph <- data.table(
    patient.id = id,
    name.col = id,
    gg.col = "HCC1954.gg.rds",
    gw.col = NA,
    annotation = list(c(
      "simple", "bfb", "chromoplexy",
      "chromothripsis", "del", "dm", "dup",
      "pyrgo", "qrdel", "qrdup", "qrp", "rigma",
      "tic", "tyfonas"
    )),
    tree = NA,
    ref = "hg19",
    cov.field = NA,
    dirpaths = file.path(paste0(db$data_folder, "/", id))
  )

  add_patients_PGV(db, table_add = new_graph, cores = 1)
}

test_delete_patient <- function() {
  id <- "TEST"
  drop_patients_PGV(db, "DEMO", delete = TRUE)
}

test_create_patient()
test_delete_patient()
