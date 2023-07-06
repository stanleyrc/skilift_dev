library(gGnome)
library(dplyr)
source("./PGVdb.R")
source("./functions.R")

db = return_PGV_db(
  datafiles.json = "/gpfs/commons/groups/imielinski_lab/pgv_content/datafiles.json", 
  data_folder = "/gpfs/commons/groups/imielinski_lab/pgv_content/data/", 
  PGV_public_dir = "/gpfs/commons/home/cxanthopoulakis/pgv/public/"
)

add_tag = data.table(patient.id = "DEMO", tags = "random tag")
db$descriptions = rbind(db$descriptions, add_tag)
db$descriptions[patient.id == "DEMO"]
db$descriptions = db$descriptions[-c(11),]
# push_PGV_db(db)

# Create a new patient called TEST which is an exact duplicate of DEMO
# Add TEST to the pgv dataset
