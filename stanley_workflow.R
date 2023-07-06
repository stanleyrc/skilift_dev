library(gGnome)
library(dplyr)
source("./PGVdb.R")
source("./functions.R")

worked_notworked = lapply(atac.pairs$Patient_ID, function(pair) {
  file.exists(filename=paste0("~/lab/pgv_content/data/", pair, "/peel.json"))
})

pairs.ran <- atac.pairs$Patient_ID[unlist(worked_notworked)!="FALSE"]

# Add peels
PGVdb = getPGV()

#create dt to add
dt.add = data.table(
  sample = pairs.ran,  type = "walk",  visible = TRUE,  patient.id = pairs.ran,  
  plot_id = paste0(pairs.ran,  "_8"),  source = "peel.json",  
  title = "Peel"
)

# rbind and push
PGVdb$plots
PGVdb$plots <- rbind(PGVdb$plots, dt.add, fill=TRUE)
pushDB(PGVdb)

###

# Centromere gGraphs
atac.pairs = readRDS("~/projects/gGnome/db/atac.pairs.rds")
hg38.bands = fread("~/DB/UCSC/hg38.cytoband.txt")
colnames(hg38.bands) = c("seqnames", "start", "end", "name", "type")
hg38.bands = hg38.bands[type=="acen", ]

# adding all patients in atac pairs
identatac.dt <- as.data.table(t(as.data.table(strsplit(atac.pairs$pair, "-"))))
identatac.dt <- identatac.dt[, 1:3]
ids.pairs <- paste0(identatac.dt$V1, "-", identatac.dt$V2, "-", identatac.dt$V3)
atac.pairs[, patient.id_all := ids.pairs]
samples = atac.pairs$patient.id_all %>% unique()
mclapply(samples, function(pair){
  dt2json(dt=hg38.bands, filename=paste0("~/lab/pgv_content/data/", pair, "/hg38_cent.json"))
}, mc.cores=30)

# add hits to pgvdb
PGVdb = getPGV()

# already added TCGA-06-A5U0 when testing so adding everything else
# create dt to add
dt.add = data.table(sample = samples, type = "genome", visible = TRUE,  patient.id = samples,  plot_id = paste0(samples, "_7"), source = "hg38_cent.json",  title = "Centromeres")

#rbind and push
PGVdb$plots <- rbind(PGVdb$plots, dt.add, fill=TRUE)
PGVdb$plots
push_PGV_db(PGVdb)
