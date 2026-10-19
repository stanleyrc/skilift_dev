#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
opt <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option("--patient_dir"), optparse::make_option("--out_dir"),
    optparse::make_option("--work_dir"), optparse::make_option("--skilift"))))
devtools::load_all(opt$skilift, quiet = TRUE, compile = FALSE)
out <- sc_export_cell_signatures(opt$patient_dir, out_dir = opt$out_dir, work_dir = opt$work_dir)
cat("Cells:", length(out$cells), "too_few:", sum(vapply(out$cells, function(x) x$too_few, logical(1))), "\n")
