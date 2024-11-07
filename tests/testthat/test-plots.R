suppressWarnings(devtools::load_all())


library(testthat)

# test <- function() { testthat::test_file("tests/testthat/test-plots.R") }

setup({
  ot_test_paths <<- list(
    oncotable = system.file('extdata/test_data/oncotable_test_data/new_oncotable/oncotable.rds', package='Skilift'),
    unit_oncotable = system.file('extdata/test_data/oncotable_test_data/new_oncotable/unit_oncotable.rds', package='Skilift'),
    jabba_simple_gg = system.file('extdata/test_data/oncotable_test_data/jabba.simple.gg.rds', package='Skilift'),
  )
})

