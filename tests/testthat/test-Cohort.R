suppressWarnings(devtools::load_all())

library(testthat)

setup({
  test_pairs <<- system.file("extdata/test_data", "test_pairs.rds", package = "Skilift")
})

op = "~/projects/Clinical_NYU/db/pairs.rds"
