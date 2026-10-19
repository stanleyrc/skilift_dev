test_that("sc_run_python_exporter runs an exporter's --help", {
  skip_if(Sys.which(Sys.getenv("GOS_PYTHON", "python3")) == "")
  out <- suppressWarnings(sc_run_python_exporter("gos_sc_tert_sv", "--help", stdout = TRUE, stderr = TRUE))
  expect_true(any(grepl("usage", out, ignore.case = TRUE)))
  expect_error(sc_run_python_exporter("no_such_exporter"), "not found")
})
