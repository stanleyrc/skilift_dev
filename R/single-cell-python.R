#' Run a gOS single-cell Python exporter shipped in inst/python/gos_single_cell
#'
#' @param name exporter script name, with or without the .py suffix (e.g. "gos_sc_tert_sv")
#' @param args character vector of command-line arguments
#' @param python python executable (env GOS_PYTHON, default python3)
#' @param ... passed to system2 (e.g. stdout = TRUE)
#' @return the system2 result (exit status, or output when stdout = TRUE)
#' @export
sc_run_python_exporter <- function(name, args = character(), python = Sys.getenv("GOS_PYTHON", "python3"), ...) {
  if (!grepl("\\.py$", name)) name <- paste0(name, ".py")
  script <- system.file("python", "gos_single_cell", name, package = "Skilift")
  if (!nzchar(script) || !file.exists(script)) {
    stop("python exporter not found in inst/python/gos_single_cell: ", name)
  }
  system2(python, c(shQuote(script), args), ...)
}
