cli_args <- commandArgs(trailingOnly = TRUE)
smoke <- "--smoke" %in% cli_args
full <- "--full" %in% cli_args

rscript <- function() {
  exe <- if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  file.path(R.home("bin"), exe)
}

ensure_dirs <- function() {
  dirs <- c("data/processed", "results/figures", "logs")
  for (d in dirs) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }
}

cleanup_rplots <- function() {
  while (!identical(names(grDevices::dev.cur()), "null device")) {
    try(grDevices::dev.off(), silent = TRUE)
  }
  if (file.exists("Rplots.pdf")) unlink("Rplots.pdf", force = TRUE)
}

require_files <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required file(s):\n  ", paste(missing, collapse = "\n  "), call. = FALSE)
  }
}

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop(
      "Missing required R package(s): ",
      paste(missing, collapse = ", "),
      "\nInstall them before running the reproduction workflow.",
      call. = FALSE
    )
  }
}

require_system_tool <- function(tool) {
  if (!nzchar(Sys.which(tool))) {
    stop("Missing required system utility: ", tool, call. = FALSE)
  }
}

run_step <- function(script, extra_args = character()) {
  require_files(script)
  cmd_args <- c(script, extra_args)
  message("Running: Rscript ", paste(cmd_args, collapse = " "))
  status <- system2(rscript(), cmd_args)
  cleanup_rplots()
  if (!identical(status, 0L)) stop("Step failed: ", script, call. = FALSE)
}

run_smoke_test <- function() {
  message("Running quick code smoke test.")
  suppressPackageStartupMessages(library(GeoModels))
  env <- new.env(parent = globalenv())
  sys.source("src/simulation_config.R", envir = env)
  sys.source("src/simulation_utils.R", envir = env)

  if (!is.data.frame(env$simulation_scenarios) || nrow(env$simulation_scenarios) != 9L) {
    stop("Smoke test failed: simulation scenarios did not load as expected.", call. = FALSE)
  }
  if (!is.function(env$parse_cli_args) || !is.function(env$write_design)) {
    stop("Smoke test failed: simulation utility functions did not load.", call. = FALSE)
  }
  if (!identical(env$simulation_scenarios$nu, c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 4L))) {
    stop("Smoke test failed: simulation scenario nu values changed unexpectedly.", call. = FALSE)
  }
  message("Smoke test passed.")
}

ensure_dirs()
cleanup_rplots()

base_packages <- c("GeoModels", "fields", "MASS", "ggplot2")
full_figure_packages <- c("cubature", "hypergeo", "mvtnorm", "gcKrig", "plot3D")
require_packages(if (smoke) base_packages else unique(c(base_packages, full_figure_packages)))
require_system_tool("pdfseparate")

require_files(c(
  "data/raw/LAS_sec_148.csv",
  "data/processed/simulation_raw_results.csv",
  "data/processed/coverage_raw_results.csv",
  "data/processed/weibull_clayton_copula_simulations.csv",
  "data/processed/weibull_gaussian_copula_simulations.csv",
  "data/processed/weibull_direct_simulations.csv"
))

if (smoke) {
  run_smoke_test()
  cleanup_rplots()
  message("Smoke workflow completed without regenerating outputs.")
  quit(save = "no", status = 0)
}

if (!full) {
  run_smoke_test()
  cleanup_rplots()
  message("No outputs were regenerated. Use --full only when intentionally rebuilding results and figures.")
  quit(save = "no", status = 0)
}

run_step("src/fit_spatial_models.R")
run_step("src/bootstrap_standard_errors.R", "--filter-only")
run_step("src/summarize_simulation_study.R")

figure_scripts <- c(
  "src/generate_diagnostics.R",
  "src/generate_combined_variogram.R",
  "src/generate_boxplots.R"
)

if (!smoke) {
  figure_scripts <- c(
    figure_scripts,
    "src/generate_contours.R",
    "src/generate_realizations.R",
    "src/generate_correlation_function.R"
  )
}

for (script in figure_scripts) {
  run_step(script)
}

cleanup_rplots()
message("Reproduction workflow completed.")
