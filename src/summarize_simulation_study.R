simulation_args_preserved <- if (exists("simulation_args", envir = .GlobalEnv, inherits = FALSE)) {
  get("simulation_args", envir = .GlobalEnv)
} else {
  NULL
}
rm(list = setdiff(ls(), "simulation_args_preserved"))
if (!is.null(simulation_args_preserved)) simulation_args <- simulation_args_preserved

source("src/simulation_config.R")
source("src/simulation_utils.R")

ensure_output_dirs()

read_chunks_or_raw <- function(chunk_dir, raw_file, dedupe_fun) {
  raw <- read_chunk_dir(chunk_dir)
  if (nrow(raw) > 0L) {
    return(dedupe_fun(raw))
  }
  if (file.exists(raw_file)) {
    message("No chunk files found in ", chunk_dir, "; using precomputed ", raw_file, ".")
    return(dedupe_fun(read.csv(raw_file, stringsAsFactors = FALSE)))
  }
  raw
}

simulation_raw <- read_chunks_or_raw(
  "data/processed/simulation_chunks",
  "data/processed/simulation_raw_results.csv",
  dedupe_simulation_raw
)
if (nrow(simulation_raw) > 0L) {
  write.csv(simulation_raw, "data/processed/simulation_raw_results.csv", row.names = FALSE)
  recovery <- parameter_recovery_summary(simulation_raw)
  selection <- model_selection_summary(simulation_raw)
  write.csv(recovery, "data/processed/simulation_parameter_recovery.csv", row.names = FALSE)
  write.csv(selection, "data/processed/simulation_model_selection.csv", row.names = FALSE)
  plot_selection_frequencies(selection, "results/figures/simulation_selection_frequencies.pdf")
  plot_parameter_recovery(recovery, "results/figures/simulation_parameter_recovery.pdf")
  message("Simulation summaries and figures written.")
} else {
  warning("No simulation chunk files found.")
}

coverage_raw <- read_chunks_or_raw(
  "data/processed/coverage_chunks",
  "data/processed/coverage_raw_results.csv",
  dedupe_coverage_raw
)
if (nrow(coverage_raw) > 0L) {
  write.csv(coverage_raw, "data/processed/coverage_raw_results.csv", row.names = FALSE)
  coverage <- coverage_summary(coverage_raw)
  write.csv(coverage, "data/processed/coverage_summary.csv", row.names = FALSE)
  write.csv(coverage, "data/processed/neighborhood_sensitivity.csv", row.names = FALSE)
  plot_coverage_calibration(coverage, "results/figures/coverage_calibration.pdf")
  message("Coverage summaries and figures written.")
} else {
  warning("No coverage chunk files found.")
}

cleanup_rplots()
