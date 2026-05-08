simulation_args_v2_preserved <- if (exists("simulation_args_v2", envir = .GlobalEnv, inherits = FALSE)) {
  get("simulation_args_v2", envir = .GlobalEnv)
} else {
  NULL
}
rm(list = setdiff(ls(), "simulation_args_v2_preserved"))
if (!is.null(simulation_args_v2_preserved)) simulation_args_v2 <- simulation_args_v2_preserved

source("src/simulation_config_v2.R")
source("src/simulation_utils_v2.R")

ensure_output_dirs_v2()

simulation_raw <- read_chunk_dir_v2("data/processed/simulation_chunks_v2")
if (nrow(simulation_raw) > 0L) {
  simulation_raw <- dedupe_simulation_raw_v2(simulation_raw)
  write.csv(simulation_raw, "data/processed/simulation_raw_results_v2.csv", row.names = FALSE)
  recovery <- parameter_recovery_summary_v2(simulation_raw)
  selection <- model_selection_summary_v2(simulation_raw)
  write.csv(recovery, "data/processed/simulation_parameter_recovery_v2.csv", row.names = FALSE)
  write.csv(selection, "data/processed/simulation_model_selection_v2.csv", row.names = FALSE)
  plot_selection_frequencies_v2(selection, "results/figures/simulation_selection_frequencies_v2.pdf")
  plot_parameter_recovery_v2(recovery, "results/figures/simulation_parameter_recovery_v2.pdf")
  message("Simulation summaries and figures written.")
} else {
  warning("No simulation chunk files found.")
}

coverage_raw <- read_chunk_dir_v2("data/processed/coverage_chunks_v2")
if (nrow(coverage_raw) > 0L) {
  coverage_raw <- dedupe_coverage_raw_v2(coverage_raw)
  write.csv(coverage_raw, "data/processed/coverage_raw_results_v2.csv", row.names = FALSE)
  coverage <- coverage_summary_v2(coverage_raw)
  write.csv(coverage, "data/processed/coverage_summary_v2.csv", row.names = FALSE)
  write.csv(coverage, "data/processed/neighborhood_sensitivity_v2.csv", row.names = FALSE)
  plot_coverage_calibration_v2(coverage, "results/figures/coverage_calibration_v2.pdf")
  message("Coverage summaries and figures written.")
} else {
  warning("No coverage chunk files found.")
}

cleanup_rplots_v2()
