simulation_args_v2_preserved <- if (exists("simulation_args_v2", envir = .GlobalEnv, inherits = FALSE)) {
  get("simulation_args_v2", envir = .GlobalEnv)
} else {
  NULL
}
rm(list = setdiff(ls(), "simulation_args_v2_preserved"))
if (!is.null(simulation_args_v2_preserved)) simulation_args_v2 <- simulation_args_v2_preserved

source("src/simulation_config_v2.R")
source("src/simulation_utils_v2.R")

args <- parse_cli_args_v2()
sim_reps_v2 <- if ("sim_reps" %in% args$provided_args || !("reps" %in% args$provided_args)) args$sim_reps else args$reps
coverage_reps_v2 <- if ("coverage_reps" %in% args$provided_args || !("reps" %in% args$provided_args)) args$coverage_reps else args$reps
sim_chunk_size_v2 <- if ("sim_chunk_size" %in% args$provided_args || !("chunk_size" %in% args$provided_args)) args$sim_chunk_size else args$chunk_size
coverage_chunk_size_v2 <- if ("coverage_chunk_size" %in% args$provided_args || !("chunk_size" %in% args$provided_args)) args$coverage_chunk_size else args$chunk_size

rscript_v2 <- function() {
  exe <- if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  file.path(R.home("bin"), exe)
}

run_step_v2 <- function(script, extra_args = character()) {
  cmd_args <- c(script, extra_args)
  message("Running: Rscript ", paste(cmd_args, collapse = " "))
  status <- system2(rscript_v2(), cmd_args)
  if (!identical(status, 0L)) stop("Step failed: ", script)
}

if (isTRUE(args$smoke)) {
  run_step_v2("src/run_simulation_study_v2.R", c("--smoke", "--overwrite"))
  run_step_v2("src/run_coverage_study_v2.R", c("--smoke", "--overwrite"))
  run_step_v2("src/summarize_simulation_study_v2.R")
} else {
  sim_args <- c(
    "--reps", as.character(sim_reps_v2),
    "--workers", as.character(args$workers),
    "--chunk-size", as.character(sim_chunk_size_v2),
    "--scenario", as.character(args$scenario),
    "--verbose", as.character(args$verbose)
  )
  cov_args <- c(
    "--reps", as.character(coverage_reps_v2),
    "--boot", as.character(args$boot),
    "--workers", as.character(args$workers),
    "--chunk-size", as.character(coverage_chunk_size_v2),
    "--scenario", ifelse(is.null(args$coverage_scenario), "S0", as.character(args$coverage_scenario)),
    "--m-values", as.character(args$m_values),
    "--verbose", as.character(args$verbose)
  )
  if (isTRUE(args$overwrite)) {
    sim_args <- c(sim_args, "--overwrite")
    cov_args <- c(cov_args, "--overwrite")
  }
  run_step_v2("src/run_simulation_study_v2.R", sim_args)
  run_step_v2("src/run_coverage_study_v2.R", cov_args)
  run_step_v2("src/summarize_simulation_study_v2.R")
}

message("All v2 steps completed.")
