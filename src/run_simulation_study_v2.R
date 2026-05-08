simulation_args_v2_preserved <- if (exists("simulation_args_v2", envir = .GlobalEnv, inherits = FALSE)) {
  get("simulation_args_v2", envir = .GlobalEnv)
} else {
  NULL
}
rm(list = setdiff(ls(), "simulation_args_v2_preserved"))
if (!is.null(simulation_args_v2_preserved)) simulation_args_v2 <- simulation_args_v2_preserved

source("src/simulation_config_v2.R")
source("src/simulation_utils_v2.R")

suppressPackageStartupMessages(library(GeoModels))

args <- parse_cli_args_v2()
if (!("reps" %in% args$provided_args) && !is.null(args$sim_reps)) {
  args$reps <- args$sim_reps
}
if (!("chunk_size" %in% args$provided_args) && !is.null(args$sim_chunk_size)) {
  args$chunk_size <- args$sim_chunk_size
}
if (isTRUE(args$smoke)) {
  args$reps <- 2L
  args$chunk_size <- 2L
  args$workers <- 1L
}

ensure_output_dirs_v2()
write_design_v2()

scenarios <- filter_scenarios_v2(args$scenario, args$smoke)
chunks <- make_chunks_v2(args$reps, args$chunk_size)
if (isTRUE(args$overwrite)) {
  clear_simulation_chunks_v2(scenarios$scenario_id)
}
tasks <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
  data.frame(
    scenario_id = scenarios$scenario_id[i],
    rep_start = chunks$rep_start,
    rep_end = chunks$rep_end,
    neighb = fixed_parameters_v2$neighb,
    seed_base = args$seed_base,
    overwrite = args$overwrite,
    verbose = args$verbose,
    stringsAsFactors = FALSE
  )
}))

message("Simulation study v2")
message("Scenarios: ", paste(scenarios$scenario_id, collapse = ", "))
message("Replicates per scenario: ", args$reps)
message("Chunk size: ", args$chunk_size)

raw <- run_chunks_parallel_v2(tasks, args$workers, kind = "simulation")
raw <- dedupe_simulation_raw_v2(raw)

write.csv(raw, "data/processed/simulation_raw_results_v2.csv", row.names = FALSE)
recovery <- parameter_recovery_summary_v2(raw)
selection <- model_selection_summary_v2(raw)
write.csv(recovery, "data/processed/simulation_parameter_recovery_v2.csv", row.names = FALSE)
write.csv(selection, "data/processed/simulation_model_selection_v2.csv", row.names = FALSE)

message("Simulation study complete.")
message("Wrote data/processed/simulation_raw_results_v2.csv")
message("Wrote data/processed/simulation_parameter_recovery_v2.csv")
message("Wrote data/processed/simulation_model_selection_v2.csv")
cleanup_rplots_v2()
