simulation_args_preserved <- if (exists("simulation_args", envir = .GlobalEnv, inherits = FALSE)) {
  get("simulation_args", envir = .GlobalEnv)
} else {
  NULL
}
rm(list = setdiff(ls(), "simulation_args_preserved"))
if (!is.null(simulation_args_preserved)) simulation_args <- simulation_args_preserved

source("src/simulation_config.R")
source("src/simulation_utils.R")

suppressPackageStartupMessages(library(GeoModels))

args <- parse_cli_args()
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

ensure_output_dirs()
write_design()

scenarios <- filter_scenarios(args$scenario, args$smoke)
chunks <- make_chunks(args$reps, args$chunk_size)
if (isTRUE(args$overwrite)) {
  clear_simulation_chunks(scenarios$scenario_id)
}
tasks <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(i) {
  data.frame(
    scenario_id = scenarios$scenario_id[i],
    rep_start = chunks$rep_start,
    rep_end = chunks$rep_end,
    neighb = fixed_parameters$neighb,
    seed_base = args$seed_base,
    overwrite = args$overwrite,
    verbose = args$verbose,
    stringsAsFactors = FALSE
  )
}))

message("Simulation study")
message("Scenarios: ", paste(scenarios$scenario_id, collapse = ", "))
message("Replicates per scenario: ", args$reps)
message("Chunk size: ", args$chunk_size)

raw <- run_chunks_parallel(tasks, args$workers, kind = "simulation")
raw <- dedupe_simulation_raw(raw)

write.csv(raw, "data/processed/simulation_raw_results.csv", row.names = FALSE)
recovery <- parameter_recovery_summary(raw)
selection <- model_selection_summary(raw)
write.csv(recovery, "data/processed/simulation_parameter_recovery.csv", row.names = FALSE)
write.csv(selection, "data/processed/simulation_model_selection.csv", row.names = FALSE)

message("Simulation study complete.")
message("Wrote data/processed/simulation_raw_results.csv")
message("Wrote data/processed/simulation_parameter_recovery.csv")
message("Wrote data/processed/simulation_model_selection.csv")
cleanup_rplots()
