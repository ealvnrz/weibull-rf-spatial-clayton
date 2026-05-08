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
if (!("reps" %in% args$provided_args) && !is.null(args$coverage_reps)) {
  args$reps <- args$coverage_reps
}
if (!("chunk_size" %in% args$provided_args) && !is.null(args$coverage_chunk_size)) {
  args$chunk_size <- args$coverage_chunk_size
}
if (isTRUE(args$smoke)) {
  args$reps <- 2L
  args$boot <- 2L
  args$chunk_size <- 2L
  args$workers <- 1L
  args$scenario <- "S0"
  args$m_values <- "2"
}
if (is.null(args$scenario) || identical(tolower(as.character(args$scenario)), "all")) {
  args$scenario <- "S0"
}

ensure_output_dirs_v2()
write_design_v2()

scenarios <- filter_scenarios_v2(args$scenario, FALSE)
if (nrow(scenarios) != 1L) {
  stop("Coverage study expects exactly one scenario. Use --scenario S0 unless intentionally changing it.")
}
m_values <- as.integer(trimws(strsplit(args$m_values, ",", fixed = TRUE)[[1]]))
chunks <- make_chunks_v2(args$reps, args$chunk_size)
if (isTRUE(args$overwrite)) {
  clear_coverage_chunks_v2(scenarios$scenario_id[1], m_values)
}

tasks <- do.call(rbind, lapply(m_values, function(m) {
  data.frame(
    scenario_id = scenarios$scenario_id[1],
    rep_start = chunks$rep_start,
    rep_end = chunks$rep_end,
    neighb = m,
    boot = args$boot,
    seed_base = args$seed_base,
    overwrite = args$overwrite,
    verbose = args$verbose,
    stringsAsFactors = FALSE
  )
}))

message("Coverage study v2")
message("Scenario: ", scenarios$scenario_id[1])
message("Replicates: ", args$reps)
message("Bootstrap replicates per fit: ", args$boot)
message("Neighborhood values: ", paste(m_values, collapse = ", "))

raw <- run_chunks_parallel_v2(tasks, args$workers, kind = "coverage")
raw <- dedupe_coverage_raw_v2(raw)
write.csv(raw, "data/processed/coverage_raw_results_v2.csv", row.names = FALSE)

coverage <- coverage_summary_v2(raw)
write.csv(coverage, "data/processed/coverage_summary_v2.csv", row.names = FALSE)
write.csv(coverage, "data/processed/neighborhood_sensitivity_v2.csv", row.names = FALSE)

message("Coverage study complete.")
message("Wrote data/processed/coverage_raw_results_v2.csv")
message("Wrote data/processed/coverage_summary_v2.csv")
message("Wrote data/processed/neighborhood_sensitivity_v2.csv")
cleanup_rplots_v2()
