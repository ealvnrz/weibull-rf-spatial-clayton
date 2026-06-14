simulation_args_preserved <- if (exists("simulation_args", envir = .GlobalEnv, inherits = FALSE)) {
  get("simulation_args", envir = .GlobalEnv)
} else {
  NULL
}
rm(list = setdiff(ls(), "simulation_args_preserved"))
if (!is.null(simulation_args_preserved)) simulation_args <- simulation_args_preserved

source("src/simulation_config.R")
source("src/simulation_utils.R")

args <- parse_cli_args()
mode <- ifelse(is.null(args$mode), "all", tolower(as.character(args$mode)))
if (!mode %in% c("all", "simulation", "coverage")) {
  stop("--mode must be one of: all, simulation, coverage")
}
sim_reps <- if ("sim_reps" %in% args$provided_args || !("reps" %in% args$provided_args)) args$sim_reps else args$reps
coverage_reps <- if ("coverage_reps" %in% args$provided_args || !("reps" %in% args$provided_args)) args$coverage_reps else args$reps
sim_chunk_size <- if ("sim_chunk_size" %in% args$provided_args || !("chunk_size" %in% args$provided_args)) args$sim_chunk_size else args$chunk_size
coverage_chunk_size <- if ("coverage_chunk_size" %in% args$provided_args || !("chunk_size" %in% args$provided_args)) args$coverage_chunk_size else args$chunk_size

max_procs <- ifelse(is.null(args$max_procs), 3L, as.integer(args$max_procs))
max_procs <- max(1L, max_procs)
poll_seconds <- ifelse(is.null(args$poll_seconds), 10L, as.integer(args$poll_seconds))
poll_seconds <- max(1L, poll_seconds)

default_sim_groups <- "S0,S1,S2;S3,S4,S5;S6,S7,S8"
sim_groups_arg <- ifelse(is.null(args$sim_groups), default_sim_groups, as.character(args$sim_groups))
sim_groups <- trimws(strsplit(sim_groups_arg, ";", fixed = TRUE)[[1]])
sim_groups <- sim_groups[nzchar(sim_groups)]

m_values <- as.integer(trimws(strsplit(args$m_values, ",", fixed = TRUE)[[1]]))

ensure_output_dirs()
if (!dir.exists("logs")) dir.create("logs", recursive = TRUE)

rscript <- function() {
  exe <- if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
  file.path(R.home("bin"), exe)
}

timestamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")
run_tag <- ifelse(is.null(args$tag), timestamp(), as.character(args$tag))

job_table <- data.frame(
  job_id = character(),
  kind = character(),
  value = character(),
  script = character(),
  stdout = character(),
  stderr = character(),
  stringsAsFactors = FALSE
)

add_job <- function(job_id, kind, value, script, script_args) {
  stdout <- file.path("logs", paste0(run_tag, "_", job_id, ".log"))
  stderr <- file.path("logs", paste0(run_tag, "_", job_id, ".err"))
  list(
    job_id = job_id,
    kind = kind,
    value = value,
    script = script,
    args = c(script, script_args),
    stdout = stdout,
    stderr = stderr
  )
}

common_args <- function(reps, chunk_size) {
  out <- c(
    "--reps", as.character(reps),
    "--workers", "1",
    "--chunk-size", as.character(chunk_size),
    "--verbose", as.character(args$verbose)
  )
  if (isTRUE(args$overwrite)) out <- c(out, "--overwrite")
  out
}

jobs <- list()

if (mode %in% c("all", "simulation")) {
  for (i in seq_along(sim_groups)) {
    group <- sim_groups[[i]]
    jobs[[length(jobs) + 1L]] <- add_job(
      job_id = paste0("simulation_group_", i),
      kind = "simulation",
      value = group,
      script = "src/run_simulation_study.R",
      script_args = c(common_args(sim_reps, sim_chunk_size), "--scenario", group)
    )
  }
}

if (mode %in% c("all", "coverage")) {
  for (m in m_values) {
    jobs[[length(jobs) + 1L]] <- add_job(
      job_id = paste0("coverage_m", m),
      kind = "coverage",
      value = as.character(m),
      script = "src/run_coverage_study.R",
      script_args = c(
        common_args(coverage_reps, coverage_chunk_size),
        "--scenario", "S0",
        "--boot", as.character(args$boot),
        "--m-values", as.character(m)
      )
    )
  }
}

if (length(jobs) == 0L) stop("No jobs to run.")

job_table <- do.call(rbind, lapply(jobs, function(job) {
  data.frame(
    job_id = job$job_id,
    kind = job$kind,
    value = job$value,
    script = job$script,
    stdout = job$stdout,
    stderr = job$stderr,
    stringsAsFactors = FALSE
  )
}))
write.csv(job_table, file.path("logs", paste0(run_tag, "_external_jobs.csv")), row.names = FALSE)

if (!requireNamespace("processx", quietly = TRUE)) {
  stop("Package 'processx' is required for external parallel execution.")
}

message("External parallel launcher")
message("Mode: ", mode)
message("Simulation replicates: ", sim_reps)
message("Coverage replicates: ", coverage_reps)
message("Coverage bootstrap replicates: ", args$boot)
message("Jobs: ", length(jobs))
message("Max concurrent processes: ", max_procs)
message("Run tag: ", run_tag)
message("Logs directory: ", normalizePath("logs", winslash = "/", mustWork = TRUE))

pending <- jobs
running <- list()
completed <- list()
failed <- list()

start_job <- function(job) {
  message(format(Sys.time(), "%H:%M:%S"), " | starting ", job$job_id, " (", job$kind, ": ", job$value, ")")
  processx::process$new(
    command = rscript(),
    args = job$args,
    stdout = job$stdout,
    stderr = job$stderr,
    supervise = TRUE
  )
}

while (length(pending) > 0L || length(running) > 0L) {
  while (length(pending) > 0L && length(running) < max_procs) {
    job <- pending[[1L]]
    pending <- pending[-1L]
    proc <- start_job(job)
    running[[job$job_id]] <- list(job = job, proc = proc, started = Sys.time())
  }

  Sys.sleep(poll_seconds)

  if (length(running) == 0L) next
  for (job_id in names(running)) {
    item <- running[[job_id]]
    if (!item$proc$is_alive()) {
      status <- item$proc$get_exit_status()
      elapsed <- round(as.numeric(difftime(Sys.time(), item$started, units = "mins")), 2)
      message(format(Sys.time(), "%H:%M:%S"), " | finished ", job_id, " | status=", status, " | elapsed=", elapsed, " min")
      if (identical(status, 0L)) {
        completed[[job_id]] <- item
      } else {
        failed[[job_id]] <- item
      }
      running[[job_id]] <- NULL
    }
  }
}

message("Completed jobs: ", length(completed))
message("Failed jobs: ", length(failed))

if (length(failed) > 0L) {
  failed_ids <- names(failed)
  message("Failed job logs:")
  for (job_id in failed_ids) {
    message("  ", job_id, ": ", failed[[job_id]]$job$stderr)
  }
  stop("At least one external job failed.")
}

message("Combining chunk outputs and regenerating figures.")
status <- system2(rscript(), "src/summarize_simulation_study.R")
if (!identical(status, 0L)) stop("Summary step failed.")

message("External parallel launcher complete.")
