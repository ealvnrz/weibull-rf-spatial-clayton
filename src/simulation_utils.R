# Utilities for the expanded simulation and coverage study.

parse_cli_args <- function(defaults = simulation_defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  if (exists("simulation_args", envir = .GlobalEnv, inherits = FALSE)) {
    source_args <- get("simulation_args", envir = .GlobalEnv)
    if (length(source_args) > 0L) args <- c(args, source_args)
  }
  out <- defaults
  provided_args <- character()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      i <- i + 1L
      next
    }
    key <- gsub("-", "_", sub("^--", "", key))
    value <- TRUE
    if (i < length(args) && !startsWith(args[[i + 1L]], "--")) {
      value <- args[[i + 1L]]
      i <- i + 1L
    }
    out[[key]] <- value
    provided_args <- c(provided_args, key)
    i <- i + 1L
  }

  if (!is.null(out$reps)) out$reps <- as.integer(out$reps)
  if (!is.null(out$sim_reps)) out$sim_reps <- as.integer(out$sim_reps)
  if (!is.null(out$coverage_reps)) out$coverage_reps <- as.integer(out$coverage_reps)
  if (!is.null(out$boot)) out$boot <- as.integer(out$boot)
  if (!is.null(out$chunk_size)) out$chunk_size <- as.integer(out$chunk_size)
  if (!is.null(out$sim_chunk_size)) out$sim_chunk_size <- as.integer(out$sim_chunk_size)
  if (!is.null(out$coverage_chunk_size)) out$coverage_chunk_size <- as.integer(out$coverage_chunk_size)
  if (!is.null(out$seed_base)) out$seed_base <- as.integer(out$seed_base)
  out$smoke <- isTRUE(out$smoke) || identical(tolower(as.character(out$smoke)), "true")
  out$overwrite <- isTRUE(out$overwrite) || identical(tolower(as.character(out$overwrite)), "true")
  out$verbose <- isTRUE(out$verbose) || identical(tolower(as.character(out$verbose)), "true")
  if (!is.null(out$quiet) && (isTRUE(out$quiet) || identical(tolower(as.character(out$quiet)), "true"))) {
    out$verbose <- FALSE
  }
  out$provided_args <- unique(provided_args)
  out
}

progress <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) {
    cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")
    flush.console()
  }
}

project_root <- function() {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

ensure_output_dirs <- function() {
  dirs <- c(
    "data/processed",
    "data/processed/simulation_chunks",
    "data/processed/coverage_chunks",
    "results/figures"
  )
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

resolve_workers <- function(workers, n_tasks = Inf) {
  if (identical(tolower(as.character(workers)), "auto")) {
    if (.Platform$OS.type == "windows") {
      message("Using 1 worker on Windows because GeoModels/progressr is not stable inside PSOCK workers.")
      workers <- 1L
    } else {
      workers <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
    }
  } else {
    workers <- as.integer(workers)
  }
  workers <- max(1L, workers)
  workers <- min(workers, max(1L, as.integer(n_tasks)))
  workers
}

filter_scenarios <- function(scenario_arg, smoke = FALSE) {
  if (isTRUE(smoke)) {
    return(simulation_scenarios[simulation_scenarios$scenario_id == "S0", , drop = FALSE])
  }
  if (is.null(scenario_arg) || identical(tolower(as.character(scenario_arg)), "all")) {
    return(simulation_scenarios)
  }
  ids <- trimws(strsplit(as.character(scenario_arg), ",", fixed = TRUE)[[1]])
  ans <- simulation_scenarios[simulation_scenarios$scenario_id %in% ids, , drop = FALSE]
  if (nrow(ans) == 0L) stop("No valid scenarios selected.")
  ans
}

scenario_index <- function(scenario_id) {
  match(scenario_id, simulation_scenarios$scenario_id)
}

scenario_by_id <- function(scenario_id) {
  ans <- simulation_scenarios[simulation_scenarios$scenario_id == scenario_id, , drop = FALSE]
  if (nrow(ans) != 1L) stop("Unknown scenario: ", scenario_id)
  ans[1, ]
}

write_design <- function(path = "data/processed/simulation_design.csv") {
  ensure_output_dirs()
  design <- simulation_scenarios
  design$beta0 <- fixed_parameters$mean
  design$beta1 <- fixed_parameters$mean1
  design$delta <- fixed_parameters$smooth
  design$mu <- 3.5
  design$power2 <- fixed_parameters$power2
  design$nugget <- fixed_parameters$nugget
  design$sill <- fixed_parameters$sill
  design$m <- fixed_parameters$neighb
  write.csv(design, path, row.names = FALSE)
}

make_chunks <- function(reps, chunk_size) {
  starts <- seq.int(1L, reps, by = chunk_size)
  data.frame(
    rep_start = starts,
    rep_end = pmin(starts + chunk_size - 1L, reps)
  )
}

create_design <- function(scenario, seed_base = simulation_defaults$seed_base) {
  set.seed(seed_base + scenario_index(scenario$scenario_id) * 1000L)
  n <- as.integer(scenario$n)
  coords <- cbind(runif(n), runif(n))
  X <- cbind(rep(1, n), runif(n))
  list(coords = coords, X = X)
}

true_params <- function(scenario) {
  c(
    mean = fixed_parameters$mean,
    mean1 = fixed_parameters$mean1,
    scale = as.numeric(scenario$alpha),
    shape = as.numeric(scenario$kappa)
  )
}

dgp_params <- function(scenario) {
  list(
    smooth = fixed_parameters$smooth,
    power2 = fixed_parameters$power2,
    mean = fixed_parameters$mean,
    mean1 = fixed_parameters$mean1,
    nu = as.numeric(scenario$nu),
    sill = fixed_parameters$sill,
    scale = as.numeric(scenario$alpha),
    nugget = fixed_parameters$nugget,
    shape = as.numeric(scenario$kappa)
  )
}

simulation_seed <- function(scenario_id, rep_id, seed_base = simulation_defaults$seed_base) {
  seed_base + scenario_index(scenario_id) * 100000L + as.integer(rep_id)
}

simulate_clayton_dataset <- function(scenario, coords, X, rep_id, seed_base) {
  set.seed(simulation_seed(scenario$scenario_id, rep_id, seed_base))
  GeoModels::GeoSimCopula(
    coordx = coords,
    corrmodel = fixed_parameters$corrmodel,
    model = fixed_parameters$model,
    param = dgp_params(scenario),
    X = X,
    copula = "Clayton"
  )[["data"]]
}

fit_args <- function(data, coords, X, scenario, model_id, neighb, sensitivity = FALSE) {
  true_param <- true_params(scenario)
  I <- fixed_parameters$upper_bound
  lower <- list(mean = -I, mean1 = -I, scale = 0.00001, shape = 0.00001)
  upper <- list(mean = I, mean1 = I, scale = I, shape = I)
  start <- list(
    mean = true_param[["mean"]],
    mean1 = true_param[["mean1"]],
    scale = true_param[["scale"]],
    shape = true_param[["shape"]]
  )

  fixed <- list(
    sill = fixed_parameters$sill,
    smooth = fixed_parameters$smooth,
    power2 = fixed_parameters$power2,
    nugget = fixed_parameters$nugget
  )
  copula <- NULL
  if (identical(model_id, "gaussian")) {
    copula <- "Gaussian"
  } else if (identical(model_id, "clayton")) {
    copula <- "Clayton"
    fixed$nu <- as.numeric(scenario$nu)
  }

  args <- list(
    data = data,
    coordx = coords,
    corrmodel = fixed_parameters$corrmodel,
    model = fixed_parameters$model,
    X = X,
    neighb = as.integer(neighb),
    likelihood = "Marginal",
    type = "Pairwise",
    optimizer = "nlminb",
    lower = lower,
    upper = upper,
    start = start,
    fixed = fixed,
    sensitivity = sensitivity
  )
  if (!is.null(copula)) args$copula <- copula
  args
}

fit_model_object <- function(data, coords, X, scenario, model_id, neighb, sensitivity = FALSE) {
  do.call(GeoModels::GeoFit, fit_args(data, coords, X, scenario, model_id, neighb, sensitivity))
}

param_or_na <- function(fit, name) {
  if (is.null(fit) || is.null(fit$param) || !(name %in% names(fit$param))) return(NA_real_)
  as.numeric(fit$param[[name]])
}

stderr_or_na <- function(fit, name) {
  if (is.null(fit) || is.null(fit$stderr) || !(name %in% names(fit$stderr))) return(NA_real_)
  as.numeric(fit$stderr[[name]])
}

fit_converged <- function(fit) {
  if (is.null(fit) || is.null(fit$param)) return(FALSE)
  vals <- as.numeric(unlist(fit$param))
  if (any(!is.finite(vals))) return(FALSE)
  conv <- fit$convergence
  if (is.null(conv)) return(TRUE)
  conv_text <- as.character(conv[1])
  if (conv_text %in% c("Successful", "0")) return(TRUE)
  conv_num <- suppressWarnings(as.numeric(conv_text))
  is.finite(conv_num) && conv_num == 0
}

fit_model_row <- function(data, coords, X, scenario, rep_id, model_id, neighb, verbose = TRUE) {
  model_label <- model_specs$model_label[match(model_id, model_specs$model_id)]
  true_param <- true_params(scenario)
  progress(
    sprintf(
      "%s rep %04d | fitting %s (m=%s)",
      scenario$scenario_id, as.integer(rep_id), model_label, as.integer(neighb)
    ),
    verbose = verbose
  )
  start_time <- proc.time()[["elapsed"]]
  fit <- NULL
  error_message <- NA_character_
  fit <- tryCatch(
    fit_model_object(data, coords, X, scenario, model_id, neighb, sensitivity = FALSE),
    error = function(e) {
      error_message <<- conditionMessage(e)
      NULL
    }
  )
  elapsed <- proc.time()[["elapsed"]] - start_time
  ok <- fit_converged(fit)
  progress(
    sprintf(
      "%s rep %04d | %s done | converged=%s | elapsed=%.1fs",
      scenario$scenario_id, as.integer(rep_id), model_label, ok, elapsed
    ),
    verbose = verbose
  )
  loglik <- if (!is.null(fit) && !is.null(fit$logCompLik)) as.numeric(fit$logCompLik) else NA_real_
  n_param <- length(parameter_names)
  claic_approx <- if (is.finite(loglik)) -2 * loglik + 2 * n_param else NA_real_
  clbic_approx <- if (is.finite(loglik)) -2 * loglik + log(as.numeric(scenario$n)) * n_param else NA_real_

  data.frame(
    scenario_id = scenario$scenario_id,
    scenario_group = scenario$scenario_group,
    rep = as.integer(rep_id),
    n = as.integer(scenario$n),
    alpha = as.numeric(scenario$alpha),
    kappa = as.numeric(scenario$kappa),
    nu = as.numeric(scenario$nu),
    fitted_model = model_id,
    fitted_model_label = model_label,
    neighb = as.integer(neighb),
    converged = ok,
    logCompLik = loglik,
    CLAIC_approx = claic_approx,
    CLBIC_approx = clbic_approx,
    elapsed_seconds = elapsed,
    estimate_mean = param_or_na(fit, "mean"),
    estimate_mean1 = param_or_na(fit, "mean1"),
    estimate_scale = param_or_na(fit, "scale"),
    estimate_shape = param_or_na(fit, "shape"),
    true_mean = true_param[["mean"]],
    true_mean1 = true_param[["mean1"]],
    true_scale = true_param[["scale"]],
    true_shape = true_param[["shape"]],
    error = ifelse(is.na(error_message), "", error_message),
    stringsAsFactors = FALSE
  )
}

run_simulation_replication <- function(scenario, rep_id, neighb, seed_base, verbose = TRUE) {
  progress(
    sprintf(
      "%s rep %04d | simulating Clayton-Weibull data (n=%d, alpha=%.3f, kappa=%.3f, nu=%.3f)",
      scenario$scenario_id, as.integer(rep_id), as.integer(scenario$n),
      as.numeric(scenario$alpha), as.numeric(scenario$kappa), as.numeric(scenario$nu)
    ),
    verbose = verbose
  )
  design <- create_design(scenario, seed_base)
  data <- simulate_clayton_dataset(scenario, design$coords, design$X, rep_id, seed_base)
  rows <- lapply(model_specs$model_id, function(model_id) {
    fit_model_row(data, design$coords, design$X, scenario, rep_id, model_id, neighb, verbose)
  })
  do.call(rbind, rows)
}

simulation_chunk_file <- function(task) {
  file.path(
    "data/processed/simulation_chunks",
    sprintf("%s_reps_%04d_%04d.csv", task$scenario_id, task$rep_start, task$rep_end)
  )
}

coverage_chunk_file <- function(task) {
  file.path(
    "data/processed/coverage_chunks",
    sprintf("%s_m%s_reps_%04d_%04d.csv", task$scenario_id, task$neighb, task$rep_start, task$rep_end)
  )
}

run_simulation_chunk <- function(task) {
  ensure_output_dirs()
  outfile <- simulation_chunk_file(task)
  if (file.exists(outfile) && !isTRUE(task$overwrite)) {
    progress(
      sprintf("Skipping existing simulation chunk %s", basename(outfile)),
      verbose = task$verbose
    )
    return(read.csv(outfile, stringsAsFactors = FALSE))
  }
  scenario <- scenario_by_id(task$scenario_id)
  progress(
    sprintf(
      "Starting simulation chunk %s | reps %04d-%04d | n=%d alpha=%.3f kappa=%.3f nu=%.3f",
      basename(outfile), task$rep_start, task$rep_end, as.integer(scenario$n),
      as.numeric(scenario$alpha), as.numeric(scenario$kappa), as.numeric(scenario$nu)
    ),
    verbose = task$verbose
  )
  rows <- lapply(seq.int(task$rep_start, task$rep_end), function(rep_id) {
    run_simulation_replication(scenario, rep_id, task$neighb, task$seed_base, task$verbose)
  })
  ans <- do.call(rbind, rows)
  write.csv(ans, outfile, row.names = FALSE)
  progress(
    sprintf("Finished simulation chunk %s | rows=%d", basename(outfile), nrow(ans)),
    verbose = task$verbose
  )
  ans
}

extract_fit_coords <- function(fit) {
  if (!is.null(fit$coordy) && length(fit$coordy) == length(fit$coordx)) {
    return(cbind(fit$coordx, fit$coordy))
  }
  if (is.matrix(fit$coordx) || is.data.frame(fit$coordx)) return(as.matrix(fit$coordx))
  fit$coordx
}

bootstrap_stderr <- function(fit, K, seed, verbose = TRUE, label = "") {
  if (!fit_converged(fit)) {
    out <- rep(NA_real_, length(parameter_names))
    names(out) <- parameter_names
    return(list(stderr = out, n_success = 0L))
  }

  coords <- extract_fit_coords(fit)
  X <- fit$X
  param_sim <- c(fit$param, fit$fixed)
  estimates <- NULL
  for (k in seq_len(K)) {
    if (k == 1L || k == K || k %% 10L == 0L) {
      progress(
        sprintf("%sbootstrap %03d/%03d", label, k, K),
        verbose = verbose
      )
    }
    set.seed(seed + k)
    sim_data <- tryCatch({
      if (is.null(fit$copula)) {
        GeoModels::GeoSim(
          coordx = coords,
          corrmodel = fit$corrmodel,
          model = fit$model,
          param = param_sim,
          X = X
        )[["data"]]
      } else {
        GeoModels::GeoSimCopula(
          coordx = coords,
          corrmodel = fit$corrmodel,
          model = fit$model,
          param = param_sim,
          X = X,
          copula = fit$copula
        )[["data"]]
      }
    }, error = function(e) NULL)
    if (is.null(sim_data)) next

    fit_args <- list(
      data = sim_data,
      coordx = coords,
      corrmodel = fit$corrmodel,
      model = fit$model,
      X = X,
      neighb = fit$neighb,
      likelihood = fit$likelihood,
      type = fit$type,
      optimizer = fit$optimizer,
      lower = fit$lower,
      upper = fit$upper,
      start = fit$param,
      fixed = fit$fixed,
      sensitivity = FALSE
    )
    if (!is.null(fit$copula)) fit_args$copula <- fit$copula

    fit_boot <- tryCatch(do.call(GeoModels::GeoFit, fit_args), error = function(e) NULL)
    if (fit_converged(fit_boot)) {
      estimates <- rbind(estimates, as.numeric(unlist(fit_boot$param[parameter_names])))
    }
  }

  if (is.null(estimates) || nrow(estimates) < 2L) {
    out <- rep(NA_real_, length(parameter_names))
    names(out) <- parameter_names
    return(list(stderr = out, n_success = ifelse(is.null(estimates), 0L, nrow(estimates))))
  }
  colnames(estimates) <- parameter_names
  list(stderr = apply(estimates, 2, sd), n_success = nrow(estimates))
}

coverage_replication <- function(scenario, rep_id, neighb, boot, seed_base, verbose = TRUE) {
  progress(
    sprintf(
      "%s coverage rep %04d | simulating and fitting Clayton (m=%d, K=%d)",
      scenario$scenario_id, as.integer(rep_id), as.integer(neighb), as.integer(boot)
    ),
    verbose = verbose
  )
  design <- create_design(scenario, seed_base)
  data <- simulate_clayton_dataset(scenario, design$coords, design$X, rep_id, seed_base)
  true_param <- true_params(scenario)
  start_time <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    fit_model_object(data, design$coords, design$X, scenario, "clayton", neighb, sensitivity = TRUE),
    error = function(e) NULL
  )
  fit_elapsed <- proc.time()[["elapsed"]] - start_time
  boot_seed <- seed_base + 5000000L + scenario_index(scenario$scenario_id) * 100000L +
    as.integer(neighb) * 10000L + as.integer(rep_id)
  start_boot <- proc.time()[["elapsed"]]
  boot_res <- bootstrap_stderr(
    fit, boot, boot_seed,
    verbose = verbose,
    label = sprintf("%s coverage rep %04d | ", scenario$scenario_id, as.integer(rep_id))
  )
  boot_elapsed <- proc.time()[["elapsed"]] - start_boot

  rows <- lapply(parameter_names, function(param) {
    estimate <- param_or_na(fit, param)
    se <- as.numeric(boot_res$stderr[[param]])
    half_width <- qnorm(0.975) * se
    lower <- estimate - half_width
    upper <- estimate + half_width
    data.frame(
      scenario_id = scenario$scenario_id,
      rep = as.integer(rep_id),
      n = as.integer(scenario$n),
      alpha = as.numeric(scenario$alpha),
      kappa = as.numeric(scenario$kappa),
      nu = as.numeric(scenario$nu),
      fitted_model = "clayton",
      neighb = as.integer(neighb),
      parameter = param,
      estimate = estimate,
      true_value = true_param[[param]],
      bootstrap_se = se,
      ci_lower = lower,
      ci_upper = upper,
      covered = is.finite(lower) && is.finite(upper) && true_param[[param]] >= lower && true_param[[param]] <= upper,
      fit_converged = fit_converged(fit),
      bootstrap_success = boot_res$n_success,
      fit_elapsed_seconds = fit_elapsed,
      bootstrap_elapsed_seconds = boot_elapsed,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

run_coverage_chunk <- function(task) {
  ensure_output_dirs()
  outfile <- coverage_chunk_file(task)
  if (file.exists(outfile) && !isTRUE(task$overwrite)) {
    progress(
      sprintf("Skipping existing coverage chunk %s", basename(outfile)),
      verbose = task$verbose
    )
    return(read.csv(outfile, stringsAsFactors = FALSE))
  }
  scenario <- scenario_by_id(task$scenario_id)
  progress(
    sprintf(
      "Starting coverage chunk %s | reps %04d-%04d | m=%d | K=%d",
      basename(outfile), task$rep_start, task$rep_end, as.integer(task$neighb), as.integer(task$boot)
    ),
    verbose = task$verbose
  )
  rows <- lapply(seq.int(task$rep_start, task$rep_end), function(rep_id) {
    coverage_replication(scenario, rep_id, task$neighb, task$boot, task$seed_base, task$verbose)
  })
  ans <- do.call(rbind, rows)
  write.csv(ans, outfile, row.names = FALSE)
  progress(
    sprintf("Finished coverage chunk %s | rows=%d", basename(outfile), nrow(ans)),
    verbose = task$verbose
  )
  ans
}

run_chunks_parallel <- function(tasks, workers, kind) {
  workers <- resolve_workers(workers, nrow(tasks))
  message("Running ", nrow(tasks), " chunks with ", workers, " worker(s).")
  task_list <- split(tasks, seq_len(nrow(tasks)))
  root <- project_root()

  if (workers <= 1L) {
    return(do.call(rbind, lapply(task_list, function(task) {
      if (identical(kind, "simulation")) run_simulation_chunk(task) else run_coverage_chunk(task)
    })))
  }

  cl <- parallel::makeCluster(workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  ans <- tryCatch(parallel::parLapplyLB(cl, task_list, function(task, root, kind) {
    setwd(root)
    source("src/simulation_config.R")
    source("src/simulation_utils.R")
    suppressPackageStartupMessages(library(GeoModels))
    if (identical(kind, "simulation")) run_simulation_chunk(task) else run_coverage_chunk(task)
  }, root = root, kind = kind), error = function(e) {
    message("Parallel execution failed: ", conditionMessage(e))
    message("Falling back to sequential execution.")
    NULL
  })
  if (is.null(ans)) {
    return(do.call(rbind, lapply(task_list, function(task) {
      if (identical(kind, "simulation")) run_simulation_chunk(task) else run_coverage_chunk(task)
    })))
  }
  do.call(rbind, ans)
}

read_chunk_dir <- function(path) {
  files <- list.files(path, pattern = "\\.csv$", full.names = TRUE)
  if (length(files) == 0L) return(data.frame())
  do.call(rbind, lapply(files, read.csv, stringsAsFactors = FALSE))
}

drop_simulation_se_columns <- function(raw) {
  se_cols <- c("se_mean", "se_mean1", "se_scale", "se_shape")
  raw[, setdiff(names(raw), se_cols), drop = FALSE]
}

dedupe_simulation_raw <- function(raw) {
  if (nrow(raw) == 0L) return(raw)
  raw <- drop_simulation_se_columns(raw)
  key <- paste(raw$scenario_id, raw$rep, raw$fitted_model, sep = "||")
  raw[!duplicated(key, fromLast = TRUE), , drop = FALSE]
}

dedupe_coverage_raw <- function(raw) {
  if (nrow(raw) == 0L) return(raw)
  key <- paste(raw$scenario_id, raw$neighb, raw$rep, raw$parameter, sep = "||")
  raw[!duplicated(key, fromLast = TRUE), , drop = FALSE]
}

clear_simulation_chunks <- function(scenario_ids) {
  ensure_output_dirs()
  for (sid in scenario_ids) {
    files <- list.files(
      "data/processed/simulation_chunks",
      pattern = paste0("^", sid, "_reps_.*\\.csv$"),
      full.names = TRUE
    )
    if (length(files) > 0L) file.remove(files)
  }
}

clear_coverage_chunks <- function(scenario_id, m_values) {
  ensure_output_dirs()
  for (m in m_values) {
    files <- list.files(
      "data/processed/coverage_chunks",
      pattern = paste0("^", scenario_id, "_m", m, "_reps_.*\\.csv$"),
      full.names = TRUE
    )
    if (length(files) > 0L) file.remove(files)
  }
}

parameter_recovery_summary <- function(raw) {
  raw <- dedupe_simulation_raw(raw)
  rows <- list()
  idx <- 1L
  for (sid in unique(raw$scenario_id)) {
    for (model_id in unique(raw$fitted_model)) {
      part <- raw[raw$scenario_id == sid & raw$fitted_model == model_id, , drop = FALSE]
      for (param in parameter_names) {
        est <- part[[paste0("estimate_", param)]]
        truth <- part[[paste0("true_", param)]]
        valid <- part$converged & is.finite(est) & is.finite(truth)
        err <- est[valid] - truth[valid]
        rows[[idx]] <- data.frame(
          scenario_id = sid,
          scenario_group = part$scenario_group[1],
          fitted_model = model_id,
          fitted_model_label = part$fitted_model_label[1],
          parameter = param,
          true_value = ifelse(any(valid), truth[valid][1], NA_real_),
          n_reps = length(unique(part$rep)),
          n_success = sum(valid),
          success_rate = sum(valid) / length(unique(part$rep)),
          mean_estimate = ifelse(any(valid), mean(est[valid]), NA_real_),
          bias = ifelse(any(valid), mean(err), NA_real_),
          mse = ifelse(any(valid), mean(err^2), NA_real_),
          empirical_sd = ifelse(sum(valid) > 1L, sd(est[valid]), NA_real_),
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
      }
    }
  }
  do.call(rbind, rows)
}

model_selection_summary <- function(raw) {
  raw <- dedupe_simulation_raw(raw)
  rows <- list()
  idx <- 1L
  for (sid in unique(raw$scenario_id)) {
    part <- raw[raw$scenario_id == sid, , drop = FALSE]
    reps <- unique(part$rep)
    models <- model_specs$model_id
    selected_claic <- setNames(rep(0L, length(models)), models)
    selected_clbic <- setNames(rep(0L, length(models)), models)
    selected_loglik <- setNames(rep(0L, length(models)), models)
    n_valid_claic <- 0L
    n_valid_clbic <- 0L
    n_valid_loglik <- 0L

    for (rep_id in reps) {
      rr <- part[part$rep == rep_id & part$converged, , drop = FALSE]
      if (nrow(rr) > 0L && any(is.finite(rr$CLAIC_approx))) {
        best <- rr$fitted_model[which.min(rr$CLAIC_approx)]
        selected_claic[[best]] <- selected_claic[[best]] + 1L
        n_valid_claic <- n_valid_claic + 1L
      }
      if (nrow(rr) > 0L && any(is.finite(rr$CLBIC_approx))) {
        best <- rr$fitted_model[which.min(rr$CLBIC_approx)]
        selected_clbic[[best]] <- selected_clbic[[best]] + 1L
        n_valid_clbic <- n_valid_clbic + 1L
      }
      if (nrow(rr) > 0L && any(is.finite(rr$logCompLik))) {
        best <- rr$fitted_model[which.max(rr$logCompLik)]
        selected_loglik[[best]] <- selected_loglik[[best]] + 1L
        n_valid_loglik <- n_valid_loglik + 1L
      }
    }

    clayton <- part[part$fitted_model == "clayton", c("rep", "CLAIC_approx", "CLBIC_approx", "logCompLik")]
    names(clayton) <- c("rep", "clayton_claic", "clayton_clbic", "clayton_loglik")
    for (model_id in models) {
      mm <- part[part$fitted_model == model_id, , drop = FALSE]
      merged <- merge(mm, clayton, by = "rep", all.x = TRUE)
      rows[[idx]] <- data.frame(
        scenario_id = sid,
        scenario_group = part$scenario_group[1],
        fitted_model = model_id,
        fitted_model_label = model_specs$model_label[match(model_id, model_specs$model_id)],
        n_reps = length(reps),
        selected_CLAIC_n = selected_claic[[model_id]],
        selected_CLAIC_pct = ifelse(n_valid_claic > 0L, selected_claic[[model_id]] / n_valid_claic, NA_real_),
        selected_CLBIC_n = selected_clbic[[model_id]],
        selected_CLBIC_pct = ifelse(n_valid_clbic > 0L, selected_clbic[[model_id]] / n_valid_clbic, NA_real_),
        selected_logLik_n = selected_loglik[[model_id]],
        selected_logLik_pct = ifelse(n_valid_loglik > 0L, selected_loglik[[model_id]] / n_valid_loglik, NA_real_),
        mean_delta_CLAIC_vs_Clayton = mean(merged$CLAIC_approx - merged$clayton_claic, na.rm = TRUE),
        mean_delta_CLBIC_vs_Clayton = mean(merged$CLBIC_approx - merged$clayton_clbic, na.rm = TRUE),
        mean_delta_logLik_vs_Clayton = mean(merged$logCompLik - merged$clayton_loglik, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  do.call(rbind, rows)
}

coverage_summary <- function(raw) {
  raw <- dedupe_coverage_raw(raw)
  rows <- list()
  idx <- 1L
  for (m in sort(unique(raw$neighb))) {
    for (param in unique(raw$parameter)) {
      part <- raw[raw$neighb == m & raw$parameter == param, , drop = FALSE]
      valid <- part$fit_converged & is.finite(part$estimate) & is.finite(part$bootstrap_se)
      err <- part$estimate[valid] - part$true_value[valid]
      empirical_sd <- ifelse(sum(valid) > 1L, sd(part$estimate[valid]), NA_real_)
      avg_se <- ifelse(any(valid), mean(part$bootstrap_se[valid]), NA_real_)
      rows[[idx]] <- data.frame(
        scenario_id = part$scenario_id[1],
        neighb = m,
        parameter = param,
        n_reps = length(unique(part$rep)),
        n_success = sum(valid),
        success_rate = sum(valid) / length(unique(part$rep)),
        true_value = ifelse(any(valid), part$true_value[valid][1], NA_real_),
        mean_estimate = ifelse(any(valid), mean(part$estimate[valid]), NA_real_),
        bias = ifelse(any(valid), mean(err), NA_real_),
        mse = ifelse(any(valid), mean(err^2), NA_real_),
        empirical_sd = empirical_sd,
        average_bootstrap_se = avg_se,
        se_sd_ratio = avg_se / empirical_sd,
        coverage_95 = ifelse(any(valid), mean(part$covered[valid]), NA_real_),
        mean_bootstrap_success = ifelse(any(valid), mean(part$bootstrap_success[valid]), NA_real_),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  do.call(rbind, rows)
}

plot_selection_frequencies <- function(selection, path) {
  pdf(path, width = 9, height = 5.5)
  on.exit(dev.off(), add = TRUE)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  scenarios <- unique(selection$scenario_id)
  models <- model_specs$model_id
  mat <- sapply(scenarios, function(sid) {
    part <- selection[selection$scenario_id == sid, , drop = FALSE]
    out <- setNames(rep(0, length(models)), models)
    out[part$fitted_model] <- part$selected_CLAIC_pct
    out
  })
  colnames(mat) <- scenarios
  barplot(
    mat,
    beside = TRUE,
    ylim = c(0, 1),
    col = c("#6B9B6B", "#5A7FA8", "#C85A5A"),
    ylab = "Selection frequency by CLAIC approximation",
    xlab = "Scenario",
    legend.text = model_specs$model_label,
    args.legend = list(x = "topright", bty = "n", cex = 0.8)
  )
}

plot_parameter_recovery <- function(recovery, path) {
  pdf(path, width = 10, height = 6)
  on.exit(dev.off(), add = TRUE)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  par(mfrow = c(1, 2), mar = c(5, 4, 2, 1))
  for (param in c("scale", "shape")) {
    part <- recovery[recovery$parameter == param, , drop = FALSE]
    scenarios <- unique(part$scenario_id)
    models <- model_specs$model_id
    mat <- sapply(scenarios, function(sid) {
      pp <- part[part$scenario_id == sid, , drop = FALSE]
      out <- setNames(rep(NA_real_, length(models)), models)
      out[pp$fitted_model] <- pp$mse
      out
    })
    colnames(mat) <- scenarios
    barplot(
      mat,
      beside = TRUE,
      col = c("#6B9B6B", "#5A7FA8", "#C85A5A"),
      ylab = "MSE",
      xlab = "Scenario",
      main = param,
      legend.text = if (identical(param, "scale")) model_specs$model_label else NULL,
      args.legend = list(x = "topright", bty = "n", cex = 0.75)
    )
  }
}

plot_coverage_calibration <- function(coverage, path) {
  pdf(path, width = 9, height = 5.5)
  on.exit(dev.off(), add = TRUE)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)
  coverage$key <- paste0("m=", coverage$neighb, " ", coverage$parameter)
  vals <- coverage$coverage_95
  names(vals) <- coverage$key
  barplot(vals, las = 2, ylim = c(0, 1), col = "#5A7FA8", ylab = "Empirical 95% coverage")
  abline(h = 0.95, col = "#C85A5A", lty = 2, lwd = 1.5)
}
