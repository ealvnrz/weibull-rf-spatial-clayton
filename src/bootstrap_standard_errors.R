# Bootstrap standard errors for parameter estimates

rm(list = ls())
library(GeoModels)

GeoVarestbootstrap_custom <- function(fit, K = 100, lower = NULL, upper = NULL, 
                                      optimizer = NULL, parallel = FALSE, 
                                      seed = NULL, method = "cholesky", 
                                      alpha = 0.95, sparse = FALSE, ...) {
  
  if (!is.null(seed)) set.seed(seed)
  
  if (length(fit$coordt) == 1) fit$coordt <- NULL
  
  if (is.null(fit$sensmat)) {
    warning("Sensitivity matrix is missing. Use sensitivity=TRUE in GeoFit for CLAIC/CLBIC calculation.")
  }
  
  if (is.null(optimizer)) {
    optimizer <- fit$optimizer
    if (is.null(lower)) lower <- fit$lower
    if (is.null(upper)) upper <- fit$upper
  }
  
  coords <- cbind(fit$coordx, fit$coordy)
  if (fit$bivariate && is.null(fit$coordx_dyn)) {
    coords <- coords[1:(length(fit$coordx)/2), ]
  }
  
  model <- fit$model
  if (fit$missp) {
    if (fit$model == "StudentT") model <- "Gaussian_misp_StudentT"
    if (fit$model == "Poisson") model <- "Gaussian_misp_Poisson"
    if (fit$model == "PoissonZIP") model <- "Gaussian_misp_PoissonZIP"
    if (fit$model == "SkewStudentT") model <- "Gaussian_misp_SkewStudenT"
    if (fit$model == "Tukeygh") model <- "Gaussian_misp_Tukeygh"
  }
  
  dimat <- fit$numtime * fit$numcoord
  tempX <- fit$X
  if (sum(fit$X[1:dimat] == 1) == dimat && !dim(fit$X)[2] > 1) {
    fit$X <- NULL
  }
  
  param_sim <- append(fit$param, fit$fixed)
  
  cat("Parametric bootstrap can be time consuming ...\n")
  
  if (is.null(fit$copula)) {
    cat("Performing", K, "simulations (non-copula model)....\n")
    if (method == "cholesky") {
      data_sim <- GeoSim(
        coordx = coords,
        coordt = fit$coordt,
        coordx_dyn = fit$coordx_dyn,
        anisopars = fit$anisopars,
        corrmodel = fit$corrmodel,
        model = fit$model,
        param = param_sim,
        grid = fit$grid,
        X = fit$X,
        n = fit$n,
        method = method,
        distance = fit$distance,
        radius = fit$radius,
        nrep = K,
        progress = FALSE
      )
    } else {
      stop("Only method='cholesky' is currently supported for non-copula models")
    }
  } else {
    cat("Performing", K, "simulations (copula model)....\n")
    if (method == "cholesky") {
      data_sim <- GeoSimCopula(
        coordx = coords,
        coordt = fit$coordt,
        coordx_dyn = fit$coordx_dyn,
        anisopars = fit$anisopars,
        corrmodel = fit$corrmodel,
        model = fit$model,
        copula = fit$copula,
        param = param_sim,
        grid = fit$grid,
        X = fit$X,
        n = fit$n,
        method = method,
        distance = fit$distance,
        radius = fit$radius,
        nrep = K
      )
    } else {
      stop("Only method='cholesky' is currently supported for copula models")
    }
  }
  
  k <- 1
  res <- NULL
  
  if (!parallel) {
    cat("Performing", K, "estimations...\n")
    
    while (k <= K) {
      start_perturbed <- fit$param
      for (param_name in names(start_perturbed)) {
        param_value <- as.numeric(start_perturbed[[param_name]])
        if (abs(param_value) < 1) {
          noise <- rnorm(1, 0, abs(param_value) * 0.01)
        } else {
          noise <- rnorm(1, 0, abs(param_value) * 0.001)
        }
        start_perturbed[[param_name]] <- param_value + noise
      }
      
      fit_args <- list(
        data = data_sim$data[[k]],
        start = start_perturbed,
        fixed = fit$fixed,
        coordx = coords,
        coordx_dyn = fit$coordx_dyn,
        copula = fit$copula,
        anisopars = fit$anisopars,
        est.aniso = fit$est.aniso,
        lower = lower,
        upper = upper,
        neighb = fit$neighb,
        corrmodel = fit$corrmodel,
        model = model,
        sparse = FALSE,
        n = fit$n,
        maxdist = fit$maxdist,
        maxtime = fit$maxtime,
        optimizer = optimizer,
        grid = fit$grid,
        likelihood = fit$likelihood,
        type = fit$type,
        X = fit$X,
        distance = fit$distance,
        radius = fit$radius
      )
      
      if (!is.null(fit$coordt)) {
        fit_args$coordt <- fit$coordt
      }
      
      fit_args <- c(fit_args, list(...))
      
      res_est <- tryCatch({
        do.call(GeoFit, fit_args)
      }, error = function(e) {
        if (k <= 3) {
          cat("  GeoFit error at iteration ", k, ": ", e$message, "\n")
        }
        return(NULL)
      })
      
      is_valid <- FALSE
      if (is.null(res_est)) {
        if (k <= 3) cat("  Iteration ", k, " failed: res_est is NULL\n")
      } else if (is.null(res_est$param)) {
        if (k <= 3) cat("  Iteration ", k, " failed: res_est$param is NULL\n")
      } else {
        param_vec <- unlist(res_est$param)
        if (any(is.na(param_vec))) {
          if (k <= 3) cat("  Iteration ", k, " failed: parameters contain NA values\n")
        } else if (any(is.infinite(param_vec))) {
          if (k <= 3) cat("  Iteration ", k, " failed: parameters contain infinite values\n")
        } else if (!is.null(res_est$logCompLik) && res_est$logCompLik >= 1.0e8) {
          if (k <= 3) cat("  Iteration ", k, " failed: logCompLik too large (", res_est$logCompLik, ")\n")
        } else {
          is_valid <- TRUE
        }
      }
      
      if (is_valid) {
        res <- rbind(res, unlist(res_est$param))
      }
      
      if (k %% 10 == 0 || k == 1) {
        cat("  Completed", k, "of", K, "estimations\n")
      }
      
      k <- k + 1
    }
  } else {
    stop("Parallel processing not yet implemented in custom function")
  }
  
  if (is.null(res) || nrow(res) == 0) {
    stop("All bootstrap iterations failed. No successful parameter estimates obtained.")
  }
  
  n_success <- nrow(res)
  if (n_success < K * 0.5) {
    warning("Only ", n_success, " out of ", K, " bootstrap iterations succeeded. Results may be unreliable.")
  }
  
  numparam <- length(fit$param)
  colnames(res) <- names(fit$param)
  
  invG <- var(res)
  
  G <- try(solve(invG), silent = TRUE)
  if (!is.matrix(G)) {
    warning("Bootstrap estimated Godambe matrix is singular. Standard errors are still computed from variance-covariance matrix.")
  }
  
  stderr <- sqrt(diag(invG))
  
  if (any(is.na(stderr)) || any(is.infinite(stderr))) {
    warning("Some standard errors are NA or infinite. This may indicate numerical issues.")
    stderr[is.na(stderr) | is.infinite(stderr)] <- NA
  }
  
  if (!is.null(fit$sensmat)) {
    if ((fit$likelihood == "Marginal" && (fit$type == "Independence" || fit$type == "Pairwise")) ||
        (fit$likelihood == "Conditional" && fit$type == "Pairwise")) {
      H <- fit$sensmat
      penalty <- sum(diag(H %*% invG))
      claic <- -2 * fit$logCompLik + 2 * penalty
      clbic <- -2 * fit$logCompLik + log(dimat) * penalty
      fit$varimat <- H %*% invG %*% H
    } else if (fit$likelihood == "Full" && fit$type == "Standard") {
      claic <- -2 * fit$logCompLik + 2 * numparam
      clbic <- -2 * fit$logCompLik + log(dimat) * 2 * numparam
    } else {
      claic <- NULL
      clbic <- NULL
    }
    fit$claic <- claic
    fit$clbic <- clbic
  }
  
  aa <- qnorm(1 - (1 - alpha) / 2) * stderr
  pp <- as.numeric(fit$param)
  low <- pp - aa
  upp <- pp + aa
  fit$conf.int <- rbind(low, upp)
  fit$pvalues <- 2 * pnorm(-abs(pp / stderr))
  
  fit$stderr <- stderr
  fit$varcov <- invG
  fit$estimates <- res
  fit$X <- tempX
  
  cat("Bootstrap completed:", nrow(res), "successful iterations out of", K, "\n")
  
  return(fit)
}

output_dir <- "data/processed"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("Loading fitted models...\n")
all_models <- readRDS(file.path(output_dir, "all_fitted_models_section148.rds"))

fit_direct <- all_models$fit_direct
fit_gaussian <- all_models$fit_gaussian
clayton_fits <- all_models$clayton_fits
nu_values <- all_models$nu_values
coords <- all_models$coords
intensity <- all_models$intensity
X <- all_models$X

nboot <- 100
neighb <- 20
seed <- 12345

cat("Computing bootstrap standard errors...\n")
cat("This may take a while...\n")
cat("NOTE: Only Gaussian copula bootstrap will be executed (other models already computed).\n\n")

se_table_file <- file.path(output_dir, "bootstrap_standard_errors_section148.csv")
if (!file.exists(se_table_file)) {
  se_table <- data.frame(
    Model = character(),
    mean_se = numeric(),
    mean1_se = numeric(),
    shape_se = numeric(),
    scale_se = numeric(),
    stringsAsFactors = FALSE
  )
  write.csv(se_table, se_table_file, row.names = FALSE)
  cat("Created new CSV file for standard errors.\n")
} else {
  cat("Using existing CSV file (will update/replace Gaussian copula row).\n")
}

add_se_to_table <- function(model_name, se_values) {
  if (is.null(se_values)) {
    cat("  Warning: No standard errors available for", model_name, "\n")
    return()
  }
  
  new_row <- data.frame(
    Model = model_name,
    mean_se = ifelse("mean" %in% names(se_values), as.numeric(se_values["mean"]), NA),
    mean1_se = ifelse("mean1" %in% names(se_values), as.numeric(se_values["mean1"]), NA),
    shape_se = ifelse("shape" %in% names(se_values), as.numeric(se_values["shape"]), NA),
    scale_se = ifelse("scale" %in% names(se_values), as.numeric(se_values["scale"]), NA),
    stringsAsFactors = FALSE
  )
  
  if (file.exists(se_table_file)) {
    se_table <- read.csv(se_table_file, stringsAsFactors = FALSE)
    model_idx <- which(se_table$Model == model_name)
    if (length(model_idx) > 0) {
      se_table[model_idx[1], ] <- new_row
      cat("  Replaced existing row for", model_name, "in CSV.\n")
    } else {
      se_table <- rbind(se_table, new_row)
      cat("  Added new row for", model_name, "to CSV.\n")
    }
  } else {
    se_table <- new_row
    cat("  Created new CSV file with", model_name, ".\n")
  }
  
  write.csv(se_table, se_table_file, row.names = FALSE)
  cat("  Standard errors saved to CSV.\n")
}

cat("Bootstrap for Gaussian copula model...\n")
bootstrap_gaussian <- GeoVarestbootstrap_custom(
  fit_gaussian,
  K = nboot,
  lower = list(mean = -100, mean1 = -100, shape = 0.00001, scale = 0),
  upper = list(mean = 100, mean1 = 100, shape = 100, scale = 100),
  optimizer = "nlminb",
  parallel = FALSE,
  seed = seed
)
se_gaussian <- bootstrap_gaussian$stderr
add_se_to_table("Gaussian copula", se_gaussian)
cat("  Completed.\n\n")

se_table_file <- file.path(output_dir, "bootstrap_standard_errors_section148.csv")
if (file.exists(se_table_file)) {
  se_table <- read.csv(se_table_file, stringsAsFactors = FALSE)
  cat("\nFinal standard errors table:\n")
  print(se_table)
} else {
  warning("CSV file not found. This should not happen if bootstrap completed successfully.")
  se_table <- NULL
}

bootstrap_results_file <- file.path(output_dir, "bootstrap_results_section148.rds")
if (file.exists(bootstrap_results_file)) {
  existing_results <- readRDS(bootstrap_results_file)
  existing_results$bootstrap_gaussian <- bootstrap_gaussian
  existing_results$se_gaussian <- se_gaussian
  if (!is.null(se_table)) {
    existing_results$se_table <- se_table
  }
  saveRDS(existing_results, bootstrap_results_file)
} else {
  saveRDS(list(
    bootstrap_gaussian = bootstrap_gaussian,
    se_gaussian = se_gaussian,
    se_table = se_table
  ), bootstrap_results_file)
}

cat("Bootstrap standard errors computed and saved.\n")
cat("Results saved to:", file.path(output_dir, "bootstrap_results_section148.rds"), "\n")
cat("Standard errors table saved to:", file.path(output_dir, "bootstrap_standard_errors_section148.csv"), "\n")
