# Generate correlation function: Clayton-Weibull random field
# Correlation function for nu=2 and kappa=2 with smooth ∈ {0,1,2}

rm(list = ls())
library(GeoModels)

# Parameters
kappa <- 2
nu <- 2
scale <- 0.1
power2 <- 1/3.5
h_max <- 0.5
n_points <- 200
smooth_values <- c(0, 1, 2)

# Helper functions
pochhammer <- function(a, k) {
  if (k == 0) return(1)
  exp(lgamma(a + k) - lgamma(a))
}

A_coeff <- function(m, n, nu) {
  nu2 <- nu / 2
  num <- pochhammer(nu2 + 1, m + n)^2
  den <- pochhammer(nu2, m) * pochhammer(1, n) * factorial(m) * factorial(n)
  return(num / den)
}

J_integral <- function(m, n, kappa, nu) {
  nu_kappa <- 1 / gamma(1 + 1/kappa)
  if (nu == 2) {
    integrand <- function(u) {
      result <- numeric(length(u))
      valid <- (u > 0) & (u < 1)
      if (any(valid)) {
        log_term <- -log(1 - u[valid])
        log_term[log_term <= 0] <- 0
        result[valid] <- nu_kappa * (log_term^(1/kappa)) * (u[valid]^m) * ((1 - u[valid])^n)
      }
      return(result)
    }
  } else {
    integrand <- function(x) {
      result <- numeric(length(x))
      valid <- (x > 0) & (x < 1)
      if (any(valid)) {
        log_term <- -log(1 - x[valid]^(nu/2))
        log_term[log_term <= 0] <- 0
        result[valid] <- nu_kappa * (nu/2) * (log_term^(1/kappa)) * 
                         (x[valid]^(m + nu/2 - 1)) * ((1 - x[valid])^n)
      }
      return(result)
    }
  }
  result <- tryCatch({
    integrate(integrand, lower = 0, upper = 1, rel.tol = 1e-8, abs.tol = 1e-10)$value
  }, error = function(e) {
    warning(paste("Integration failed for m=", m, ", n=", n, ", kappa=", kappa, ", nu=", nu, ": ", e$message))
    return(0)
  })
  return(result)
}

compute_cross_moment <- function(rho, kappa, nu, max_terms = 50, J_cache = NULL) {
  if (abs(rho - 1) < 1e-10) {
    nu_kappa <- 1 / gamma(1 + 1/kappa)
    variance <- gamma(1 + 2/kappa) * nu_kappa^2 - 1
    return(variance + 1)
  }
  
  rho2 <- rho^2
  if (rho2 > 0.95) {
    max_terms <- max(max_terms, 150)
  } else if (rho2 > 0.9) {
    max_terms <- max(max_terms, 100)
  } else if (rho2 > 0.7) {
    max_terms <- max(max_terms, 70)
  }
  
  if (rho2 >= 1) rho2 <- 0.99999
  
  prefactor <- (1 - rho2)^(nu/2 + 1)
  
  if (is.null(J_cache)) J_cache <- list()
  
  get_J <- function(m, n) {
    key <- paste(m, n, sep = ",")
    if (key %in% names(J_cache)) {
      return(J_cache[[key]])
    } else {
      return(J_integral(m, n, kappa, nu))
    }
  }
  
  total <- 0
  prev_total <- 0
  converged <- FALSE
  
  for (m in 0:max_terms) {
    row_sum <- 0
    max_term_in_row <- 1e-10
    
    for (n in 0:max_terms) {
      if (m + n > max_terms) next
      A_val <- A_coeff(m, n, nu)
      J_val <- get_J(m, n)
      term <- A_val * (rho2^(m + n)) * (J_val^2)
      if (is.na(term) || is.infinite(term) || term == 0) {
        if (n > 10) break
        next
      }
      row_sum <- row_sum + term
      max_term_in_row <- max(max_term_in_row, abs(term))
      if (n > 10 && max_term_in_row > 0 && abs(term) < 1e-12 * max_term_in_row) break
    }
    
    total <- total + row_sum
    if (m > 10) {
      rel_change <- abs(total - prev_total) / (abs(total) + 1e-10)
      if (rel_change < 1e-10) {
        converged <- TRUE
        break
      }
    }
    prev_total <- total
  }
  
  if (!converged && rho2 > 0.9) {
    warning(sprintf("Series may not have converged for rho=%.4f (max_terms=%d)", rho, max_terms))
  }
  
  return(prefactor * total)
}

clayton_weibull_correlation <- function(rho, kappa, nu, max_terms = 50, J_cache = NULL) {
  if (abs(rho - 1) < 1e-10) return(1.0)
  
  if (rho > 0.95) {
    rho_ref <- 0.95
    nu_kappa <- 1 / gamma(1 + 1/kappa)
    variance <- gamma(1 + 2/kappa) * nu_kappa^2 - 1
    cross_moment_ref <- compute_cross_moment(rho_ref, kappa, nu, max_terms, J_cache)
    corr_ref <- (cross_moment_ref - 1) / variance
    alpha <- (1 - rho) / (1 - rho_ref)
    return(1 - alpha * (1 - corr_ref))
  }
  
  nu_kappa <- 1 / gamma(1 + 1/kappa)
  variance <- gamma(1 + 2/kappa) * nu_kappa^2 - 1
  cross_moment <- compute_cross_moment(rho, kappa, nu, max_terms, J_cache)
  return((cross_moment - 1) / variance)
}

get_underlying_correlation <- function(h, smooth, scale, power2) {
  corrmodel <- "GenWend_Matern"
  param <- list(smooth = smooth, scale = scale, power2 = power2, sill = 1, nugget = 0)
  corr_result <- GeoCorrFct(x = h, corrmodel = corrmodel, param = param, 
                            distance = "Eucl", covariance = FALSE)
  return(corr_result$corr)
}

# Distance vector
h <- seq(0, h_max, length.out = n_points)

# Precompute J integrals
J_cache <- list()
max_J_terms <- 150
for (m in 0:max_J_terms) {
  for (n in 0:max_J_terms) {
    if (m + n <= max_J_terms) {
      key <- paste(m, n, sep = ",")
      J_cache[[key]] <- J_integral(m, n, kappa, nu)
    }
  }
}

# Compute correlations for each smooth parameter
results <- list()
for (smooth in smooth_values) {
  rho_h <- get_underlying_correlation(h, smooth, scale, power2)
  corr_clayton <- numeric(length(h))
  for (i in seq_along(h)) {
    adaptive_max_terms <- if (rho_h[i] > 0.95) 150 else if (rho_h[i] > 0.9) 100 else if (rho_h[i] > 0.7) 70 else 50
    corr_clayton[i] <- clayton_weibull_correlation(rho_h[i], kappa, nu, 
                                                   max_terms = adaptive_max_terms, 
                                                   J_cache = J_cache)
  }
  results[[paste0("smooth_", smooth)]] <- list(h = h, rho = rho_h, correlation = corr_clayton)
}

# Output directory
output_dir <- "results/figures"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Create plot
filename <- "clayton_weibull_correlation_function.pdf"
filepath <- file.path(output_dir, filename)

pdf(filepath, width = 12, height = 4)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 2.5, 1), mgp = c(2.5, 1, 0))

colors <- c("black", "red", "blue")
for (i in seq_along(smooth_values)) {
  smooth <- smooth_values[i]
  key <- paste0("smooth_", smooth)
  
  plot(results[[key]]$h, results[[key]]$rho, 
       type = "l", lwd = 2, col = "gray", lty = 2,
       xlab = "Distance h", ylab = "Correlation", ylim = c(0, 1))
  
  lines(results[[key]]$h, results[[key]]$correlation,
        col = colors[i], lty = 1, lwd = 2)
}

dev.off()

cat(sprintf("Generated: %s\n", filename))
