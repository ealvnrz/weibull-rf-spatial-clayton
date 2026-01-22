# Simulate Weibull random fields with Gaussian copula
# Parameter recovery study with 1000 simulation runs

rm(list = ls())
library(GeoModels)

# Model configuration
model <- "Weibull"
copula <- "Gaussian"

# Spatial coordinates
NN <- 800
x <- runif(NN)
y <- runif(NN)
coords <- cbind(x, y)

# Regression matrix
X <- cbind(rep(1, NN), runif(NN))

# Parameters
mean <- 0.2
mean1 <- -0.2
shape <- 1
corrmodel <- "GenWend_Matern"
scale <- 0.06
power2 <- 1/3.5
smooth <- 0
sill <- 1
nugget <- 0

# Parameter list
param <- list(
  smooth = smooth,
  power2 = power2,
  mean = mean,
  mean1 = mean1,
  sill = sill,
  scale = scale,
  nugget = nugget,
  shape = shape
)

# Simulation settings
set.seed(7689)
n_sim <- 1000
I <- 20

# Storage for results
res <- NULL
k <- 1
n_failed <- 0

while(k <= n_sim) {
  cat("Simulation", k, "of", n_sim)
  if (n_failed > 0) {
    cat(" (", n_failed, "failed)")
  }
  cat("\n")
  
  # Simulate data
  data <- GeoSimCopula(
    coordx = coords,
    corrmodel = corrmodel,
    model = model,
    param = param,
    X = X,
    copula = copula
  )$data
  
  # Estimation setup
  start <- list(mean = mean, mean1 = mean1, scale = scale, shape = shape)
  lower <- list(mean = -I, mean1 = -I, scale = 0, shape = 0)
  upper <- list(mean = I, mean1 = I, scale = I, shape = I)
  fixed <- list(sill = sill, smooth = smooth, power2 = power2, nugget = nugget)
  
  # Fit model
  fit <- GeoFit(
    data = data,
    coordx = coords,
    corrmodel = corrmodel,
    model = model,
    X = X,
    neighb = 2,
    likelihood = "Marginal",
    type = "Pairwise",
    optimizer = "nlminb",
    lower = lower,
    upper = upper,
    copula = copula,
    start = start,
    fixed = fixed
  )
  
  # Store results (only if fit was successful)
  if (!is.null(fit$param) && !any(is.na(fit$param))) {
    res <- rbind(res, fit$param)
    k <- k + 1
  } else {
    n_failed <- n_failed + 1
    cat("  Warning: Optimization failed, retrying...\n")
  }
}

# Process results
res1 <- as.numeric(res[, 1])
res2 <- as.numeric(res[, 2])
res3 <- as.numeric(res[, 3])
res4 <- as.numeric(res[, 4])

valid_sims <- complete.cases(res1, res2, res3, res4)
res1 <- res1[valid_sims]
res2 <- res2[valid_sims]
res3 <- res3[valid_sims]
res4 <- res4[valid_sims]

res <- cbind(res1, res2, res3, res4)
colnames(res) <- c("mean", "mean1", "scale", "shape")

# Calculate summary statistics
true_params <- c(mean, mean1, scale, shape)
param_names <- c("beta0", "beta1", "scale", "shape")
bias <- colMeans(res) - true_params
mse <- apply(res, 2, var) + bias^2

results_summary <- data.frame(
  Parameter = param_names,
  True_Value = true_params,
  Mean_Estimate = colMeans(res),
  Bias = bias,
  MSE = mse,
  SD = apply(res, 2, sd)
)

# Output directory
output_dir <- "data/processed"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Save results as CSV
write.csv(res, file = file.path(output_dir, "weibull_gaussian_copula_simulations.csv"), 
          row.names = FALSE)
write.csv(results_summary, file = file.path(output_dir, "weibull_gaussian_copula_summary.csv"), 
          row.names = FALSE)

cat("\nSimulation complete.\n")
cat("Successful simulations:", nrow(res), "\n")
if (n_failed > 0) {
  cat("Failed simulations:", n_failed, "\n")
}
cat("Results saved to:\n")
cat("- data/processed/weibull_gaussian_copula_simulations.csv\n")
cat("- data/processed/weibull_gaussian_copula_summary.csv\n")
