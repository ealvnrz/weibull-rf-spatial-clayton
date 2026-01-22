# Fit spatial models: Weibull random fields with different constructions
# Comparing χ² transformation, Gaussian copula, and Clayton copula

rm(list = ls())
library(GeoModels)

data_file <- "data/raw/LAS_sec_148.csv"
data <- read.csv(data_file)
coords <- cbind(data$x, data$y)
intensity <- data$intensity_scaled

output_dir <- "data/processed"
figs_dir <- "results/figures"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
if (!dir.exists(figs_dir)) dir.create(figs_dir, recursive = TRUE)

# Model configuration
model <- "Weibull"
corrmodel <- "GenWend_Matern"
power2 <- 1/4
nugget <- 0
smooth <- 0
sill <- 1
neighb <- 20

X <- cbind(rep(1, nrow(data)), data$b)

# Starting values
mean <- 0
mean1 <- 0.5
shape <- 3
scale <- 10

# Parameter bounds
I <- 100
lower <- list(mean = -I, mean1 = -I, shape = 0.00001, scale = 0)
upper <- list(mean = I, mean1 = I, shape = I, scale = I)

# Independence model
lower_ind <- list(mean = -I, shape = 0.00001)
upper_ind <- list(mean = I, shape = I)
fixed_ind <- list(power2 = power2, nugget = nugget, smooth = smooth, scale = 0.0001)
start_ind <- list(mean = mean, shape = shape)

fit_ind <- GeoFit(
  data = intensity,
  coordx = coords,
  model = model,
  corrmodel = corrmodel,
  start = start_ind,
  fixed = fixed_ind,
  lower = lower_ind,
  upper = upper_ind,
  optimizer = "nlminb",
  likelihood = "Marginal",
  type = "Independence",
  sensitivity = TRUE
)

# Weibull RF via χ² transformation
fixed_direct <- list(power2 = power2, nugget = nugget, smooth = smooth)
start_direct <- list(mean = mean, mean1 = mean1, shape = shape, scale = scale)

fit_direct <- GeoFit2(
  data = intensity,
  X = X,
  coordx = coords,
  model = model,
  corrmodel = corrmodel,
  start = start_direct,
  fixed = fixed_direct,
  lower = lower,
  upper = upper,
  optimizer = "nlminb",
  likelihood = "Marginal",
  type = "Pairwise",
  sensitivity = TRUE,
  neighb = neighb
)

# Weibull RF via Gaussian copula
copula_gaussian <- "Gaussian"
fixed_gaussian <- list(power2 = power2, nugget = nugget, sill = sill, smooth = smooth)
start_gaussian <- list(mean = mean, mean1 = mean1, shape = shape, scale = scale)

fit_gaussian <- GeoFit2(
  data = intensity,
  X = X,
  coordx = coords,
  model = model,
  corrmodel = corrmodel,
  copula = copula_gaussian,
  start = start_gaussian,
  fixed = fixed_gaussian,
  lower = lower,
  upper = upper,
  optimizer = "nlminb",
  likelihood = "Marginal",
  type = "Pairwise",
  sensitivity = TRUE,
  neighb = neighb
)

# Weibull RF via Clayton copula (multiple nu values)
copula_clayton <- "Clayton"
nu_values <- c(1, 2, 4, 6)
clayton_fits <- list()

for (i in seq_along(nu_values)) {
  nu <- nu_values[i]
  fixed_clayton <- list(power2 = power2, nugget = nugget, nu = nu, sill = sill, smooth = smooth)
  start_scale <- ifelse(nu >= 4, 15, scale)
  start_clayton <- list(mean = mean, mean1 = mean1, shape = shape, scale = start_scale)
  
  fit_clayton <- GeoFit2(
    data = intensity,
    X = X,
    coordx = coords,
    model = model,
    corrmodel = corrmodel,
    copula = copula_clayton,
    start = start_clayton,
    fixed = fixed_clayton,
    lower = lower,
    upper = upper,
    optimizer = "nlminb",
    likelihood = "Marginal",
    type = "Pairwise",
    sensitivity = TRUE,
    neighb = neighb
  )
  
  clayton_fits[[i]] <- fit_clayton
}

# Model comparison
n <- length(intensity)
model_comparison <- data.frame(
  Model = c("Independence", "χ² transformation", 
            paste("Clayton (nu=", nu_values, ")", sep = ""),
            "Gaussian copula"),
  LogLik = c(fit_ind$logCompLik, fit_direct$logCompLik,
             sapply(clayton_fits, function(x) x$logCompLik),
             fit_gaussian$logCompLik),
  AIC = c(2*2 - 2*fit_ind$logCompLik, 2*4 - 2*fit_direct$logCompLik,
          sapply(clayton_fits, function(x) 2*4 - 2*x$logCompLik),
          2*4 - 2*fit_gaussian$logCompLik),
  BIC = c(log(n)*2 - 2*fit_ind$logCompLik, log(n)*4 - 2*fit_direct$logCompLik,
          sapply(clayton_fits, function(x) log(n)*4 - 2*x$logCompLik),
          log(n)*4 - 2*fit_gaussian$logCompLik),
  mean = c(as.numeric(fit_ind$param["mean"]), as.numeric(fit_direct$param["mean"]),
           sapply(clayton_fits, function(x) as.numeric(x$param["mean"])),
           as.numeric(fit_gaussian$param["mean"])),
  mean1 = c(NA, as.numeric(fit_direct$param["mean1"]),
            sapply(clayton_fits, function(x) as.numeric(x$param["mean1"])),
            as.numeric(fit_gaussian$param["mean1"])),
  shape = c(as.numeric(fit_ind$param["shape"]), as.numeric(fit_direct$param["shape"]),
            sapply(clayton_fits, function(x) as.numeric(x$param["shape"])),
            as.numeric(fit_gaussian$param["shape"])),
  scale = c(NA, as.numeric(fit_direct$param["scale"]),
            sapply(clayton_fits, function(x) as.numeric(x$param["scale"])),
            as.numeric(fit_gaussian$param["scale"])),
  Convergence = c(fit_ind$convergence == 0 || fit_ind$convergence == "Successful",
                  fit_direct$convergence == 0 || fit_direct$convergence == "Successful",
                  sapply(clayton_fits, function(x) x$convergence == 0 || x$convergence == "Successful"),
                  fit_gaussian$convergence == 0 || fit_gaussian$convergence == "Successful")
)

write.csv(model_comparison, 
          file.path(output_dir, "model_comparison_section148.csv"), 
          row.names = FALSE)

parameter_estimates <- model_comparison[, c("Model", "mean", "mean1", "shape", "scale", 
                                            "LogLik", "AIC", "BIC", "Convergence")]
write.csv(parameter_estimates,
          file.path(output_dir, "parameter_estimates_section148.csv"),
          row.names = FALSE)

# Save fitted models for bootstrap analysis
saveRDS(fit_ind, file.path(output_dir, "fit_independence_section148.rds"))
saveRDS(fit_direct, file.path(output_dir, "fit_direct_section148.rds"))
saveRDS(fit_gaussian, file.path(output_dir, "fit_gaussian_section148.rds"))
for (i in seq_along(clayton_fits)) {
  saveRDS(clayton_fits[[i]], file.path(output_dir, paste0("fit_clayton_nu", nu_values[i], "_section148.rds")))
}
saveRDS(list(fit_ind = fit_ind, fit_direct = fit_direct, fit_gaussian = fit_gaussian, 
             clayton_fits = clayton_fits, nu_values = nu_values, coords = coords, 
             intensity = intensity, X = X),
        file.path(output_dir, "all_fitted_models_section148.rds"))

# Calculate variograms for plotting
res_direct <- GeoResiduals(fit_direct)
vario_w <- GeoVariogram(data = res_direct$data, coordx = coords, maxdist = 8)

res_gaussian <- GeoResiduals(fit_gaussian)
vario_g <- GeoVariogram(data = res_gaussian$data, coordx = coords, maxdist = 8)

clayton_residuals <- list()
clayton_variograms <- list()
for (i in seq_along(clayton_fits)) {
  clayton_residuals[[i]] <- GeoResiduals(clayton_fits[[i]])
  clayton_variograms[[i]] <- GeoVariogram(data = clayton_residuals[[i]]$data, coordx = coords, maxdist = 8)
}

# Select best Clayton model (lowest AIC)
clayton_aic <- sapply(clayton_fits, function(x) 2*4 - 2*x$logCompLik)
best_clayton_idx <- which.min(clayton_aic)
best_clayton_nu <- nu_values[best_clayton_idx]

# Helper function for plotting theoretical variogram
GeoCovariogram2 <- function(fitted, vario, variogram = TRUE, ...) {
  if(!inherits(fitted, "GeoFit")) stop("Enter an object obtained of class GeoFit\n")
  if(!inherits(vario, "GeoVariogram")) stop("Enter an object of class GeoVariogram\n")
  space <- !fitted$bivariate && !fitted$spacetime
  h <- seq(9e-5, 8, 0.01)
  if(space) {
    cc <- GeoCorrFct(x = h, corrmodel = fitted$corrmodel, covariance = TRUE, variogram = variogram,
                     param = append(fitted$param, fitted$fixed), model = fitted$model)
    plot(cc, type = "l", xlab = "Distance", ylab = "Semivariogram", col = "red", xlim = c(0, 8))
    box()
  }
}

# Combined variogram plot (best Clayton model only)
pdf(file.path(figs_dir, "combined_variograms_section148.pdf"), width = 10, height = 6)
par(bg = "transparent")
plot(1, 1, type = "n", xlab = "", ylab = "", xlim = c(0, 8), ylim = c(0, 0.00028))
GeoCovariogram2(res_gaussian, vario = vario_g, pch = 20, ylim = c(0, 0.00025))
GeoCovariogram(res_direct, show.vario = TRUE, vario = vario_w, pch = 20, 
               col = "blue", ylim = c(0, 0.00028), add.vario = TRUE)
GeoCovariogram(clayton_residuals[[best_clayton_idx]], show.vario = TRUE, 
               vario = clayton_variograms[[best_clayton_idx]], 
               pch = 20, col = "black", add.vario = TRUE)
box()
dev.off()

# Empirical variogram (points only)
pdf(file.path(figs_dir, "empirical_variogram_section148.pdf"), width = 6, height = 6)
par(bg = "transparent")
plot(1, 1, type = "n", xlab = "Distance", ylab = "Semivariogram", xlim = c(0, 8), ylim = c(0, 0.00025))
points(vario_g$centers, vario_g$variograms, pch = 20)
box()
dev.off()
