# Generate realizations: Clayton-Weibull random field
# Nine realizations with nu ∈ {1,2,4} and smooth ∈ {0,1,2}

rm(list = ls())
library(GeoModels)
library(fields)

# Model configuration
model <- "Weibull"
copula <- "Clayton"
NuisParam("Weibull", num_betas = 1, copula = "Clayton")
CorrParam("GenWend_Matern")

# Spatial coordinates
set.seed(22601)
NN <- 10000
x <- runif(NN)
y <- runif(NN)
coords <- cbind(x, y)

# Fixed parameters
mean <- 0.2
shape <- 5
corrmodel <- "GenWend_Matern"
scale <- 0.15
power2 <- 1/4
nugget <- 0
sill <- 1

# Parameter combinations: smooth ∈ {0,1,2}, nu ∈ {1,2,4}
# Order: rows by smooth, columns by nu (matching caption)
smooth_values <- c(0, 1, 2)
nu_values <- c(1, 2, 4)

# Output directory
output_dir <- "results/figures"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Generate all 9 realizations
for (smooth in smooth_values) {
  for (nu in nu_values) {
    # Simulation parameters
    param <- list(
      smooth = smooth,
      power2 = power2,
      mean = mean,
      nu = nu,
      scale = scale,
      nugget = nugget,
      shape = shape,
      sill = sill
    )
    
    # Simulate Clayton-Weibull random field
    data <- GeoSimCopula(
      coordx = coords,
      corrmodel = corrmodel,
      model = model,
      param = param,
      copula = copula,
      sparse = TRUE
    )$data
    
    # Generate filename
    filename <- sprintf("nu%d_smooth%d.pdf", nu, smooth)
    filepath <- file.path(output_dir, filename)
    
    # Create plot
    pdf(filepath, width = 6, height = 6)
    quilt.plot(coords, data)
    dev.off()
    
    cat(sprintf("Generated: %s\n", filename))
  }
}

cat("\nAll 9 realizations generated successfully.\n")
