# Generate boxplots: parameter estimates comparison
# Boxplots comparing parameter estimates across three Weibull constructions

rm(list = ls())

# Input data paths
data_dir <- "data/processed"
clayton_file <- file.path(data_dir, "weibull_clayton_copula_simulations.csv")
gaussian_file <- file.path(data_dir, "weibull_gaussian_copula_simulations.csv")
direct_file <- file.path(data_dir, "weibull_direct_simulations.csv")

# Load simulation results
clayton_results <- read.csv(clayton_file)
gaussian_results <- read.csv(gaussian_file)
direct_results <- read.csv(direct_file)

# Output directory
output_dir <- "results/figures"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Create PDF
filename <- "combined_parameter_boxplots.pdf"
filepath <- file.path(output_dir, filename)
pdf(filepath, width = 12, height = 8)

# Set up plot layout: 2x2 grid
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1), oma = c(0, 0, 0, 0))

# Parameters and labels
parameters <- c("beta0", "beta1", "scale", "shape")
param_labels <- c(expression(beta[0]), expression(beta[1]), 
                 expression(alpha), expression(kappa))
true_values <- c(0.2, -0.2, 0.06, 1.0)

# Colors and model names (muted but distinguishable colors)
colors <- c("#C85A5A", "#5A7FA8", "#6B9B6B")  # Muted red, muted blue, muted green
model_names <- c("Clayton", "Gaussian", expression(chi^2))

# Create boxplot for each parameter
for (i in 1:4) {
  param <- parameters[i]
  param_label <- param_labels[[i]]
  true_val <- true_values[i]
  
  # Extract data for this parameter
  if (param == "beta0") {
    clayton_data <- clayton_results$mean
    gaussian_data <- gaussian_results$mean
    direct_data <- direct_results$mean
  } else if (param == "beta1") {
    clayton_data <- clayton_results$mean1
    gaussian_data <- gaussian_results$mean1
    direct_data <- direct_results$mean1
  } else if (param == "scale") {
    clayton_data <- clayton_results$scale
    gaussian_data <- gaussian_results$scale
    direct_data <- direct_results$scale
  } else if (param == "shape") {
    clayton_data <- clayton_results$shape
    gaussian_data <- gaussian_results$shape
    direct_data <- direct_results$shape
  }
  
  # Combine data for boxplot
  data_list <- list(clayton_data, gaussian_data, direct_data)
  
  # Create boxplot
  boxplot(data_list, 
          col = colors,
          names = model_names,
          main = "",
          xlab = "",
          ylab = "",
          axes = FALSE,
          border = "black",
          lwd = 0.5)
  
  # Add true value line
  abline(h = true_val, col = "red", lty = 2, lwd = 0.8)
  
  # Add axes
  axis(1, at = 1:3, labels = model_names, cex.axis = 1)
  axis(2, cex.axis = 1)
  
  # Add parameter label
  mtext(param_label, side = 3, line = 0.5, cex = 1.2, font = 2)
  
  # Add box
  box()
}

dev.off()

cat(sprintf("Generated: %s\n", filename))
