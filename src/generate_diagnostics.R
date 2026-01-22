# Generate diagnostic plots: quilt plot, histogram with Weibull fit, spatial scatterplot, and empirical variogram
# Exploratory analysis for the application section

rm(list = ls())
library(GeoModels)
library(fields)
library(MASS)
library(ggplot2)

# Load data
data_file <- "data/raw/LAS_sec_148.csv"
data <- read.csv(data_file)

# Extract coordinates and intensity
coords <- cbind(data$x, data$y)
intensity <- data$intensity_scaled

# Filter valid observations
valid <- is.finite(coords[, 1]) & is.finite(coords[, 2]) & 
         is.finite(intensity) & (intensity > 0)
coords <- coords[valid, ]
intensity <- intensity[valid]

cat("Loaded", length(intensity), "valid observations\n")

# Output directory
output_dir <- "results/figures"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Quilt plot
cat("Generating quilt plot...\n")
pdf(file.path(output_dir, "quilt_plot_section148.pdf"), width = 6, height = 6)
par(mar = c(5.1, 4.1, 4.1, 2.1))
quilt.plot(coords, intensity, 
           xaxt = "n", yaxt = "n",
           #nrow = 200, ncol = 200,
           col = tim.colors(64))
dev.off()
cat("Saved: quilt_plot_section148.pdf\n")

# Histogram with Weibull fit
cat("Fitting Weibull distribution...\n")
intensity_positive <- intensity[intensity > 0]

fit_weibull <- MASS::fitdistr(intensity_positive, "weibull")
shape <- as.numeric(fit_weibull$estimate["shape"])
scale <- as.numeric(fit_weibull$estimate["scale"])
cat("Weibull parameters: shape =", round(shape, 3), ", scale =", round(scale, 3), "\n")

# Extract breaks from ggplot to match exactly
# Create temporary ggplot to extract breaks
data_plot <- data.frame(intensity = intensity_positive)
p_temp <- ggplot(data_plot, aes(x = intensity)) +
  geom_histogram(aes(y = after_stat(density)), bins = 16, 
                fill = "gray85", color = "black", linewidth = 0.3)
# Extract breaks from ggplot
gg_build <- ggplot_build(p_temp)
gg_breaks <- gg_build$data[[1]]$xmin
gg_breaks <- c(gg_breaks, max(gg_build$data[[1]]$xmax))  # Add the last break

# R base version (matching ggplot style exactly)
# Use the exact breaks from ggplot
breaks <- gg_breaks

x_seq <- seq(min(intensity_positive), max(intensity_positive), length.out = 200)
density_weibull <- dweibull(x_seq, shape = shape, scale = scale)

# Calculate ylim to include both histogram and density curve
h <- hist(intensity_positive, breaks = breaks, freq = FALSE, plot = FALSE)
max_hist <- max(h$density, na.rm = TRUE)
max_density <- max(density_weibull, na.rm = TRUE)
ylim_max <- max(max_hist, max_density) * 1  # Add 5% padding

pdf(file.path(output_dir, "histogram_weibull_fit_section148.pdf"), width = 6, height = 6)
par(mar = c(5.1, 4.1, 4.1, 2.1))
hist(intensity_positive, breaks = breaks, freq = FALSE,
     col = "gray85", border = "black", lwd = 0.3,
     xlab = "Intensity", ylab = "Density",
     main = "", axes = FALSE, ylim = c(0, ylim_max))
# Custom axis labels - show fewer values
x_range <- range(intensity_positive)
x_ticks <- pretty(x_range, n = 4)  # ~4 ticks on x-axis
y_ticks <- pretty(c(0, ylim_max), n = 4)  # ~4 ticks on y-axis
axis(1, at = x_ticks, labels = x_ticks)
axis(2, at = y_ticks, labels = y_ticks)
lines(x_seq, density_weibull, col = "black", lwd = 1.2)
box()
dev.off()
cat("Saved: histogram_weibull_fit_section148.pdf\n")

# Spatial scatterplot to check for asymmetry (in Gaussian scale)
cat("Creating spatial scatterplot to check for asymmetry...\n")
# Fit independence model to get Gaussian transformation
model_temp <- "Weibull"
corrmodel_temp <- "GenWend_Matern"
power2_temp <- 1/4
nugget_temp <- 0
smooth_temp <- 0

temp_fit <- GeoFit(
  data = intensity_positive,
  coordx = coords[intensity > 0, ],
  model = model_temp,
  corrmodel = corrmodel_temp,
  start = list(mean = 0, shape = shape),
  fixed = list(power2 = power2_temp, nugget = nugget_temp, smooth = smooth_temp, scale = 0.0001),
  optimizer = "nlminb",
  lower = list(mean = -10, shape = 0.5),
  upper = list(mean = 10, shape = 20),
  likelihood = "Marginal",
  type = "Independence"
)

# Get Gaussian-transformed data
gaussian_data <- GeoPit(temp_fit, type = "Gaussian")

# High resolution PDF
pdf(file.path(output_dir, "spatial_scatterplot_section148.pdf"), 
    width = 6, height = 6, 
    pointsize = 12,
    useDingbats = FALSE)
par(mar = c(4, 4, 2, 2))
GeoScatterplot(gaussian_data$data, coords[intensity > 0, ], 
               neighb = 4,
               xlim = c(-4, 4),
               ylim = c(-4, 4))
dev.off()
cat("Saved: spatial_scatterplot_section148.pdf\n")

