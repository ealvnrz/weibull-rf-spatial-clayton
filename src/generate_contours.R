# Generate contour plots: comparison of three Weibull constructions
# Nine contour plots with rho ∈ {0.1,0.5,0.9} and nu ∈ {1,2,6}

rm(list = ls())
library(cubature)
library(hypergeo)
library(mvtnorm)
library(gcKrig)
library(plot3D)

# Parameters
kappa <- 5
mui <- 1
muj <- 1
x_range <- seq(-3, 3, 0.1)
rho_values <- c(0.1, 0.5, 0.9)
nu_values <- c(1, 2, 6)

# Helper functions
bivariate_gaussian <- function(zi, zj, rho) {
  aa <- 1 / (sqrt(1 - rho^2) * 2 * pi)
  tst <- exp(-0.5 * (1 - rho^2)^(-1) * (zi^2 + zj^2 - 2 * rho * zi * zj))
  return(aa * tst)
}

clayton_uniform <- function(zi, zj, rho, nu) {
  k <- 0
  a1 <- 0
  res0 <- 0.0
  RR <- 0.0
  pp1 <- 0.0
  
  x <- zi^(2/nu)
  y <- zj^(2/nu)
  nu2 <- nu/2
  cc <- nu2 + 1
  c <- ((1 - rho^2)^(-cc))
  aux <- (rho^2) * x * y
  aux1 <- (rho^2) * (1 - x) * (1 - y)
  
  while (k <= 10000) {
    pp1 <- (nu2 - 2*(cc + k)) * log1p(-aux) + 
           log(Re(hypergeo(nu2 - (cc + k), nu2 - (cc + k), nu2, aux)))
    bb1 <- pp1 + k * log(aux1) + 2 * (lgamma(cc + k) - lgamma(cc)) - 
           lgamma(k + 1) - lgamma(1 + k)
    a1 <- a1 + exp(bb1)
    RR <- a1 / c
    if (abs(RR - res0) < 1e-40) break
    else res0 <- RR
    k <- k + 1
  }
  return(RR)
}

clayton_weibull_copula <- function(zi, zj, rho, nu, kappa) {
  const <- (gamma(1 + 1/kappa))^(-1)
  ep1 <- (-(zi/const)^kappa)
  ep2 <- (-(zj/const)^kappa)
  x <- 1 - exp(ep1)
  y <- 1 - exp(ep2)
  p1 <- kappa * zi^(kappa - 1) * exp(ep1) / (const^kappa)
  p2 <- kappa * zj^(kappa - 1) * exp(ep2) / (const^kappa)
  return(clayton_uniform(x, y, rho, nu) * p1 * p2)
}

weibull_copula_uniform <- function(wi, wj, rho, nu, kappa) {
  k <- (gamma(1 + 1/kappa))^(-1)
  zi <- qweibull(wi, kappa, k)
  zj <- qweibull(wj, kappa, k)
  a <- clayton_weibull_copula(zi, zj, rho, nu, kappa)
  c1 <- k * (-log(1 - wi))^(1/kappa - 1) / (kappa * (1 - wi))
  c2 <- k * (-log(1 - wj))^(1/kappa - 1) / (kappa * (1 - wj))
  return(a * abs(c1 * c2))
}

clayton_weibull_gaussian_scale <- function(zi, zj, rho, nu, kappa) {
  x <- pnorm(zi)
  y <- pnorm(zj)
  p1 <- dnorm(zi)
  p2 <- dnorm(zj)
  return(weibull_copula_uniform(x, y, rho, nu, kappa) * p1 * p2)
}

chi2_weibull_bivariate <- function(y1, y2, mui, muj, shape, rho) {
  k <- (gamma(1 + 1/shape))^(-1)
  ci <- mui
  cj <- muj
  ui <- y1/ci
  uj <- y2/cj
  a <- 1 - rho^2
  z <- 2 * rho * (ui * uj)^(shape/2) * k^(-shape) / a
  
  A <- 2 * log(shape) + (-2 * shape) * log(k) + 
       (shape - 1) * log(y1 * y2) - (shape) * log(ci * cj) - log(a)
  B <- -k^(-shape) * (ui^(shape) + uj^(shape)) / a
  res <- A + B + log(besselI(z, 0))
  return(exp(res))
}

chi2_weibull_uniform <- function(w1, w2, mui, muj, shape, rho) {
  k <- (gamma(1 + 1/shape))^(-1)
  y1 <- qweibull(w1, shape, k * mui)
  y2 <- qweibull(w2, shape, k * muj)
  a <- chi2_weibull_bivariate(y1, y2, mui, muj, shape, rho)
  c1 <- mui * k * (-log(1 - w1))^(1/shape - 1) / (shape * (1 - w1))
  c2 <- muj * k * (-log(1 - w2))^(1/shape - 1) / (shape * (1 - w2))
  return(a * abs(c1 * c2))
}

chi2_weibull_gaussian_scale <- function(z1, z2, mui, muj, shape, rho) {
  res <- chi2_weibull_uniform(pnorm(z1), pnorm(z2), mui, muj, shape, rho) * 
         (dnorm(z1) * dnorm(z2))
  return(res)
}

# Output directory
output_dir <- "results/figures"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Generate all 9 contour plots
for (rho in rho_values) {
  for (nu in nu_values) {
    # Vectorize functions
    bivariate_gaussian_vec <- Vectorize(bivariate_gaussian, c("zi", "zj"))
    clayton_weibull_gauss_vec <- Vectorize(clayton_weibull_gaussian_scale, c("zi", "zj"))
    chi2_weibull_gauss_vec <- Vectorize(chi2_weibull_gaussian_scale, c("z1", "z2"))
    
    # Generate density matrices
    gaussian_cop <- outer(x_range, x_range, bivariate_gaussian_vec, rho)
    clayton_cop <- outer(x_range, x_range, clayton_weibull_gauss_vec, rho, nu, kappa)
    chi2_cop <- outer(x_range, x_range, chi2_weibull_gauss_vec, mui, muj, kappa, rho)
    
    # Generate filename
    filename <- sprintf("weibull_contour_rho%.1f_nu%d.pdf", rho, nu)
    filepath <- file.path(output_dir, filename)
    
    # Create plot
    pdf(filepath, width = 8, height = 6)
    par(mar = c(4, 4, 1, 1), oma = c(0, 0, 0, 0))
    
    colors <- c("blue", "red", "green")
    contour(x_range, x_range, gaussian_cop, 
            lwd = 2, col = colors[1], 
            levels = seq(0.02, 0.4, 0.04),
            xlab = "", ylab = "", axes = FALSE)
    contour(x_range, x_range, clayton_cop, 
            add = TRUE, lwd = 2, col = colors[2], 
            levels = seq(0.02, 0.4, 0.04))
    contour(x_range, x_range, chi2_cop, 
            add = TRUE, lwd = 2, col = colors[3], 
            levels = seq(0.02, 0.4, 0.04))
    
    axis(1, at = seq(-3, 3, 1), labels = seq(-3, 3, 1))
    axis(2, at = seq(-3, 3, 1), labels = seq(-3, 3, 1))
    box()
    
    dev.off()
    
    cat(sprintf("Generated: %s\n", filename))
  }
}

cat("\nAll 9 contour plots generated successfully.\n")
