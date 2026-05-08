# Regenerate Figure 6 from saved fitted models.
# This avoids refitting the spatial models in fit_spatial_models.R.

rm(list = ls())
library(GeoModels)

fits_file <- "data/processed/all_fitted_models_section148.rds"
figs_dir <- "results/figures"

if (!file.exists(fits_file)) {
  stop("Saved fitted models not found: ", fits_file,
       "\nRun src/fit_spatial_models.R first, or restore data/processed/all_fitted_models_section148.rds.")
}

if (!dir.exists(figs_dir)) {
  dir.create(figs_dir, recursive = TRUE)
}

fits <- readRDS(fits_file)

fit_direct <- fits$fit_direct
fit_gaussian <- fits$fit_gaussian
clayton_fits <- fits$clayton_fits
nu_values <- fits$nu_values
coords <- fits$coords

res_direct <- GeoResiduals(fit_direct)
vario_w <- GeoVariogram(data = res_direct$data, coordx = coords, maxdist = 8)

res_gaussian <- GeoResiduals(fit_gaussian)
vario_g <- GeoVariogram(data = res_gaussian$data, coordx = coords, maxdist = 8)

clayton_residuals <- lapply(clayton_fits, GeoResiduals)
clayton_variograms <- lapply(clayton_residuals, function(res) {
  GeoVariogram(data = res$data, coordx = coords, maxdist = 8)
})

clayton_aic <- sapply(clayton_fits, function(fit) 2 * 4 - 2 * fit$logCompLik)
best_clayton_idx <- which.min(clayton_aic)
best_clayton_nu <- nu_values[best_clayton_idx]

semivariogram_label <- "Semivariogram (scaled intensity squared)"

extract_pdf_page <- function(input, page, output) {
  pdfseparate <- Sys.which("pdfseparate")
  if (!nzchar(pdfseparate)) {
    stop("The 'pdfseparate' command is required to extract the plotted page.")
  }

  page_pattern <- tempfile(pattern = "combined_variogram_page_", fileext = "_%d.pdf")
  status <- system2(
    pdfseparate,
    args = c(
      "-f", as.character(page),
      "-l", as.character(page),
      shQuote(normalizePath(input, winslash = "/", mustWork = TRUE)),
      shQuote(normalizePath(page_pattern, winslash = "/", mustWork = FALSE))
    )
  )

  extracted <- sub("%d", as.character(page), page_pattern, fixed = TRUE)
  if (status != 0 || !file.exists(extracted)) {
    stop("Could not extract page ", page, " from ", input)
  }

  file.copy(extracted, output, overwrite = TRUE)
  unlink(extracted)
}

GeoCovariogram2 <- function(fitted, vario, variogram = TRUE, ...) {
  if (!inherits(fitted, "GeoFit")) {
    stop("Enter an object obtained of class GeoFit\n")
  }
  if (!inherits(vario, "GeoVariogram")) {
    stop("Enter an object of class GeoVariogram\n")
  }
  space <- !fitted$bivariate && !fitted$spacetime
  h <- seq(9e-5, 8, 0.01)
  if (space) {
    cc <- GeoCorrFct(
      x = h,
      corrmodel = fitted$corrmodel,
      covariance = TRUE,
      variogram = variogram,
      param = append(fitted$param, fitted$fixed),
      model = fitted$model
    )
    plot(
      cc,
      type = "l",
      xlab = "Distance (m)",
      ylab = semivariogram_label,
      col = "red",
      xlim = c(0, 8)
    )
    box()
  }
}

combined_variogram_file <- file.path(figs_dir, "combined_variograms_section148.pdf")
combined_variogram_tmp <- file.path(figs_dir, "combined_variograms_section148_twopage.pdf")

pdf(combined_variogram_tmp, width = 10, height = 6)
par(bg = "transparent")
plot(1, 1, type = "n", xlab = "", ylab = "", xlim = c(0, 8),
     ylim = c(0, 0.00028))
GeoCovariogram2(res_gaussian, vario = vario_g, pch = 20,
                ylim = c(0, 0.00025))
GeoCovariogram(res_direct, show.vario = TRUE, vario = vario_w, pch = 20,
               col = "blue", ylim = c(0, 0.00028), add.vario = TRUE)
GeoCovariogram(clayton_residuals[[best_clayton_idx]], show.vario = TRUE,
               vario = clayton_variograms[[best_clayton_idx]],
               pch = 20, col = "black", add.vario = TRUE)
box()
dev.off()

extract_pdf_page(combined_variogram_tmp, page = 2, output = combined_variogram_file)
unlink(combined_variogram_tmp)

cat("Saved: ", combined_variogram_file, "\n", sep = "")
cat("Best Clayton nu: ", best_clayton_nu, "\n", sep = "")
