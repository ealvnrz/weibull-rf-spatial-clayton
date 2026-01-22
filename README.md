# Weibull Random Fields with Copula Constructions

This repository contains code and results for fitting Weibull random fields using different constructions: χ² transformation, Gaussian copula, and Clayton copula.

## Structure

```
.
├── src/                # R scripts for model fitting and analysis
├── data/
│   ├── raw/            # Raw data files
│   └── processed/      # Processed data and fitted models
└── results/
    └── figures/        # Generated figures (PDF)
```

## Main Scripts

- `fit_spatial_models.R`: Fits all spatial models (Independence, χ² transformation, Gaussian copula, Clayton copula with different nu values)
- `bootstrap_standard_errors.R`: Computes bootstrap standard errors for parameter estimates
- `generate_diagnostics.R`: Generates exploratory diagnostic plots
- `generate_contours.R`: Generates contour plots comparing different constructions
- `generate_realizations.R`: Generates realizations of Clayton-Weibull random fields
- `generate_boxplots.R`: Generates boxplots comparing parameter estimates
- `generate_correlation_function.R`: Generates correlation function plots
- `simulate_*.R`: Simulation studies for parameter recovery

## Data

The main data file is `data/raw/LAS_sec_148.csv`, containing LiDAR intensity data with spatial coordinates.

## Results

- Fitted models are saved in `data/processed/` as RDS files
- Model comparison and parameter estimates are saved as CSV files in `data/processed/`
- Figures are saved in `results/figures/` as PDF files

## Requirements

- R (>= 4.0)
- GeoModels package
- Additional packages: fields, MASS, ggplot2, cubature, hypergeo, mvtnorm, gcKrig, plot3D

## Usage

1. Run `fit_spatial_models.R` to fit all models
2. Run `bootstrap_standard_errors.R` to compute standard errors (requires fitted models)
3. Run individual `generate_*.R` scripts to create specific figures

## Note

This repository is a complement to the associated research paper. All results have been pre-computed and are included in the repository.
