# Weibull random fields for spatially dependent positive data

Code, data, and precomputed outputs for the accepted article:

> Weibull random fields for spatially dependent positive data: an application to mining haul roads

Authors: Eloy Alvarado, Christian Caamaño-Carrillo, and David R. Godoy.

This repository is intended to reproduce the computational material used in the
accepted manuscript. It contains R scripts, the LiDAR haul-road data extract,
fitted model objects, processed tables, and generated figures.

## Repository Structure

```text
.
|-- reproduce_all.R        # Reproduction workflow
|-- src/                   # R scripts for application, simulations, and figures
|-- data/
|   |-- raw/               # LiDAR section data
|   `-- processed/         # Fitted objects, tables, and precomputed summaries
`-- results/
    `-- figures/           # Generated PDF figures
```

## Data

The real-data application uses `data/raw/LAS_sec_148.csv`, a LiDAR intensity
extract for haul-road section 148. See `data/README.md` for provenance notes,
column definitions, scaling details, and data-use cautions.

## Dependencies

Use R 4.0 or later. The scripts require these R packages:

- `GeoModels`
- `fields`
- `MASS`
- `ggplot2`
- `cubature`
- `hypergeo`
- `mvtnorm`
- `gcKrig`
- `plot3D`
- `processx` for `src/run_parallel_external.R`

## Application Workflow

The application workflow is centered on haul-road section 148:

```powershell
Rscript src/fit_spatial_models.R
Rscript src/bootstrap_standard_errors.R --filter-only
Rscript src/generate_diagnostics.R
Rscript src/generate_combined_variogram.R
```

The active application model-comparison code follows the accepted manuscript
and uses Clayton copula values `nu = 1, 2, 4` only.

Primary application outputs:

- `data/processed/model_comparison_section148.csv`
- `data/processed/parameter_estimates_section148.csv`
- `data/processed/bootstrap_standard_errors_section148.csv`
- `data/processed/all_fitted_models_section148.rds`
- `results/figures/combined_variograms_section148.pdf`
- `results/figures/empirical_variogram_section148.pdf`
- `results/figures/quilt_plot_section148.pdf`
- `results/figures/histogram_weibull_fit_section148.pdf`
- `results/figures/spatial_scatterplot_section148.pdf`

## Simulation Workflow

The simulation study uses a one-factor-at-a-time design around a baseline
Clayton-Weibull random field. Scenario definitions are in
`src/simulation_config.R`, and shared utilities are in `src/simulation_utils.R`.

Run the full Monte Carlo study:

```powershell
Rscript src/run_simulation_study.R --reps 1000 --workers auto
```

Run the coverage and neighborhood-sensitivity study:

```powershell
Rscript src/run_coverage_study.R --reps 100 --boot 50 --workers auto
```

Regenerate summaries and simulation figures from chunk outputs, or from the
precomputed raw CSV outputs included in this repository:

```powershell
Rscript src/summarize_simulation_study.R
```

On Windows, the external launcher can run independent `Rscript` jobs:

```powershell
Rscript src/run_parallel_external.R --mode all --sim-reps 1000 --coverage-reps 100 --boot 50 --sim-chunk-size 25 --coverage-chunk-size 5 --max-procs 8 --sim-groups "S0;S1;S2;S3;S4;S5;S6;S7;S8" --overwrite
```

The full Monte Carlo and bootstrap workflows are computationally expensive.
The repository includes final processed summaries so manuscript tables and
figures can be reproduced without rerunning all replicates.

Simulation outputs:

- `data/processed/simulation_design.csv`
- `data/processed/simulation_raw_results.csv`
- `data/processed/simulation_parameter_recovery.csv`
- `data/processed/simulation_model_selection.csv`
- `data/processed/coverage_raw_results.csv`
- `data/processed/coverage_summary.csv`
- `data/processed/neighborhood_sensitivity.csv`
- `results/figures/simulation_selection_frequencies.pdf`
- `results/figures/simulation_parameter_recovery.pdf`
- `results/figures/coverage_calibration.pdf`

## Figure Workflows

The figure scripts intentionally use different Clayton parameter sets depending
on the manuscript figure:

- Application model comparison: `nu = 1, 2, 4`
- Contour-comparison figure: `nu = 1, 2, 6`
- Random-field realizations figure: `nu = 1, 2, 4`

Run contour and realization figures directly when needed:

```powershell
Rscript src/generate_contours.R
Rscript src/generate_realizations.R
Rscript src/generate_correlation_function.R
Rscript src/generate_boxplots.R
```

The contour figure preserves `nu = 6` by design and is separate from the
real-data application model comparison.

## Reproduction

Run a quick non-mutating smoke check:

```powershell
Rscript reproduce_all.R --smoke
```

This checks required files and packages and verifies that the simulation
configuration and utility code load. It does not regenerate outputs.

Run the broader reproduction workflow only when intentionally rebuilding
outputs, including contour, realization, and correlation figures:

```powershell
Rscript reproduce_all.R --full
```
