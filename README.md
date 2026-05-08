# Weibull Random Fields with Copula Constructions

This repository contains code, data, and precomputed results for fitting and
comparing Weibull random fields for spatially dependent positive data. The
application is based on LiDAR intensity measurements from a mining haul road,
and the simulation code reproduces the expanded study used in the manuscript.

The three Weibull random-field constructions are:

- chi-square transformation of Gaussian random fields;
- Gaussian copula Weibull random field;
- Clayton copula Weibull random field with a reflection-asymmetry parameter.

## Structure

```text
.
|-- src/                # R scripts for fitting, simulation, summaries, figures
|-- data/
|   |-- raw/            # Raw LiDAR data
|   `-- processed/      # Fitted models, summaries, and final simulation outputs
`-- results/
    `-- figures/        # Generated PDF figures
```

## Main Application Scripts

- `src/fit_spatial_models.R`: fits the application models and generates model
  comparison outputs.
- `src/bootstrap_standard_errors.R`: computes bootstrap standard errors for the
  application estimates.
- `src/generate_diagnostics.R`: generates exploratory LiDAR diagnostic plots.
- `src/generate_combined_variogram.R`: regenerates the combined semivariogram
  comparison figure from fitted models.
- `src/generate_contours.R`: generates contour plots comparing dependence
  constructions.
- `src/generate_realizations.R`: generates Clayton-Weibull random-field
  realizations.
- `src/generate_boxplots.R`: generates the original parameter-recovery boxplots.
- `src/generate_correlation_function.R`: generates Clayton-Weibull correlation
  function plots.

## Expanded Simulation Study

The expanded simulation uses a one-factor-at-a-time design around the baseline
scenario:

```text
S0: n = 800,  alpha = 0.06, kappa = 1,  nu = 1
S1: n = 400,  alpha = 0.06, kappa = 1,  nu = 1
S2: n = 1600, alpha = 0.06, kappa = 1,  nu = 1
S3: n = 800,  alpha = 0.08, kappa = 1,  nu = 1
S4: n = 800,  alpha = 0.12, kappa = 1,  nu = 1
S5: n = 800,  alpha = 0.06, kappa = 3,  nu = 1
S6: n = 800,  alpha = 0.06, kappa = 10, nu = 1
S7: n = 800,  alpha = 0.06, kappa = 1,  nu = 2
S8: n = 800,  alpha = 0.06, kappa = 1,  nu = 4
```

In each scenario, data are generated from the Clayton-Weibull random field and
three models are fitted: chi-square transformation, Gaussian copula, and
Clayton copula.

The main v2 scripts are:

- `src/simulation_config_v2.R`: scenario definitions and default parameters.
- `src/simulation_utils_v2.R`: shared simulation, fitting, bootstrap, summary,
  and plotting utilities.
- `src/run_simulation_study_v2.R`: main simulation study.
- `src/run_coverage_study_v2.R`: bootstrap standard-error and 95% coverage
  study.
- `src/summarize_simulation_study_v2.R`: combines chunk outputs and generates
  summary tables/figures.
- `src/run_all_v2.R`: convenience wrapper for the full v2 workflow.
- `src/run_parallel_external_v2.R`: Windows-friendly external parallel launcher
  using independent `Rscript` processes.

## Requirements

- R >= 4.0
- `GeoModels`
- Additional R packages used by the scripts include `fields`, `MASS`, `ggplot2`,
  `cubature`, `hypergeo`, `mvtnorm`, `gcKrig`, and `plot3D`.

## Usage

Run a smoke test first:

```powershell
Rscript src/run_all_v2.R --smoke
```

Run the final simulation study with 1000 Monte Carlo replicates:

```powershell
Rscript src/run_simulation_study_v2.R --reps 1000 --workers auto
```

Run the coverage and neighborhood-sensitivity study with 100 Monte Carlo
replicates and 50 bootstrap replicates per fitted model:

```powershell
Rscript src/run_coverage_study_v2.R --reps 100 --boot 50 --workers auto
```

Regenerate all summary CSV files and PDF figures from available chunk outputs:

```powershell
Rscript src/summarize_simulation_study_v2.R
```

For a short timing run:

```powershell
Rscript src/run_simulation_study_v2.R --scenario S0 --reps 10 --workers 1 --overwrite
Rscript src/run_coverage_study_v2.R --scenario S0 --reps 10 --boot 5 --workers 1 --m-values 2 --overwrite
Rscript src/summarize_simulation_study_v2.R
```

## Windows-Friendly External Parallel Runs

On Windows, the internal `parallel::parLapplyLB` backend can be unstable with
the current `GeoModels`/`progressr` stack. For real parallel execution on
Windows, use the external launcher. It starts independent `Rscript` processes,
each with `workers = 1`, and writes logs to `logs/`.

Simulation study split by scenario groups:

```powershell
Rscript src/run_parallel_external_v2.R --mode simulation --sim-reps 1000 --sim-chunk-size 25 --max-procs 8 --sim-groups "S0;S1;S2;S3;S4;S5;S6;S7;S8" --overwrite
```

Coverage study split by neighborhood size:

```powershell
Rscript src/run_parallel_external_v2.R --mode coverage --coverage-reps 100 --boot 50 --m-values 2,5,10 --coverage-chunk-size 5 --max-procs 3 --overwrite
```

Full external run:

```powershell
Rscript src/run_parallel_external_v2.R --mode all --sim-reps 1000 --coverage-reps 100 --boot 50 --sim-chunk-size 25 --coverage-chunk-size 5 --max-procs 8 --sim-groups "S0;S1;S2;S3;S4;S5;S6;S7;S8" --overwrite
```

From an interactive R console, pass arguments through `simulation_args_v2`:

```r
setwd("<path-to-repo>")
simulation_args_v2 <- c(
  "--mode", "all",
  "--sim-reps", "1000",
  "--coverage-reps", "100",
  "--boot", "50",
  "--sim-chunk-size", "25",
  "--coverage-chunk-size", "5",
  "--max-procs", "8",
  "--sim-groups", "S0;S1;S2;S3;S4;S5;S6;S7;S8",
  "--overwrite"
)
source("src/run_parallel_external_v2.R")
rm(simulation_args_v2)
```

## Outputs

Final generated CSV files are written to `data/processed/`:

- `simulation_design_v2.csv`
- `simulation_raw_results_v2.csv`
- `simulation_parameter_recovery_v2.csv`
- `simulation_model_selection_v2.csv`
- `coverage_raw_results_v2.csv`
- `coverage_summary_v2.csv`
- `neighborhood_sensitivity_v2.csv`

Final generated figures are written to `results/figures/`:

- `simulation_selection_frequencies_v2.pdf`
- `simulation_parameter_recovery_v2.pdf`
- `coverage_calibration_v2.pdf`

Intermediate chunk files are saved under
`data/processed/simulation_chunks_v2/` and
`data/processed/coverage_chunks_v2/`. These folders are ignored by git because
they are reproducible intermediate outputs.

## Notes

- The article-scale defaults are `B = 1000` for the main simulation study and
  `B = 100`, `K = 50` for the coverage study.
- Coverage intervals use `estimate +/- qnorm(0.975) * SE`; there is no extra
  multiplicative factor.
- Precomputed application outputs and final simulation summaries are included
  so the manuscript tables and figures can be reproduced without rerunning the
  full Monte Carlo study.
