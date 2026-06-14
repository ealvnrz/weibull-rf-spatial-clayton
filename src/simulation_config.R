# Configuration for the expanded simulation study.

simulation_defaults <- list(
  reps = 1000L,
  sim_reps = 1000L,
  coverage_reps = 100L,
  boot = 50L,
  workers = "auto",
  chunk_size = 25L,
  sim_chunk_size = 25L,
  coverage_chunk_size = 5L,
  seed_base = 7689L,
  overwrite = FALSE,
  smoke = FALSE,
  scenario = "all",
  m_values = "2,5,10",
  verbose = TRUE
)

fixed_parameters <- list(
  model = "Weibull",
  corrmodel = "GenWend_Matern",
  mean = 0.2,
  mean1 = -0.2,
  smooth = 0,
  power2 = 1 / 3.5,
  nugget = 0,
  sill = 1,
  neighb = 2,
  lower_bound = 20,
  upper_bound = 20
)

simulation_scenarios <- data.frame(
  scenario_id = paste0("S", 0:8),
  scenario_group = c(
    "baseline",
    "sample_size", "sample_size",
    "spatial_dependence", "spatial_dependence",
    "marginal_skewness", "marginal_skewness",
    "clayton_asymmetry", "clayton_asymmetry"
  ),
  description = c(
    "Baseline",
    "Smaller sample size",
    "Larger sample size",
    "Stronger spatial dependence",
    "Even stronger spatial dependence",
    "Lower marginal skewness",
    "Much lower marginal skewness",
    "Symmetric Clayton case",
    "Opposite/asymmetric tail behavior"
  ),
  n = c(800L, 400L, 1600L, 800L, 800L, 800L, 800L, 800L, 800L),
  alpha = c(0.06, 0.06, 0.06, 0.08, 0.12, 0.06, 0.06, 0.06, 0.06),
  kappa = c(1, 1, 1, 1, 1, 3, 10, 1, 1),
  nu = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 4L),
  stringsAsFactors = FALSE
)

model_specs <- data.frame(
  model_id = c("direct", "gaussian", "clayton"),
  model_label = c("Chi-square transformation", "Gaussian copula", "Clayton copula"),
  stringsAsFactors = FALSE
)

parameter_names <- c("mean", "mean1", "scale", "shape")
