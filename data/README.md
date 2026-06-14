# Data Documentation

## Provenance

`data/raw/LAS_sec_148.csv` contains a LiDAR-derived point extract for mining
haul-road section 148. The file is the real-data input used in the application
workflow for the accepted article:

> Weibull random fields for spatially dependent positive data: an application to mining haul roads

The data are included so the fitted application models, tables, and figures can
be reproduced.

## Column Dictionary

- `x`: raw horizontal coordinate from the LiDAR extract.
- `y`: raw vertical coordinate from the LiDAR extract.
- `x_scaled`: scaled version of `x`, used for normalized spatial summaries.
- `y_scaled`: scaled version of `y`, used for normalized spatial summaries.
- `z_scaled`: scaled elevation or height coordinate.
- `intensity`: raw LiDAR return intensity.
- `intensity_scaled`: scaled positive intensity response used in the Weibull
  random-field application.
- `r`: scaled red channel associated with the point.
- `g`: scaled green channel associated with the point.
- `b`: scaled blue channel associated with the point.
- `Cluster`: haul-road section or cluster identifier.
- `Cluster_R`: companion cluster identifier retained from the preprocessing
  workflow.

The raw coordinates `x` and `y` preserve the coordinate scale of the LiDAR
extract. The scaled coordinate columns are normalized derivatives used for
inspection and preprocessing; the application fitting scripts use `x` and `y`
as spatial coordinates and `intensity_scaled` as the response.
