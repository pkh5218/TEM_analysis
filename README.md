# TEM_analysis

This repository contains MATLAB scripts used to compute the **Transformed Eulerian Mean (TEM)** diagnostics based on the **Kuo–Eliassen equation**, using ERA5 reanalysis data on pressure levels.

## Contents

The repository includes the following MATLAB files:

| File Name               | Description |
|------------------------|-------------|
| `main_KE_solver_daily.m`      | **Main script** to solve the full TEM computation |
| `solve_TEM_KE_daily.m`        | Solves Kuo–Eliassen equation with interpolated data |
| `solve_PDE_LU_sparse.m`       | Solves the matrix system using LU decomposition (sparse) |
| `compute_TEM_ERA5.m`          | Obtain TEM streamfunction by integrating v* through ERA5 data. |
| `compute_var_3days.m`         | Calculates 3-day data (for computing F and J tendency) |
| `compute_J.m`                 | Computes diabatic heating term |
| `compute_F.m`                 | Computes frictional forcing |
| `compute_df.m` / `dfdp` / `dfdy` | Derivative operators |
| `compute_int.m`               | Compute interpolation of data on even-spaced pressure levels|
| `compute_S2.m`                | Computes static stability (S²) |

## ▶️ How to Run

1. Prepare ERA5 daily data: `u`, `v`, `omega`, and `T`, on pressure levels.
2. Modify paths and filenames as needed in `main_KE_solver_daily.m`.
3. Run the main script in MATLAB:

```matlab
main_KE_solver_daily
