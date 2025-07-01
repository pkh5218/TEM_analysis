# TEM_KE_Daily

This repository contains a MATLAB implementation of the **Kuo–Eliassen equation** in the Transformed Eulerian Mean (TEM) framework, designed to calculate the daily secondary circulation and forcing contributions from ERA5 reanalysis data.

##  Overview

The core script `solve_TEM_KE_daily.m` solves the quasi-diagnostic Kuo–Eliassen equation using zonally averaged heating, momentum, EP flux divergence, and G-vector divergence as forcings. The method handles daily resolution and layered vertical interpolation.

##  File Descriptions

| File | Description |
|------|-------------|
| `main_TEM_solver.m` | Main script. Iterates over time and layers to call the solver. Saves outputs to NetCDF. |
| `solve_TEM_KE_daily.m` | Core function that calculates all forcing terms and solves the Kuo–Eliassen equation via LU decomposition. |
| `save_nc_block.m` | Helper function to write output variables (with dimensions lat × lev × day) to NetCDF file. |
| `compute_TEM_ERA5.m` | Calculates TEM streamfunction from zonal-mean `v` and `v′θ′`. Outputs both EUM and TEM forms. |
| `compute_var_3days.m` | Extracts 3-day sequences of input variables for central-difference temporal derivatives. |
| `compute_S2.m` | Computes static stability (S²) from temperature and pressure. |
| `compute_J.m` | Estimates diabatic heating rate J via full temperature tendency equation including advection and radiative cooling. |
| `compute_int.m` | Vertically interpolates a 2D field to new pressure levels using spline method. |
| `compute_F.m` | Estimates frictional forcing from momentum budget terms and time derivatives. |
| `compute_dfdp.m` | Computes ∂f/∂p using variable grid spacing. |
| `compute_dfdy.m` | Computes ∂f/∂y using latitude-dependent spacing. |
| `compute_df.m` | Computes 1D finite difference spacing array for non-uniform grids. |
| `solve_PDE_LU_sparse.m` | Sparse matrix LU solver for elliptic PDE (Kuo–Eliassen form). Supports mixed derivatives and boundary conditions. |

---

##  Dependencies

- MATLAB (R2021a or later recommended)
- Input variables: ERA5 daily data of `u`, `v`, `w`, `T`, and `theta` on pressure levels
- Preprocessing: Zonal averaging must be done beforehand

## How to Use

1. **Prepare Inputs**:
    - ERA5 variables (`u`, `v`, `w`, `T`, `theta`) with size `[lon x lat x level x day]`
    - Latitude, pressure, and coordinate maps

2. **Run the driver**:
    ```matlab
    main_TEM_driver
    ```

3. **Outputs**:
    - All forcings (e.g., `FJ`, `FF`, `FEPphi1`, `FGp2`, etc.)
    - Daily `FTEM` streamfunction: the sum of all forcing responses
    - Saved to NetCDF using `save_nc_block`

## 📤 Output Variables

Each output variable has dimensions `[lat x lev x day]`, including:

- `FTEM`: Total streamfunction response
- `FJ`, `FF`: Diabatic heating & friction components
- `FEP*`, `FG*`: EP flux and G-vector components
- `tempsi`, `eumpsi`: Streamfunction from integrating v*

## 📞 Contact

Maintainer: [Pin-Chun Huang](mailto:pkh5218@psu.edu)  
GitHub: [https://github.com/pkh5218/TEM_analysis](https://github.com/pkh5218/TEM_analysis)
