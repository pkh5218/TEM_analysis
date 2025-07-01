# TEM_KE_Daily

This repository contains a MATLAB implementation of the **Kuo–Eliassen equation** in the Transformed Eulerian Mean (TEM) framework, designed to calculate the daily secondary circulation and forcing contributions from ERA5 reanalysis data.

## 📌 Overview

The core script `solve_TEM_KE_daily.m` solves the quasi-diagnostic Kuo–Eliassen equation using zonally averaged heating, momentum, EP flux divergence, and G-vector divergence as forcings. The method handles daily resolution and layered vertical interpolation.

## 🗂 File Structure

```
.
├── main_TEM_driver.m          # Main script to load data, loop over days, and save outputs
├── solve_TEM_KE_daily.m       # Main function solving the Kuo–Eliassen equation for one day
├── save_nc_block.m            # Utility to save output variables into NetCDF
├── src/
│   ├── compute_TEM_ERA5.m
│   ├── compute_var_3days.m
│   ├── compute_S2.m
│   ├── compute_J.m
│   ├── compute_int.m
│   ├── compute_F.m
│   ├── compute_dfdp.m
│   ├── compute_dfdy.m
│   └── compute_df.m
```

## ⚙️ Dependencies

- MATLAB (R2021a or later recommended)
- Input variables: ERA5 daily data of `u`, `v`, `w`, `T`, and `theta` on pressure levels
- Preprocessing: Zonal averaging must be done beforehand

## 🚀 How to Use

1. **Prepare Inputs**:
    - ERA5 variables (`u`, `v`, `w`, `T`, `theta`) with size `[lon x lat x level x day]`
    - Latitude, pressure, and coordinate maps (`phi`, `p_map`, `cosphi_map`, etc.)

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

---

**Note**: Please cite ERA5 appropriately if you use this code for publication.
