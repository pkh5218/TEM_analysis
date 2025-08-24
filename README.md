# TEM–KE Solver (ERA5)

This repository contains a MATLAB implementation of the **Kuo–Eliassen (KE) Transformed Eulerian Mean (TEM) solver**, adapted for daily ERA5 reanalysis data.  
The code computes the TEM and EUM streamfunctions, Eliassen–Palm (EP) fluxes and their divergence and other forcing terms of the KE equation.  
The solver is optimized for **uniform pressure grids** to improve accuracy.

---

## Features

- Interpolates ERA5 daily data (`u`, `v`, `w`, `T`) to a uniform pressure grid (default: 12.5 hPa spacing).
- Computes:
  - TEM and EUM streamfunctions from integration of v* (`tempsi`, `eumpsi`)
  - EP fluxes (`EP*`) and their divergences (`EPD*`)
  - G-vector terms (`G*`) and their divergences (`GD*`)
  - Diabatic heating (`J`) and friction terms (`X`)
  - LU-decomposed KE equation solutions for each forcing term (`F*`)

---

## Code Structure

- **Main script:** `main_TEM_KE_daily.m`  
  Controls the workflow: reading ERA5, interpolation, flux/budget calculation, solving KE equation, and saving outputs.

- **Core functions:**
  - `compute_TEM_ERA5.m` – TEM & EUM streamfunctions from integration of v*
  - `compute_J.m` – compute diabatic heating as an residual
  - `compute_X.m` – compute friction as an residual
  - `compute_dfdp.m`, `compute_dfdy.m` – vertical and meridional derivatives
  - `compute_int.m` – vertical interpolation
  - `solve_PDE_LU_uniform.m` – LU decomposition for KE solver

---

## Output Variables

Each yearly output NetCDF file contains the following variables:

### 1. TEM-related (`TEM_YYYY.nc`)
- **FTEM** – Total TEM streamfunction from solving KE equation
- **FEUM** – Total EUM streamfunction from solving KE equation
- **tempsi** – TEM streamfunction from integration of v*
- **eumpsi** – Eulerian mean streamfunction from integration of v
- **FJ** – Contribution of diabatic heating to TEM streamfunction
- **FF** – Contribution of friction to TEM streamfunction
- **FEPDphi_QG**, **FEPDphi1** – Contribution of EP flux divergence (φ-component) to TEM streamfunction 
- **FEPDp_QG**, **FEPDp1**, **FEPDp2** – Contribution of EP flux divergence (p-component) to TEM streamfunction
- **FGDphi1**, **FGDphi2** – Contribution of G-vector divergence (φ-component) to TEM streamfunction 
- **FGDp1**, **FGDp2** – Contribution of G-vector divergence (p-component) to TEM streamfunction 

### 2. RHS Forcings (`RHS_YYYY.nc`)
- **HJ** – Forcing from diabatic heating  
- **HF** – Forcing from friction  
- **HEPDphi_QG**, **HEPDphi1** – Forcing from EP flux divergence (φ-component)  
- **HEPDp_QG**, **HEPDp1**, **HEPDp2** – Forcing from EP flux divergence (p-component)  
- **HGDphi1**, **HGDphi2** – Forcing from G-vector divergence (φ-component)  
- **HGDp1**, **HGDp2** – Forcing from G-vector divergence (p-component)  

### 3. Forcing Fields (`Forcing_YYYY.nc`)
- **J** – Diabatic heating (zonal mean)  
- **X** – Momentum forcing (zonal mean)  
- **EPDphi_QG**, **EPDphi1** – EP flux divergence (φ-component)  
- **EPDp_QG**, **EPDp1**, **EPDp2** – EP flux divergence (p-component)  
- **GDphi1**, **GDphi2** – G-vector divergence (φ-component)  
- **GDp1**, **GDp2** – G-vector divergence (p-component)  

### 4. Flux Variables (`Flux_YYYY.nc`)
- **EPphi_QG**, **EPphi1** – EP flux φ-components  
- **EPp_QG**, **EPp1**, **EPp2** – EP flux p-components  
- **Gphi1**, **Gphi2** – G-vector φ-components  
- **Gp1**, **Gp2** – G-vector p-components  
- **EMF** – Meridional eddy momentum flux (u'v')
- **EHF** – Meridional eddy heat flux (v'θ')
- **VEMF** – Vertical eddy momentum flux (u'w')
- **VEHF** – Vertical eddy heat flux (w'θ')

---

## Usage Notes

- Default pressure grid: **12.5 hPa spacing** from 1000 → 1 hPa.  
- Interpolation method: **pchip** 

Maintainer: [Pin-Chun Huang](mailto:pkh5218@psu.edu)  
GitHub: [https://github.com/pkh5218/TEM_analysis](https://github.com/pkh5218/TEM_analysis)
