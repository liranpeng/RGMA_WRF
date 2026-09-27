
# High-Resolution WRF-Chem Simulation: Domain Nesting and Aerosol Initialization Workflow

This repository documents the step-by-step procedure for conducting a nested WRF-Chem simulation over the ARM Eastern North Atlantic (ENA) site. The workflow includes three main stages:

- **Step 1:** Spin-up simulation for nested domains 1, 2, and 3 (no chemistry)
- **Step 2:** High-resolution domain 4 initialization and ndown step
- **Step 3:** Aerosol initialization via boundary and restart file modification

---

## Step 1: Domain 1–2–3 Spin-up (No Chemistry)

1. Copy base directory and navigate into it:
   ```bash
   cp -r save_em_real Domain123_20170713_nochem
   cd Domain123_20170713_nochem
   ```

2. Copy WPS output (`met_em*`) into working directory:
   ```bash
   cp /pscratch/sd/h/heroplr/wrf_scratch/WPS/met_em* .
   ```

3. Generate WRF boundary and initial conditions:
   ```bash
   sbatch submit_real.sh
   ```

4. Submit the containerized WRF simulation:
   ```bash
   sbatch submit_container.sh
   ```

---

## Step 2: Domain 4 High-Resolution Initialization

### 2.1 Chemistry-on Initialization

1. Create and enter the new domain directory:
   ```bash
   cp -r save_em_real Domain4_20170713_chem
   cd Domain4_20170713_chem
   ```

2. Copy high-resolution met_em files and rename:
   ```bash
   cp /pscratch/sd/h/heroplr/wrf_scratch/WPS/met_em.d03.2017* .
   nohup ./rename_met_em_d03_to_d01.sh >& rename_met_em_d03_to_d01.out &
   ```

3. Prepare `namelist.input`:
   ```bash
   cp namelist.input_domain4_init namelist.input
   ```

4. Copy `wrfout` files from Domain 123 simulation and rename:
   ```bash
   cp ../Domain123_20170713_nochem/wrfout_d03* .
   nohup ./rename_wrfout.sh >& rename_wrfout.out &
   ```

5. Run real.exe for Domain 4:
   ```bash
   sbatch submit_real.sh
   ```

6. Prepare for `ndown.exe`:
   ```bash
   cp wrfinput_d02 wrfndi_d02
   sbatch submit_ndown.sh
   ```

   This step creates `wrfbdy_d02`.

### 2.2 No-Aerosol Run for Restart File

1. Create new directory for clean-aerosol simulation:
   ```bash
   cd ..
   cp -r save_em_real Domain4_20170713_chem_noaero
   cd Domain4_20170713_chem_noaero
   ```

2. Copy boundary and initial condition files from the previous step:
   ```bash
   cp ../Domain4_20170713_chem/wrfinput_d02 wrfinput_d01
   cp ../Domain4_20170713_chem/wrfbdy_d02 wrfbdy_d01
   cp ../namelist.input_d04 namelist.input
   ```

3. Run simulation to generate restart files:
   ```bash
   sbatch submit_container.sh
   ```

   This generates files like:
   ```
   wrfrst_d01_2017-07-13_04:00:00
   ```

---

## Step 3: Aerosol Initialization via Restart and Boundary Update

*Optional step for clarity: perform this step in a new folder.*

1. Copy template directory:
   ```bash
   cp -r save_em_real Domain4_20170713_chem_aero
   cd Domain4_20170713_chem_aero
   ```

2. Copy necessary files:
   ```bash
   cp ../Domain4_20170713_chem_noaero/wrfbdy_d01 .
   cp ../Domain4_20170713_chem_noaero/wrfrst_d01_2017-07-13_04:00:00 .
   ```

3. Inject aerosols into boundary and restart files:
   ```bash
   ./update_wrfbdy.py
   cp wrfbdy_d01_so4added.nc wrfbdy_d01

   ./update_wrfrst.py
   cp wrfrst_d01_2017-07-13_04:00:00_so4.nc wrfrst_d01_2017-07-13_04:00:00
   ```

4. Set namelist and launch the final chemistry-aware simulation:
   ```bash
   cp namelist.input_domain4_restart namelist.input
   sbatch submit_container.sh
   ```

---

## Notes
- To go through this document, please first finish WPS 
- `rename_*` saved under Script/
- `update_*` saved under Run_WRF/AddAerosol/ 
- Each step is designed to preserve clarity, reproducibility, and modularity across different experiment branches.
- Adjust `namelist.input` files accordingly for different simulation stages (`init`, `restart`, `noaero`, etc.).
- Use this setup steps, I have finished four cases, 2016-07-01, 2017-07-01, 2017-07-13, and 2017-07-18

# 📂 WRF Data Backup – Directory Structure

This directory contains backup data for WRF simulations on multiple dates.  
Each date folder contains **two subfolders**:  

- `Domain123` – wrf output restart namelist input boudnary files for domains 1, 2, and 3  
- `Domain4` – wrf output restart namelist input boudnary files domain 4  

---

## 📁 Path Map

```
/pscratch/sd/h/heroplr/Data_To_Yan/WRF_Data_Backup/
│
├── 20170630/
│   ├── Domain123/
│   └── Domain4/
│
├── 20170701/
│   ├── Domain123/
│   └── Domain4/
│
├── 20170713/
│   ├── Domain123/
│   └── Domain4/
│
└── 20170718/
    ├── Domain123/
    └── Domain4/
```

---

# 📊 Paper Figures – ARG vs ML Aerosol Activation

`Paper_Figures/` holds one self-contained bundle per figure of the ARG-vs-ML
activation study (2017-07-14 cases, `1aer`/`3aer` × `org`/`mid`). Each bundle has
its own `README.md`, the plotting script(s) unchanged from the producing run, a
`run.sh` wrapper that picks the conda python, the inputs the script reads, an
`original_figure/` reference copy, and a `provenance/` folder with the scripts
and logs that generated the inputs.

| Bundle | Figure | Inputs shipped |
|--------|--------|----------------|
| `arg_relerr_map_filtered/` | 2×4 map of the ARG scheme's median relative error in activated fraction vs. the Fortran parcel model, binned in (w, N_accum), (w, N_Aitken), (N_accum, N_Aitken), (w, p); v2/v3 variants split Aitken / accumulation. | per-case `online/*.npz` (**not in git**, see below) |
| `parcel_vs_ARG_vs_ML_ALLCASES/` | 2×3 pooled hexbin: parcel-model truth vs. the ML fraction WRF applied (FN11/FN21) vs. the ARG fraction WRF computed (FN13/FN23). | per-case `online/*.npz` (**not in git**) |
| `arg_vs_ml_consequences_rwp_v8/` | 5×2 time series 09:00–13:00 UTC 15 Jul 2017, ML-ACT solid vs ARG-ACT (WRF_dm_v2) dashed: N_acc, N_ait, CCN(0.05%), N_d, r_v, LWP, rain-water path, cloudy-column fraction. v12 variant included. | `*_series_CCN2.npz` domain-mean cache (in git) |
| `dCRE_combined_4panel_MLvsARGv2/` | 2×2 shortwave dCRE decomposition, ML-ACT vs ARG-ACT, org / mid: half-hourly stacked components and time-mean bars (total, CF, A_Nc, A_LWP, A_re, A_cov). v3/v4 variants and 30-min ARG−ML tables included. | `dcre_csv/*.csv` (in git) |

Replotting: `cd Paper_Figures/<bundle> && ./run.sh`. The two dCRE and
consequences bundles are pure replots from the shipped CSV / cache files and run
in seconds on a login node. The two parcel-model bundles need the per-case
`20170714_*/online/{parcel_inputs,parcel_shard_0,diag_fn_ml_offline}.npz`
files (~650 MB per case, 2.6 GB total), which exceed GitHub's file-size limit
and are therefore **not committed**. They live on Stampede3 at

    /scratch/07088/tg863871/RGMA_WRF/Paper_Figures/<bundle>/20170714_*/online/
    /scratch/07088/tg863871/Perlm_Backup/WRF_stam3_ml_mixout_bce_lr1e-4_cosine_more_fix/test/<case>/online/

Copy them beside the script before running those two bundles. Large rasterised
EPS exports are likewise left out; the PNG (and PDF where produced) versions of
every figure are committed.

Source model runs (multi-TB, not included):

    ML-ACT : /scratch/07088/tg863871/Perlm_Backup/WRF_stam3_ml_mixout_bce_lr1e-4_cosine_more_fix/test/20170714_{1,3}aer_{org,mid}
    ARG-ACT: /scratch/07088/tg863871/WRF_dm_v2/test/20170714_{org,mid}_{1,3}aer_ARG

The corresponding model source is preserved as the `WRF_ML` and `WRF_ARG`
branches of https://github.com/liranpeng/WRF.

---

## Authors

This simulation pipeline was developed by Liran Peng and is part of RGMA project.
