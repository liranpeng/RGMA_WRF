
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
---

## Authors

This simulation pipeline was developed by Liran Peng and is part of RGMA project.
