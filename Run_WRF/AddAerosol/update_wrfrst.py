from datetime import datetime
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.units import units
import numpy as np
import xarray as xr
from scipy import interpolate

# Constants
rho_air = 1.2               # air density in kg/m³
rho_sulfate = 1800          # kg/m³
pi = np.pi
Ntot_to_Nmix = 1 / rho_air  # convert #/m³ to #/kg_air

# Geometric diameters (in meters)
Dg_ai = 10e-9               # Aitken mode
Dg_aj = 70e-9               # Accumulation mode

# Geometric standard deviations
sigma_ai = 1.7
sigma_aj = 2.0

# Compute volume mean diameters (lognormal)
Dv_ai = Dg_ai * np.exp(1.5 * (np.log(sigma_ai))**2)
Dv_aj = Dg_aj * np.exp(1.5 * (np.log(sigma_aj))**2)

# Particle volume (m³)
Vp_ai = (pi/6) * Dv_ai**3
Vp_aj = (pi/6) * Dv_aj**3

# Conversion factor from number to mass mixing ratio (kg/kg)
conversion_ai = Ntot_to_Nmix * rho_sulfate * Vp_ai
conversion_aj = Ntot_to_Nmix * rho_sulfate * Vp_aj

# Mass per particle (kg)
so4_i = rho_sulfate * Vp_ai
so4_j = rho_sulfate * Vp_aj

print(f"so4_i (Aitken particle mass): {so4_i:.4e} kg")
print(f"so4_j (Accumulation particle mass): {so4_j:.4e} kg")

# Open WRF domain files
path = '/pscratch/sd/h/heroplr/wrf_scratch/WPS/20170713/'
fdomain01 = xr.open_dataset(path + 'met_em.d01.2017-07-13_07:10:00.nc')
fdomain02 = xr.open_dataset(path + 'met_em.d02.2017-07-13_07:10:00.nc')
fdomain03 = xr.open_dataset(path + 'met_em.d03.2017-07-13_07:10:00.nc')
fdomain04 = xr.open_dataset(path + 'met_em.d04.2017-07-13_07:10:00.nc')
# Extract latitude and longitude for each domain
LAT01 = fdomain01.XLAT_M.values  # Shape: (1, lat, lon)
LON01 = fdomain01.XLONG_M.values # Shape: (1, lat, lon)
LAT02 = fdomain02.XLAT_M.values  # Shape: (1, lat, lon)
LON02 = fdomain02.XLONG_M.values # Shape: (1, lat, lon)
LAT03 = fdomain03.XLAT_M.values  # Shape: (1, lat, lon)
LON03 = fdomain03.XLONG_M.values # Shape: (1, lat, lon)
LAT04 = fdomain04.XLAT_M.values  # Shape: (1, lat, lon)
LON04 = fdomain04.XLONG_M.values # Shape: (1, lat, lon)

path = '/pscratch/sd/h/heroplr/wrf_scratch/WRF_dm/test/Domain4_20170713_chem_aer1/'
casename = 'WRF'
file_str='wrfrst_d01_2017-07-13_04:00:00'
fin = xr.open_dataset(path+file_str)

nSO4_xro = xr.open_dataset("nSO4_interp_domain4_2017-07-13_domain4.nc")["nSO4_interp_domain4_2017-07-13"]
so4i = fin["so4ai"]
so4j = fin["so4aj"]
#seas = fin["seas"]

so4i_new = so4i.copy()
so4j_new = so4j.copy()
so4_i = 0.8e9 
so4_j = 0.2e9
MERRA_time = 1 # MERRA data every 10 min, total 144 time steps in 24 hours
for i in range(nSO4_xro.shape[3]):
    for j in range(nSO4_xro.shape[2]):
        # Sulfate
        so4i_new[0,:,j,i] = (nSO4_xro[MERRA_time,:,j,i].values) * so4_i 
        so4j_new[0,:,j,i] = (nSO4_xro[MERRA_time,:,j,i].values) * so4_j 

fin["so4ai"] = so4i_new
fin["so4aj"] = so4j_new
#fin["seas"] = seas_new

fin.to_netcdf(path+'wrfrst_d01_2017-07-13_04:00:00_so4.nc')

