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
file_str='wrfbdy_d01'
fin = xr.open_dataset(path+file_str)
refer02 = xr.open_dataset(path+'wrfinput_d01')
#refer02 = xr.open_dataset(path+'wrfinput_d02')
P = refer02["P"][0] + refer02["PB"][0]
time = fin['Time']
lats = refer02['XLAT'][0]
lons = refer02['XLONG'][0]
nSO4_xro = xr.open_dataset("nSO4_interp_domain4_2017-07-13.nc")["nSO4_interp"]
so4ai_BXS = fin["so4ai_BXS"]
so4ai_BXE = fin["so4ai_BXE"]
so4ai_BYS = fin["so4ai_BYS"]
so4ai_BYE = fin["so4ai_BYE"]
so4ai_BTXS = fin["so4ai_BTXS"]
so4ai_BTXE = fin["so4ai_BTXE"]
so4ai_BTYS = fin["so4ai_BTYS"]
so4ai_BTYE = fin["so4ai_BTYE"]
so4aj_BXS = fin["so4aj_BXS"]
so4aj_BXE = fin["so4aj_BXE"]
so4aj_BYS = fin["so4aj_BYS"]
so4aj_BYE = fin["so4aj_BYE"]
so4aj_BTXS = fin["so4aj_BTXS"]
so4aj_BTXE = fin["so4aj_BTXE"]
so4aj_BTYS = fin["so4aj_BTYS"]
so4aj_BTYE = fin["so4aj_BTYE"]
so4ai_BXS_new = so4ai_BXS.copy()
so4ai_BXE_new = so4ai_BXE.copy()
so4aj_BXS_new = so4aj_BXS.copy()
so4aj_BXE_new = so4aj_BXE.copy()
so4ai_BTXS_new = so4ai_BTXS.copy()
so4ai_BTXE_new = so4ai_BTXE.copy()
so4aj_BTXS_new = so4aj_BTXS.copy()
so4aj_BTXE_new = so4aj_BTXE.copy()
so4ai_BYS_new = so4ai_BYS.copy()
so4ai_BYE_new = so4ai_BYE.copy()
so4aj_BYS_new = so4aj_BYS.copy()
so4aj_BYE_new = so4aj_BYE.copy()
so4ai_BTYS_new = so4ai_BTYS.copy()
so4ai_BTYE_new = so4ai_BTYE.copy()
so4aj_BTYS_new = so4aj_BTYS.copy()
so4aj_BTYE_new = so4aj_BTYE.copy()

for i in range(549):

    # Sulfate
    so4ai_BYS_new[:,0,:,i] = (nSO4_xro[:,:,0,i].values) * 0.8e9 * 5.0
    so4ai_BYE_new[:,0,:,i] = (nSO4_xro[:,:,-1,i].values) * 0.8e9 * 5.0
    so4aj_BYS_new[:,0,:,i] = (nSO4_xro[:,:,0,i].values) * 0.2e9 * 5.0
    so4aj_BYE_new[:,0,:,i] = (nSO4_xro[:,:,-1,i].values) * 0.2e9 * 5.0


for i in range(529):

    # Sulfate
    so4ai_BXS_new[:,0,:,i] = (nSO4_xro[:,:,i,0].values) * 0.8e9 * 5.0
    so4ai_BXE_new[:,0,:,i] = (nSO4_xro[:,:,i,-1].values) * 0.8e9 * 5.0
    so4aj_BXS_new[:,0,:,i] = (nSO4_xro[:,:,i,0].values) * 0.2e9 * 5.0
    so4aj_BXE_new[:,0,:,i] = (nSO4_xro[:,:,i,-1].values) * 0.2e9 * 5.0


dt = 1./600. 
# Tendency X bounday
for t in range(len(time)-1):
    so4ai_BTXS_new[t,0,:,:] = (so4ai_BXS_new[t+1,0,:,:]- so4ai_BXS_new[t,0,:,:])*dt
    so4ai_BTXE_new[t,0,:,:] = (so4ai_BXE_new[t+1,0,:,:]- so4ai_BXE_new[t,0,:,:])*dt
    so4aj_BTXS_new[t,0,:,:] = (so4aj_BXS_new[t+1,0,:,:]- so4aj_BXS_new[t,0,:,:])*dt
    so4aj_BTXE_new[t,0,:,:] = (so4aj_BXE_new[t+1,0,:,:]- so4aj_BXE_new[t,0,:,:])*dt


# Tendency Y bounday
for t in range(len(time)-1):
    so4ai_BTYS_new[t,0,:,:] = (so4ai_BYS_new[t+1,0,:,:]- so4ai_BYS_new[t,0,:,:])*dt
    so4ai_BTYE_new[t,0,:,:] = (so4ai_BYE_new[t+1,0,:,:]- so4ai_BYE_new[t,0,:,:])*dt
    so4aj_BTYS_new[t,0,:,:] = (so4aj_BYS_new[t+1,0,:,:]- so4aj_BYS_new[t,0,:,:])*dt
    so4aj_BTYE_new[t,0,:,:] = (so4aj_BYE_new[t+1,0,:,:]- so4aj_BYE_new[t,0,:,:])*dt

# Tendency X bounday
so4ai_BTXS_new[-1,0,:,:] = so4ai_BTXS_new[-2,0,:,:]
so4ai_BTXE_new[-1,0,:,:] = so4ai_BTXE_new[-2,0,:,:]
so4aj_BTXS_new[-1,0,:,:] = so4aj_BTXS_new[-2,0,:,:]
so4aj_BTXE_new[-1,0,:,:] = so4aj_BTXE_new[-2,0,:,:]

# X bounday
for i in range(1,5):
    so4ai_BXS_new[:,i,:,:] = so4ai_BXS_new[:,0,:,:]
    so4ai_BXE_new[:,i,:,:] = so4ai_BXE_new[:,0,:,:]
    so4aj_BXS_new[:,i,:,:] = so4aj_BXS_new[:,0,:,:]
    so4aj_BXE_new[:,i,:,:] = so4aj_BXE_new[:,0,:,:]

# X bounday
for i in range(1,5):
    so4ai_BTXS_new[:,i,:,:] = so4ai_BTXS_new[:,0,:,:]
    so4ai_BTXE_new[:,i,:,:] = so4ai_BTXE_new[:,0,:,:]
    so4aj_BTXS_new[:,i,:,:] = so4aj_BTXS_new[:,0,:,:]
    so4aj_BTXE_new[:,i,:,:] = so4aj_BTXE_new[:,0,:,:]


# Tendency X bounday
so4ai_BTYS_new[-1,0,:,:] = so4ai_BTYS_new[-2,0,:,:]
so4ai_BTYE_new[-1,0,:,:] = so4ai_BTYE_new[-2,0,:,:]
so4aj_BTYS_new[-1,0,:,:] = so4aj_BTYS_new[-2,0,:,:]
so4aj_BTYE_new[-1,0,:,:] = so4aj_BTYE_new[-2,0,:,:]

# Y bounday
for i in range(1,5):
    so4ai_BYS_new[:,i,:,:] = so4ai_BYS_new[:,0,:,:]
    so4ai_BYE_new[:,i,:,:] = so4ai_BYE_new[:,0,:,:]
    so4aj_BYS_new[:,i,:,:] = so4aj_BYS_new[:,0,:,:]
    so4aj_BYE_new[:,i,:,:] = so4aj_BYE_new[:,0,:,:]

# Y bounday
for i in range(1,5):
    so4ai_BTYS_new[:,i,:,:] = so4ai_BTYS_new[:,0,:,:]
    so4ai_BTYE_new[:,i,:,:] = so4ai_BTYE_new[:,0,:,:]
    so4aj_BTYS_new[:,i,:,:] = so4aj_BTYS_new[:,0,:,:]
    so4aj_BTYE_new[:,i,:,:] = so4aj_BTYE_new[:,0,:,:]


fin["so4ai_BXS"] = so4ai_BXS_new
fin["so4ai_BXE"] = so4ai_BXE_new
fin["so4ai_BYS"] = so4ai_BYS_new
fin["so4ai_BYE"] = so4ai_BYE_new
fin["so4aj_BXS"] = so4aj_BXS_new
fin["so4aj_BXE"] = so4aj_BXE_new
fin["so4aj_BYS"] = so4aj_BYS_new
fin["so4aj_BYE"] = so4aj_BYE_new

fin["so4ai_BTXS"] = so4ai_BTXS_new
fin["so4ai_BTXE"] = so4ai_BTXE_new
fin["so4ai_BTYS"] = so4ai_BTYS_new
fin["so4ai_BTYE"] = so4ai_BTYE_new
fin["so4aj_BTXS"] = so4aj_BTXS_new
fin["so4aj_BTXE"] = so4aj_BTXE_new
fin["so4aj_BTYS"] = so4aj_BTYS_new
fin["so4aj_BTYE"] = so4aj_BTYE_new


nu0_BXS = fin["nu0_BXS"]
nu0_BXE = fin["nu0_BXE"]
nu0_BYS = fin["nu0_BYS"]
nu0_BYE = fin["nu0_BYE"]
nu0_BTXS = fin["nu0_BTXS"]
nu0_BTXE = fin["nu0_BTXE"]
nu0_BTYS = fin["nu0_BTYS"]
nu0_BTYE = fin["nu0_BTYE"]

ac0_BXS = fin["ac0_BXS"]
ac0_BXE = fin["ac0_BXE"]
ac0_BYS = fin["ac0_BYS"]
ac0_BYE = fin["ac0_BYE"]
ac0_BTXS = fin["ac0_BTXS"]
ac0_BTXE = fin["ac0_BTXE"]
ac0_BTYS = fin["ac0_BTYS"]
ac0_BTYE = fin["ac0_BTYE"]

nu0_BXS_new = nu0_BXS.copy()
nu0_BXE_new = nu0_BXE.copy()
nu0_BYS_new = nu0_BYS.copy()
nu0_BYE_new = nu0_BYE.copy()
nu0_BTXS_new = nu0_BTXS.copy()
nu0_BTXE_new = nu0_BTXE.copy()
nu0_BTYS_new = nu0_BTYS.copy()
nu0_BTYE_new = nu0_BTYE.copy()

ac0_BXS_new = ac0_BXS.copy()
ac0_BXE_new = ac0_BXE.copy()
ac0_BYS_new = ac0_BYS.copy()
ac0_BYE_new = ac0_BYE.copy()
ac0_BTXS_new = ac0_BTXS.copy()
ac0_BTXE_new = ac0_BTXE.copy()
ac0_BTYS_new = ac0_BTYS.copy()
ac0_BTYE_new = ac0_BTYE.copy()

#so4_i = 1.1845179250904482e-13
#so4_j = 2.246196954245592e-12

dt = 1.e-1
for i in range(len(time)):  
    nu0_BXS_new[i,:,:,:] = so4ai_BXS_new[i,:,:,:]*dt/so4_i 
    nu0_BXE_new[i,:,:,:] = so4ai_BXE_new[i,:,:,:]*dt/so4_i 
    nu0_BYS_new[i,:,:,:] = so4ai_BYS_new[i,:,:,:]*dt/so4_i 
    nu0_BYE_new[i,:,:,:] = so4ai_BYE_new[i,:,:,:]*dt/so4_i
    nu0_BTXS_new[i,:,:,:] = so4ai_BTXS_new[i,:,:,:]*dt/so4_i
    nu0_BTXE_new[i,:,:,:] = so4ai_BTXE_new[i,:,:,:]*dt/so4_i 
    nu0_BTYS_new[i,:,:,:] = so4ai_BTYS_new[i,:,:,:]*dt/so4_i
    nu0_BTYE_new[i,:,:,:] = so4ai_BTYE_new[i,:,:,:]*dt/so4_i

for i in range(len(time)):  
    ac0_BXS_new[i,:,:,:] = so4aj_BXS_new[i,:,:,:]*dt/so4_j 
    ac0_BXE_new[i,:,:,:] = so4aj_BXE_new[i,:,:,:]*dt/so4_j 
    ac0_BYS_new[i,:,:,:] = so4aj_BYS_new[i,:,:,:]*dt/so4_j 
    ac0_BYE_new[i,:,:,:] = so4aj_BYE_new[i,:,:,:]*dt/so4_j
    ac0_BTXS_new[i,:,:,:] = so4aj_BTXS_new[i,:,:,:]*dt/so4_j
    ac0_BTXE_new[i,:,:,:] = so4aj_BTXE_new[i,:,:,:]*dt/so4_j 
    ac0_BTYS_new[i,:,:,:] = so4aj_BTYS_new[i,:,:,:]*dt/so4_j
    ac0_BTYE_new[i,:,:,:] = so4aj_BTYE_new[i,:,:,:]*dt/so4_j    

fin.to_netcdf(path+'wrfbdy_d01_so4added.nc')


