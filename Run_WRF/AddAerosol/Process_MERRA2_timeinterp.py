import xarray as xr
import numpy as np
import pandas as pd

# Define interpolation function + save to file
def interp_and_save(var_stack, lon, lat, varname, domain_id):
    time_3hr = np.arange(0, 25, 3)  # [0, 3, ..., 24] → 9 points
    time_10min = np.arange(0, 24, 1/6)

    # Pad the time dimension if only 8 time steps are given
    if var_stack.shape[0] == 8:
        last_step = var_stack[-1:, :, :, :]  # shape: (1, lev, lat, lon)
        var_stack = np.concatenate([var_stack, last_step], axis=0)  # shape: (9, ...)

    lev_dim, lat_dim, lon_dim = var_stack.shape[1:]

    da = xr.DataArray(
        var_stack,
        dims=["time", "lev", "lat", "lon"],
        coords={
            "time": time_3hr,
            "lev": np.arange(lev_dim),
            "lat": np.arange(lat_dim),
            "lon": np.arange(lon_dim),
        },
        name=varname
    )

    # Interpolate
    interp_da = da.interp(time=time_10min)

    # Replace lat/lon coordinates with 2D arrays
    interp_da = interp_da.assign_coords({
        "lat": (("lat", "lon"), lat),
        "lon": (("lat", "lon"), lon)
    })

    # Set proper datetime
    interp_da["time"] = pd.date_range("2017-07-18", periods=len(time_10min), freq="10min")

    # Save to NetCDF
    filename = f"{varname}_domain{domain_id}.nc"
    interp_da.to_netcdf(filename)
    print(f"Saved: {filename}")


obspath_d04 = '/pscratch/sd/h/heroplr/wrf_scratch/MERRA2_Data/domain4/'
obsfile_t0 = 'regridded_merra2_d04_to_wrf_logv3_20170718_t0.nc'
obsfile_t1 = 'regridded_merra2_d04_to_wrf_logv3_20170718_t1.nc'
obsfile_t2 = 'regridded_merra2_d04_to_wrf_logv3_20170718_t2.nc'
obsfile_t3 = 'regridded_merra2_d04_to_wrf_logv3_20170718_t3.nc'
obsfile_t4 = 'regridded_merra2_d04_to_wrf_logv3_20170718_t4.nc'
obsfile_t5 = 'regridded_merra2_d04_to_wrf_logv3_20170718_t5.nc'
obsfile_t6 = 'regridded_merra2_d04_to_wrf_logv3_20170718_t6.nc'
obsfile_t7 = 'regridded_merra2_d04_to_wrf_logv3_20170718_t7.nc'
obsin_t40 = xr.open_mfdataset(obspath_d04+obsfile_t0)
obsin_t41 = xr.open_mfdataset(obspath_d04+obsfile_t1)
obsin_t42 = xr.open_mfdataset(obspath_d04+obsfile_t2)
obsin_t43 = xr.open_mfdataset(obspath_d04+obsfile_t3)
obsin_t44 = xr.open_mfdataset(obspath_d04+obsfile_t4)
obsin_t45 = xr.open_mfdataset(obspath_d04+obsfile_t5)
obsin_t46 = xr.open_mfdataset(obspath_d04+obsfile_t6)
obsin_t47 = xr.open_mfdataset(obspath_d04+obsfile_t7)

# === DOMAIN 4 ===

# Stack arrays
nBC_stack = np.stack([obsin_t40.nBC_interp, obsin_t41.nBC_interp, obsin_t42.nBC_interp, obsin_t43.nBC_interp,
                      obsin_t44.nBC_interp, obsin_t45.nBC_interp, obsin_t46.nBC_interp, obsin_t47.nBC_interp], axis=0)

nOC_stack = np.stack([obsin_t40.nOC_interp, obsin_t41.nOC_interp, obsin_t42.nOC_interp, obsin_t43.nOC_interp,
                      obsin_t44.nOC_interp, obsin_t45.nOC_interp, obsin_t46.nOC_interp, obsin_t47.nOC_interp], axis=0)

nSO2_stack = np.stack([obsin_t40.nSO2_interp, obsin_t41.nSO2_interp, obsin_t42.nSO2_interp, obsin_t43.nSO2_interp,
                       obsin_t44.nSO2_interp, obsin_t45.nSO2_interp, obsin_t46.nSO2_interp, obsin_t47.nSO2_interp], axis=0)

nSO4_stack = np.stack([obsin_t40.nSO4_interp, obsin_t41.nSO4_interp, obsin_t42.nSO4_interp, obsin_t43.nSO4_interp,
                       obsin_t44.nSO4_interp, obsin_t45.nSO4_interp, obsin_t46.nSO4_interp, obsin_t47.nSO4_interp], axis=0)

# Get coordinate grids
lon_d04 = obsin_t40.lon.values
lat_d04 = obsin_t40.lat.values

# Interpolate + Save
interp_and_save(nBC_stack, lon_d04, lat_d04, "nBC_interp_domain4_2017-07-18", domain_id=4)
interp_and_save(nOC_stack, lon_d04, lat_d04, "nOC_interp_domain4_2017-07-18", domain_id=4)
interp_and_save(nSO2_stack, lon_d04, lat_d04, "nSO2_interp_domain4_2017-07-18", domain_id=4)
interp_and_save(nSO4_stack, lon_d04, lat_d04, "nSO4_interp_domain4_2017-07-18", domain_id=4)
