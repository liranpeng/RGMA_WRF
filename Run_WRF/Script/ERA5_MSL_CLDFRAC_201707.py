import os
import glob
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt

# === Paths ===
era5_dir = "/pscratch/sd/h/heroplr/jupyternotebook/era5_july2017"
wrf_dir = "/pscratch/sd/h/heroplr/wrf_scratch/MERRA2_Data/wrf_domain"

# === Load WRF files and store lat/lon for each domain ===
wrf_files = sorted(glob.glob(os.path.join(wrf_dir, "wrfout_d0?_2016-07-01_00:17:00")))
wrf_bounds = {}
for f in wrf_files:
    domain = os.path.basename(f).split('_')[1]
    ds = xr.open_dataset(f)
    lat = ds["XLAT"].isel(Time=0)
    lon = ds["XLONG"].isel(Time=0)
    wrf_bounds[domain] = {"lat": lat, "lon": lon}
    ds.close()

def wrap_lon(lon):
    return lon % 360

# WRF Domain 1 bounds
lat_d01 = wrf_bounds["d01"]["lat"].values
lon_d01 = wrf_bounds["d01"]["lon"].values
lat_min = lat_d01.min()
lat_max = lat_d01.max()
lon_min = lon_d01.min()
lon_max = lon_d01.max()
lon_min_wrapped = wrap_lon(lon_min)
lon_max_wrapped = wrap_lon(lon_max)

print(f"Domain 1 lat: {lat_min}–{lat_max}")
print(f"Domain 1 lon: {lon_min}–{lon_max} (wrapped: {lon_min_wrapped}–{lon_max_wrapped})")

# === ERA5 files ===
era5_files = sorted(glob.glob(os.path.join(era5_dir, "*.nc")))

domain_color = {
    "d01": "blue",
    "d02": "green",
    "d03": "orange",
    "d04": "black"
}

for era5_file in era5_files:
    print(f"Processing {era5_file}")
    ds = xr.open_dataset(era5_file)

    # Extract variables
    msl = ds["msl"] / 100.0
    u10 = ds["u10"]
    v10 = ds["v10"]
    lcc = ds["lcc"]  # Low cloud cover (%)

    # Confirm longitude convention
    is_era5_wrapped = ds["longitude"].values.max() > 180

    # Subset
    if is_era5_wrapped:
        if lon_min_wrapped < lon_max_wrapped:
            ds_subset = ds.sel(
                latitude=slice(lat_max, lat_min),
                longitude=slice(lon_min_wrapped, lon_max_wrapped)
            )
        else:
            ds_east = ds.sel(
                latitude=slice(lat_max, lat_min),
                longitude=slice(lon_min_wrapped, 360)
            )
            ds_west = ds.sel(
                latitude=slice(lat_max, lat_min),
                longitude=slice(0, lon_max_wrapped)
            )
            ds_subset = xr.concat([ds_east, ds_west], dim="longitude")
    else:
        ds_subset = ds.sel(
            latitude=slice(lat_max, lat_min),
            longitude=slice(lon_min, lon_max)
        )

    # Extract subset arrays
    msl_sub = ds_subset["msl"] / 100.0
    u10_sub = ds_subset["u10"]
    v10_sub = ds_subset["v10"]
    lcc_sub = ds_subset["lcc"]

    lats_sub = ds_subset["latitude"].values
    lons_sub = ds_subset["longitude"].values
    lon2d, lat2d = np.meshgrid(lons_sub, lats_sub)

    for t_idx in range(msl_sub.shape[0]):
        time = str(msl_sub["valid_time"].values[t_idx])

        fig, ax = plt.subplots(figsize=(10, 8))

        # === Low cloud fraction shading ===
        lcc_data = lcc_sub.isel(valid_time=t_idx).values
        # Mask zeros
        lcc_masked = np.ma.masked_where(lcc_data == 0, lcc_data)
        cf_lcc = ax.pcolormesh(
            lon2d,
            lat2d,
            lcc_masked,
            cmap="Greys",
            shading="auto",
            alpha=0.5
        )
        cbar = plt.colorbar(cf_lcc, ax=ax, orientation="vertical", label="Low Cloud Fraction (%)")

        # === MSL contours ===
        cs = ax.contour(
            lon2d,
            lat2d,
            msl_sub.isel(valid_time=t_idx),
            levels=np.arange(990, 1030, 2),
            colors="k",
            linewidths=0.5
        )
        ax.clabel(cs, inline=True, fontsize=8)

        # === Wind vectors ===
        skip = (slice(None, None, 5), slice(None, None, 5))
        ax.quiver(
            lon2d[skip],
            lat2d[skip],
            u10_sub.isel(valid_time=t_idx).values[skip],
            v10_sub.isel(valid_time=t_idx).values[skip],
            scale=500
        )

        # === WRF domain boundaries ===
        for dom, coords in wrf_bounds.items():
            lat = coords["lat"]
            lon = coords["lon"].copy()
            if is_era5_wrapped:
                lon = wrap_lon(lon)

            color = domain_color.get(dom, "gray")
            label = f"{dom} Edge" if dom == "d01" else None
            ax.plot(lon.isel(south_north=0), lat.isel(south_north=0), color=color, linewidth=1, label=label)
            ax.plot(lon.isel(south_north=-1), lat.isel(south_north=-1), color=color, linewidth=1)
            ax.plot(lon.isel(west_east=0), lat.isel(west_east=0), color=color, linewidth=1)
            ax.plot(lon.isel(west_east=-1), lat.isel(west_east=-1), color=color, linewidth=1)

        ax.set_title(f"{os.path.basename(era5_file)}\nTime index {t_idx} ({time})")
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")

        # Domain extent
        if is_era5_wrapped:
            ax.set_xlim([lon_min_wrapped, lon_max_wrapped])
        else:
            ax.set_xlim([lon_min, lon_max])
        ax.set_ylim([lat_min, lat_max])

        ax.legend(loc="lower left", fontsize="small", frameon=True)
        ax.grid(True)
        plt.tight_layout()

        # Save figure
        outname = f"msl_lcc_wind_{os.path.basename(era5_file).replace('.nc','')}_t{t_idx:02d}.png"
        plt.savefig(outname, dpi=150)
        plt.close()

    ds.close()

print("All plots completed.")

