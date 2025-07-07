import os
import pandas as pd
import cdsapi

# === Your date list ===
flight_dates = pd.to_datetime([
    # 2017 Summer
    "2017-06-21", "2017-06-24", "2017-06-25", "2017-06-26", "2017-06-29",
    "2017-06-30", "2017-07-01", "2017-07-03", "2017-07-04", "2017-07-06",
    "2017-07-08", "2017-07-09", "2017-07-11", "2017-07-12", "2017-07-13",
    "2017-07-15", "2017-07-16", "2017-07-18", "2017-07-19",
    # 2018 Winter
    "2018-01-19", "2018-01-23", "2018-01-24", "2018-01-26", "2018-01-27",
    "2018-01-28", "2018-01-30", "2018-01-31", "2018-02-01", "2018-02-03",
    "2018-02-08", "2018-02-09", "2018-02-11", "2018-02-12", "2018-02-13",
    "2018-02-16", "2018-02-17", "2018-02-18", "2018-02-19"
])

# === Output directory ===
DATADIR = "./era5_soil_downloads/"
os.makedirs(DATADIR, exist_ok=True)

# === Initialize the cdsapi client ===
c = cdsapi.Client()

# === Optional: kwargs (remove if you don't use them)
# cdsapi_kwargs = {}

# === Loop over dates and download ===
for date in flight_dates:
    year = date.strftime("%Y")
    month = date.strftime("%m")
    day = date.strftime("%d")

    download_file = f"{DATADIR}ERA5_wrf_soil-{year}_{month}_{day}.grib"

    if not os.path.isfile(download_file):
        print(f"\nDownloading data for {year}-{month}-{day}...")

        try:
            c.retrieve(
                'reanalysis-era5-single-levels',
                {
                    'product_type': 'reanalysis',
                    'variable': [
                        "soil_temperature_level_1",
                        "soil_temperature_level_2",
                        "soil_temperature_level_3",
                        "soil_temperature_level_4",
                        "soil_type",
                        "volumetric_soil_water_layer_1",
                        "volumetric_soil_water_layer_2",
                        "volumetric_soil_water_layer_3",
                        "volumetric_soil_water_layer_4",
                    ],
                    'year': year,
                    'month': [month],
                    'day': [day],
                    'time': [
                        "00:00", "01:00", "02:00", "03:00", "04:00", "05:00",
                        "06:00", "07:00", "08:00", "09:00", "10:00", "11:00",
                        "12:00", "13:00", "14:00", "15:00", "16:00", "17:00",
                        "18:00", "19:00", "20:00", "21:00", "22:00", "23:00"
                    ],
                    'data_format': 'grib'
                },
            ).download(download_file)

            print(f"✔ Data saved as {download_file}")

        except Exception as e:
            print(f"✘ Failed to download data for {year}-{month}-{day}: {str(e)}")

    else:
        print(f"✔ File for {year}-{month}-{day} already exists, skipping download.")

