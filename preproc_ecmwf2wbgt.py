#!/usr/bin/env python3
"""
=====================================================================
Prepare 6-hourly WBGT input fields (day/night heat stress) from ECMWF
open-data forecasts.
Developed by: @ICPAC
=====================================================================
"""

import glob
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime

# ----------------------------------------------------------
# User settings
# ----------------------------------------------------------
# domain extent
lat_min, lat_max = -12.5, 23.5
lon_min, lon_max = 21, 52

# ---------------------------------------------------------
today = datetime.now().strftime("%Y%m%d")

print(f"Processing input forecast data for {today}...")

step_hours = 6


def preprocess(ds):

    vt = pd.to_datetime(ds.valid_time.values)
    ds = ds.expand_dims(dim={"time": [vt]})
    ds = ds.drop_vars(["valid_time", "step"], errors="ignore")

    return ds

def open_level(file_list, filter_keys):
    return xr.open_mfdataset(
        file_list,
        engine="cfgrib",
        combine="nested",
        concat_dim="time",
        preprocess=preprocess,
        parallel=True,
        decode_timedelta=True,
        coords="minimal",
        compat="override",
        backend_kwargs={
            "filter_by_keys": filter_keys,
            "indexpath": ""
        }
    )


def prepare_ecmwf(ds):
    ds = ds.sortby("time")

    # Subset domain early to keep memory low
    ds = ds.sel(
        latitude=slice(lat_max, lat_min),
        longitude=slice(lon_min, lon_max)
    )

    return ds


def calc_rh(T, Td):
    T = T - 273.15
    Td = Td - 273.15

    es = 6.112 * np.exp((17.67 * T) / (T + 243.5))
    e = 6.112 * np.exp((17.67 * Td) / (Td + 243.5))

    return np.clip(100 * (e / es), 0, 100)


def step_of(fname):
    # ECMWF_sfc_YYYYMMDDHH_SSS.grib2 -> SSS
    return int(fname.rsplit("_", 1)[-1].split(".")[0])


data_dir = f"ecmwf_forecasts_{today}"
files = sorted(glob.glob(f"{data_dir}/ECMWF_sfc_*.grib2"))
files_wk12 = [f for f in files if 3 <= step_of(f) <= 144]   # 3-hourly segment
files_wk34 = [f for f in files if step_of(f) >= 150]        # 6-hourly segment
print(f"Opening {len(files)} GRIB2 files "
      f"({len(files_wk12)} 3-hourly + {len(files_wk34)} 6-hourly steps)...")

print("Reading surface variables...")
ds_surface = prepare_ecmwf(open_level(files, {"typeOfLevel": "surface"}))

print("Reading 10m variables...")
ds_10m = prepare_ecmwf(open_level(
    files, {"typeOfLevel": "heightAboveGround", "level": 10}))

print("Reading 2m instantaneous variables...")
ds_2m = prepare_ecmwf(open_level(
    files, {"typeOfLevel": "heightAboveGround", "level": 2,
            "stepType": "instant"}))

print("Reading 2m max/min temperature (3-h windows, 0-144 h)...")
ds_max3 = prepare_ecmwf(open_level(
    files_wk12, {"typeOfLevel": "heightAboveGround", "stepType": "max"}))
ds_min3 = prepare_ecmwf(open_level(
    files_wk12, {"typeOfLevel": "heightAboveGround", "stepType": "min"}))

print("Reading 2m max/min temperature (6-h windows, 150-360 h)...")
ds_max6 = prepare_ecmwf(open_level(
    files_wk34, {"typeOfLevel": "heightAboveGround", "stepType": "max"}))
ds_min6 = prepare_ecmwf(open_level(
    files_wk34, {"typeOfLevel": "heightAboveGround", "stepType": "min"}))

t0 = ds_surface.time.values[0]  # initialization (step 0)
hours_from_init = ((ds_surface.time.values - t0)
                   / np.timedelta64(1, "h")).astype(int)
time6_full = ds_surface.time.values[hours_from_init % step_hours == 0]

# (step 0 has no preceding interval, so it is dropped)
time6 = time6_full[1:]
h6 = ((time6 - t0) / np.timedelta64(1, "h")).astype(int)
time6_wk12 = time6[h6 <= 144]   # covered by paired 3-h windows
time6_wk34 = time6[h6 > 144]    # covered by native 6-h windows
print(f"Output: {len(time6)} 6-hourly steps "
      f"({len(time6_wk12)} from 3-h pairs + {len(time6_wk34)} native 6-h)")


def clean(da):
    return da.drop_vars(["surface", "heightAboveGround"], errors="ignore")


def to_6h(da):
    return clean(da).chunk({"time": -1}).interp(time=time6, method="linear")

mx3 = clean(ds_max3["mx2t3"]).chunk({"time": -1})
mn3 = clean(ds_min3["mn2t3"]).chunk({"time": -1})

tmax_wk12 = mx3.rolling(time=2, min_periods=1).max().sel(time=time6_wk12)
tmin_wk12 = mn3.rolling(time=2, min_periods=1).min().sel(time=time6_wk12)

tmax_wk34 = clean(ds_max6["mx2t6"]).sel(time=time6_wk34)
tmin_wk34 = clean(ds_min6["mn2t6"]).sel(time=time6_wk34)

tmax = xr.concat([tmax_wk12, tmax_wk34], dim="time").rename("t2max")
tmin = xr.concat([tmin_wk12, tmin_wk34], dim="time").rename("t2min")

pressure = to_6h(ds_surface["sp"]) # in Pa / 100   # Pa -> hPa
u10 = to_6h(ds_10m["u10"])
v10 = to_6h(ds_10m["v10"])
d2m = to_6h(ds_2m["d2m"])

ssrd_accum = clean(ds_surface["ssrd"]).sel(time=time6_full)
dt_seconds = (ssrd_accum["time"].diff("time") / np.timedelta64(1, "s")).astype("float64")
if (dt_seconds != step_hours * 3600).any():
    print("WARNING: gaps in the 6-hourly axis detected - "
          "check for missing downloads!")

solar = ssrd_accum.diff("time").clip(min=0)  # J m-2 per 6-h interval
wind = np.sqrt(u10**2 + v10**2)
rh_tmax = calc_rh(tmax, d2m)
rh_tmin = calc_rh(tmin, d2m)

print("Merging variables...")

ds = xr.Dataset(
    {
        "t2max": tmax.assign_attrs({
                "long_name": "Maximum 2m air temperature over the previous 6 hours",
                "standard_name": "air_temperature",
                "units": "K",
                "cell_methods": "time: maximum (interval: 6 hours)",
                "description": "maximum near-surface temperature"
            }),
        "t2min": tmin.assign_attrs({
                "long_name": "Minimum 2m air temperature over the previous 6 hours",
                "standard_name": "air_temperature",
                "units": "K",
                "cell_methods": "time: minimum (interval: 6 hours)",
                "description": "minimum near-surface temperature"
            }),
        "d2m": d2m.assign_attrs({
                "long_name": "2m dew point temperature",
                "standard_name": "dew_point_temperature",
                "units": "K",
                "description": "Instantaneous value at the end of each 6-h interval"
            }),
        "rh_tmax": rh_tmax.assign_attrs({
                "long_name": "Relative humidity at the 6-h maximum temperature",
                "standard_name": "relative_humidity",
                "units": "%",
                "description": "Relative humidity during the hottest conditions of the interval"
            }),
        "rh_tmin": rh_tmin.assign_attrs({
                "long_name": "Relative humidity at the 6-h minimum temperature",
                "standard_name": "relative_humidity",
                "units": "%",
                "description": ("Relative humidity during the coolest conditions of the interval")
            }),
        "U": wind.assign_attrs({
                "long_name": "10m wind speed",
                "standard_name": "wind_speed",
                "units": "m s-1",
                "description": ("Resultant 10m wind speed magnitude")
            }),
        "sp": pressure.assign_attrs({
                "long_name": "Surface pressure",
                "standard_name": "surface_air_pressure",
                "units": "Pa"
            }),
        "ssrd": solar.assign_attrs({
                "long_name": "Surface solar radiation downwards",
                "standard_name": "surface_downwelling_shortwave_flux_in_air",
                "units": "J m-2",
                "description": "Downward shortwave radiation accumulated over the 6 hours"
            }),
    },
    attrs={
            "title": "6-hourly meteorological variables for day/night heat stress analysis",
            "institution": "IGAD Climate Prediction and Application Center - ICPAC",
            "source": "ECMWF real-time forecast data from the IFS (open-data, 0.25 deg)",
            "history": f"Created on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} UTC",
            "references": "https://www.ecmwf.int/en/forecasts/datasets/open-data; https://doi.org/10.21957/open-data",
            "Conventions": "CF-1.8",
            "time_note": "Times are in UTC",
        }
)

drop_vars = [v for v in ["surface", "heightAboveGround", "valid_time"] if v in ds.variables]
ds = ds.drop_vars(drop_vars, errors="ignore")

out_file = f"ecmwf_forecasts4wbgt_{today}.nc"

print("Writing output file...")
encoding = {v: {"zlib": True, "complevel": 4} for v in ds.data_vars}
ds.to_netcdf(out_file, encoding=encoding)

print(f"\n Completed: {out_file}")
