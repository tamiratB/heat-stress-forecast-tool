#!/usr/bin/env python3
# ============================================================================
# calc_spatial_WBGT_forecast.py
#
# Driver that computes gridded Wet Bulb Globe Temperature (WBGT) from the
# 6-hourly ECMWF forecast fields (ecmwf_forecasts4wbgt_YYYYMMDD.nc).
#
# Uses the Liljegren et al. (2008) model in wbgt_functions.py, solving twice
# per step: warm extreme (t2max, rh_tmax) -> daytime WBGT, and cool extreme
# (t2min, rh_tmin) -> nighttime WBGT. Writes WBGT_forecast_output_YYYY-MM-DD.nc
# with WBGT_tmax, WBGT_tmin and their Tw/Tg components (degC).
#
# Developed by: @ICPAC
# ============================================================================

import xarray as xr
import numpy as np
from datetime import datetime, timezone
import wbgt_functions as wbgt
import os

# solver iteration size, MAX_ITER, should be
# adjust based on numerical convergence requirement that highly depend on the physical process
# 200 or less is enough in most cases

MAX_ITER = 100
SAVE_DIAGNOSTICS = None  # place holder for future development
DIAG_LOCATIONS = None    # same


def compute_gridded_wbgt(ds, nc_out, max_iter, save_diagnostics=True,
                         diag_locations=None):

    sp = ds.sp / 100
    ssrd = ds.ssrd
    Uk = 4.87 / np.log(67.8 * 10 - 5.42)
    U = ds.U * Uk
    d2m = ds.d2m

    S_instant, fdir_instant, theta_instant = wbgt.compute_instantaneous_solar(
            ssrd,
            debug=False
        )

    def _solve_wbgt(t2m, rh, label):

        td = xr.where(d2m > t2m, t2m, d2m)
        Tw, Tg, WBGT = wbgt.compute_WBGT(
                td,
                t2m,
                rh,
                sp,
                U,
                S_instant, fdir_instant, theta_instant,
                debug=False,
                max_iter=max_iter,
                test_year=None,
                save_diagnostics=save_diagnostics,
                diag_locations=diag_locations,
            )

        return Tw - 273.15, Tg - 273.15, WBGT - 273.15

    print("\nComputing daytime heat stress...")
    Tw_max, Tg_max, WBGT_max = _solve_wbgt(ds.t2max, ds.rh_tmax, "tmax")
    print("\nComputing nighttime heat stress...")
    Tw_min, Tg_min, WBGT_min = _solve_wbgt(ds.t2min, ds.rh_tmin, "tmin")

    out_ds = xr.Dataset(
        {
            "WBGT_tmax": WBGT_max.assign_attrs({
                "long_name": "Daytime WBGT Forecast",
                "units": "degC",
                "standard_name": "wet_bulb_globe_temperature",
                "description": "daytime/warm-extreme heat stress"
            }),
            "WBGT_tmin": WBGT_min.assign_attrs({
                "long_name": "Nighttime WBGT Forecast",
                "units": "degC",
                "standard_name": "wet_bulb_globe_temperature",
                "description": "nighttime/cool-extreme heat stress"
            }),
            "Tw_tmax": Tw_max.assign_attrs({
                "long_name": "Natural Wet Bulb Temperature (daytime, from t2max)",
                "units": "degC",
                "standard_name": "wet_bulb_temperature",
                "description": "Natural wet bulb temperature used in the WBGT_tmax calculation"
            }),
            "Tg_tmax": Tg_max.assign_attrs({
                "long_name": "Globe Temperature (daytime, from t2max)",
                "units": "degC",
                "standard_name": "black_globe_temperature",
                "description": "Black globe temperature used in the WBGT_tmax calculation"
            }),
            "Tw_tmin": Tw_min.assign_attrs({
                "long_name": "Natural Wet Bulb Temperature (nighttime, from t2min)",
                "units": "degC",
                "standard_name": "wet_bulb_temperature",
                "description": "Natural wet bulb temperature used in the WBGT_tmin calculation"
            }),
            "Tg_tmin": Tg_min.assign_attrs({
                "long_name": "Globe Temperature (nighttime, from t2min)",
                "units": "degC",
                "standard_name": "black_globe_temperature",
                "description": "Black globe temperature used in the WBGT_tmin calculation"
            }),
        },
        coords=ds.coords,
        attrs={
            "title": "Wet Bulb Globe Temperature (WBGT) Forecast - day (t2max) and night (t2min)",
            "institution": "IGAD Climate Prediction and Application Center - ICPAC",
            "source": "Computed using Liljegren WBGT model",
            "history": f"Created on {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S UTC')}",
            "references": "Liljegren et al. (2008), DOI: 10.1080/15459620802310770",
            "Conventions": "CF-1.8"
        }
    )

    for coord in out_ds.coords:
        out_ds[coord].attrs = ds[coord].attrs

    out_ds.to_netcdf(
        nc_out,
        engine="netcdf4",
        encoding=wbgt.nc_write_encoding(
            out_ds, compress=True, least_significant_digit=2
        )
    )
    print(f"Saved: {nc_out}")


if __name__ == "__main__":

    print(f"Processing forecasts initialized on {datetime.now().strftime('%Y-%m-%d')}...")
    nc_in = f"ecmwf_forecasts4wbgt_{datetime.now().strftime('%Y%m%d')}.nc"

    ds = xr.open_dataset(
            nc_in,
        )

    if "valid_time" in ds.coords or "valid_time" in ds.dims:
        ds = ds.rename({"valid_time": "time"})
    if "latitude" in ds.coords or "longitude" in ds.dims:
        ds = ds.rename({"latitude": "lat"})
        ds = ds.rename({"longitude": "lon"})

    os.makedirs('forecasts', exist_ok=True)
    nc_out = f"forecasts/WBGT_forecast_output_{datetime.now().strftime("%Y-%m-%d")}.nc"

    compute_gridded_wbgt(ds,
                         nc_out,
                         MAX_ITER,
                         save_diagnostics=SAVE_DIAGNOSTICS,
                         diag_locations=DIAG_LOCATIONS
                         )

    print(f"\nForecast computation successfully completed! {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}")

