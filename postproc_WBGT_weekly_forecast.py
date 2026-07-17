#!/usr/bin/env python3
""""
================================================================================
Weekly WBGT heat-stress summary maps from the 6-hourly WBGT forecast.

Input : WBGT_forecast_output_YYYY-MM-DD.nc  (from calc_spatial_WBGT_forecast.py)
        6-hourly, UTC, variables WBGT_tmax (daytime) / WBGT_tmin (nighttime).

Produces, per forecast week:
  A. Average of the MAXIMA  (typical peak heat stress)
  B. MAXIMUM WBGT                 (worst-case peak in the week)
  C. FREQUENCY WBGT > threshold   (% of 6-hourly intervals exceeding)

Also produces an operational heat-stress CATEGORY map per week (the weekly
peak WBGT classified into standard risk bands with work/rest guidance) -
the recommended layout for briefing non-specialist decision makers.

Times are handled in East Africa Time (UTC+3) so that "daily" maxima and
the calendar weeks line up with the local day.
Developed by: @ICPAC
===============================================================================
"""

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from datetime import datetime, timedelta
import geopandas as gpd
import regionmask
import cmaps
import os

# -------------
# User settings
# --------------------------------------------------------------------------------
# load shapefiles
ADMIN0 = gpd.read_file("/data/shapefiles/gha/gha_admin0.shp").to_crs("EPSG:4326")
ADMIN1 = gpd.read_file("/data/shapefiles/gha/gha_admin1.shp").to_crs("EPSG:4326")

# exceedance threshold for freq. metric C  [degC]
THRESHOLD = 29.0

# daytime heat stress (WBGT_tmax) or nighttime stress (WBGT_tmin)
VAR = "WBGT_tmin"

# timeshift
TZ_SHIFT_HOURS = 3          # UTC -> East Africa Time

# domain extent
EXTENT = [21, 52, -12.5, 23.5]

# operational heat-stress categories (WBGT, degC). Bands follow the widely
# used occupational / sport WBGT flag system (ISO 7243-style reference
# limits)
CAT_BOUNDS = [0, 25, 28, 30, 32, 60]

# -------------------------------------------------------------------------------
today = datetime.now()
print(f"\nPost-processing forecasts for WBGT_forecast_output_{today.strftime('%Y-%m-%d')}.nc...")

# input file (defaults to today's run; override as needed)
NC_FILE = f"WBGT_forecast_output_{today.strftime('%Y-%m-%d')}.nc"

# forecast weeks, as inclusive EAT calendar-date ranges
week1_start = today
week1_end = today + timedelta(days=7)

week2_start = today + timedelta(days=8)
week2_end = today + timedelta(days=14)

WEEKS = {
    f"Week 1 ({week1_start:%d %b}-{week1_end:%d %b})": (
            week1_start.strftime("%Y-%m-%d"),
            week1_end.strftime("%Y-%m-%d"),
    ),
    f"Week 2 ({week2_start:%d %b}-{week2_end:%d %b})": (
            week2_start.strftime("%Y-%m-%d"),
            week2_end.strftime("%Y-%m-%d"),
    ),
}

# control number of data points in day to avoid next day cold bias
MIN_STEPS_PER_DAY = 3

proj = ccrs.PlateCarree()

CAT_LABELS = [
    "Low (<25): normal activity",
    "Moderate (25-28): stay hydrated",
    "High (28-30): limit exertion",
    "Very high (30-32): curtail outdoor work",
    "Extreme (>=32): avoid outdoor work",
]
CAT_COLORS = ["#4daf4a", "#ffff33", "#ff7f00", "#e41a1c", "#7f0000"]

HEAT_CMAP = cmaps.gui_default(np.linspace(0, 1, 22))
FREQ_CMAP = cmaps.prcp_1(np.linspace(0, 1, 17))

# remove first and last colors
FREQ_CMAP = ListedColormap(FREQ_CMAP[:-3])
HEAT_CMAP = ListedColormap(HEAT_CMAP[2:])


def add_map_features(ax):
    ax.set_extent(EXTENT, crs=proj)
    ax.add_geometries(ADMIN1['geometry'],
                      crs=ccrs.PlateCarree(),
                      facecolor='none',
                      edgecolor='black',
                      linewidth=1.5
                      )
    ax.add_feature(cfeature.LAKES,
                   linewidth=0.3,
                   edgecolor="black",
                   facecolor="none"
                   )
    gl = ax.gridlines(draw_labels=True,
                      linewidth=0.3,
                      linestyle="--",
                      alpha=0.4
                      )
    gl.top_labels = False
    gl.right_labels = False
    gl.xlocator = mticker.FixedLocator(np.arange(25, 56, 10))
    gl.ylocator = mticker.FixedLocator(np.arange(-10, 21, 10))


def _build_region_path():
    geom = ADMIN0.union_all()
    polys = list(geom.geoms) if geom.geom_type == "MultiPolygon" else [geom]
    verts, codes = [], []
    for poly in polys:
        for ring in [poly.exterior, *poly.interiors]:
            xy = np.asarray(ring.coords)
            if len(xy) < 3:
                continue
            verts.append(xy)
            codes.append(np.concatenate((
                [MplPath.MOVETO],
                np.full(len(xy) - 2, MplPath.LINETO),
                [MplPath.CLOSEPOLY],
            )))
    return MplPath(np.concatenate(verts), np.concatenate(codes))


_REGION_PATH = _build_region_path()


def clip_to_region(ax, artist):
    patch = PathPatch(_REGION_PATH,
                      transform=ax.transData,
                      facecolor="none",
                      edgecolor="none"
                      )
    ax.add_patch(patch)
    artist.set_clip_path(patch)


def weekly_metrics(da_week):
    daily_max = da_week.resample(time="1D").max()
    steps_per_day = (xr.ones_like(da_week.isel(lat=0, lon=0, drop=True)).resample(time="1D").sum())
    complete = steps_per_day >= MIN_STEPS_PER_DAY
    avg_daily_max = daily_max.where(complete, drop=True).mean("time")

    week_max = da_week.max("time")
    freq = 100.0 * (da_week > THRESHOLD).mean("time")

    return avg_daily_max, week_max, freq


ds = xr.open_dataset(NC_FILE)
ds["time"] = pd.to_datetime(ds.time.values) + pd.Timedelta(hours=TZ_SHIFT_HOURS)
ghaMask = regionmask.mask_geopandas(ADMIN0,
                                    ds.lon,
                                    ds.lat
                                    )
ds = ds.where(~ghaMask.isnull())
da = ds[VAR]

week_data = {name: da.sel(time=slice(start, end))
             for name, (start, end) in WEEKS.items()}

metrics = {name: weekly_metrics(sub) for name, sub in week_data.items()}

A_levels = np.linspace(23, 33, 11)
B_levels = np.linspace(23, 33, 11)
freq_levels = np.arange(0, 101, 10)

horizon = {
    "WBGT_tmax": "Daytime",
    "WBGT_tmin": "Nighttime"
    }

COLS = [
    (f"{horizon[VAR]} (weekly mean)", HEAT_CMAP, A_levels, "WBGT (°C)"),
    (f"{horizon[VAR]} (maximum)", HEAT_CMAP, B_levels, "WBGT (°C)"),
    (f"Frequency WBGT > {THRESHOLD:.0f}°C", FREQ_CMAP, freq_levels, "% of time"),
]

os.makedirs('plots', exist_ok=True)

nrows, ncols = len(WEEKS), 3
fig = plt.figure(figsize=(5.2 * ncols, 4.2 * nrows))
ims = [None, None, None]

for r, (wname, m) in enumerate(metrics.items()):
    for c in range(3):
        ax = plt.subplot(nrows, ncols, r * ncols + c + 1, projection=proj)
        title, cmap, levels, _ = COLS[c]
        ims[c] = ax.contourf(da.lon,
                             da.lat,
                             m[c],
                             levels=levels,
                             cmap=cmap,
                             extend="max" if c == 2 else "both",
                             transform=proj
                             )
        add_map_features(ax)
        if r == 0:
            ax.set_title(title, fontsize=13, fontweight="bold")
        if c == 0:
            ax.text(-0.18, 0.5, wname, transform=ax.transAxes,
                    rotation=90, va="center", ha="center",
                    fontsize=13, fontweight="bold")

for c in range(3):
    left = 0.09 + c * 0.315
    cax = fig.add_axes([left, 0.055, 0.22, 0.014])
    cb = plt.colorbar(ims[c],
                      cax=cax,
                      orientation="horizontal"
                      )
    cb.set_label(COLS[c][3], fontsize=11)

fig.suptitle("Weekly WBGT Heat Stress Summary",
             fontsize=17,
             fontweight="bold",
             y=0.99
             )
plt.subplots_adjust(left=0.06,
                    right=0.98,
                    top=0.93,
                    bottom=0.10,
                    wspace=0.12,
                    hspace=0.15
                    )

out1 = f"plots/{horizon[VAR]}_weekly_forecast_summary_{today.strftime('%Y-%m-%d')}.png"
plt.savefig(out1,
            dpi=200,
            bbox_inches="tight"
            )
print(f"Saved: {out1}")

fig2 = plt.figure(figsize=(6.2 * len(WEEKS), 6.4))
for i, (wname, m) in enumerate(metrics.items()):
    ax = plt.subplot(1, len(WEEKS), i + 1, projection=proj)
    pm = ax.contourf(da.lon,
                     da.lat,
                     m[1],          # m[1] = weekly max WBGT
                     levels=CAT_BOUNDS,
                     colors=CAT_COLORS,
                     extend="neither",
                     transform=proj
                     )
    clip_to_region(ax, pm)
    add_map_features(ax)
    ax.set_title(wname, fontsize=14, fontweight="bold")

# shared categorical legend
cbar = fig2.colorbar(pm,
                     ax=fig2.axes,
                     orientation="horizontal",
                     fraction=0.05, pad=0.08,
                     ticks=[(CAT_BOUNDS[i] + CAT_BOUNDS[i + 1]) / 2 for i in range(len(CAT_LABELS))]
                     )
cbar.ax.set_xticklabels(CAT_LABELS, rotation=15, ha="right", fontsize=9)
cbar.set_label("Peak WBGT heat stress category", fontsize=11)

fig2.suptitle("Weekly WBGT Heat Stress Outlooks",
              fontsize=16, fontweight="bold")

out2 = f"plots/{horizon[VAR]}_weekly_forecast_categories_{today.strftime('%Y-%m-%d')}.png"
plt.savefig(out2, dpi=200, bbox_inches="tight")
print(f"Saved: {out2}")
