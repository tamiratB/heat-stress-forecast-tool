#!/usr/bin/env python3
""""
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
Developed by: Tamirat. B. Jimma
"""

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch, Patch
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

# updated to postprocess both extremes at once; drop one from
# the list to process it alone.
VARIABLES = ["WBGT_tmax", "WBGT_tmin"]

# timeshift
TZ_SHIFT_HOURS = 3          # UTC -> East Africa Time

# domain extent
EXTENT = [21, 52, -12.5, 23.5]

# operational heat-stress categories (WBGT, degC). Bands follow the widely
# used occupational / sport WBGT flag system (ISO 7243-style reference
# limits); adjust to your local guidance if needed.
CAT_BOUNDS = [0, 25, 28, 30, 32, 60]

# -------------------------------------------------------------------------------
today = datetime.now()

print(f"\nPost-processing forecasts for WBGT_forecast_output_{today.strftime('%Y-%m-%d')}.nc...")

# input file (defaults to today's run; override as needed)
NC_FILE = f"WBGT_forecast_output_{today.strftime('%Y-%m-%d')}.nc"

week1_start = today
week1_end = today + timedelta(days=6)

week2_start = today + timedelta(days=7)
week2_end = today + timedelta(days=13)

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

MIN_STEPS_PER_DAY = 3

proj = ccrs.PlateCarree()

CATEGORIES = [
    ("Low",       "< 25 °C",     "normal activity",       "#4daf4a"),
    ("Moderate",  "25 - 28 °C",  "stay hydrated",         "#ffff33"),
    ("High",      "28 - 30 °C",  "limit exertion",        "#ff7f00"),
    ("Very high", "30 - 32 °C",  "curtail outdoor work",  "#e41a1c"),
    ("Extreme",   "≥ 32 °C",     "avoid outdoor work",    "#7f0000"),
]

CAT_NAMES = [c[0] for c in CATEGORIES]
CAT_RANGES = [c[1] for c in CATEGORIES]
CAT_ADVICE = [c[2] for c in CATEGORIES]
CAT_COLORS = [c[3] for c in CATEGORIES]

assert len(CATEGORIES) == len(CAT_BOUNDS) - 1, \
    "CATEGORIES and CAT_BOUNDS are out of step"

HEAT_CMAP = cmaps.gui_default(np.linspace(0, 1, 22))
FREQ_CMAP = cmaps.prcp_1(np.linspace(0, 1, 17))
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
    steps_per_day = (xr.ones_like(da_week.isel(lat=0, lon=0, drop=True))
                     .resample(time="1D").sum())
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

A_levels = np.linspace(23, 33, 11)
B_levels = np.linspace(23, 33, 11)
freq_levels = np.arange(0, 101, 10)

horizon = {
    "WBGT_tmax": "Daytime",
    "WBGT_tmin": "Nighttime"
    }


def plot_summary(var, da, metrics):
    """3-metric x n-week overview figure (mean, maximum, exceedance frequency)."""
    cols = [
        (f"{horizon[var]} (weekly mean)", HEAT_CMAP, A_levels, "WBGT (°C)"),
        (f"{horizon[var]} (maximum)", HEAT_CMAP, B_levels, "WBGT (°C)"),
        (f"Frequency WBGT > {THRESHOLD:.0f}°C", FREQ_CMAP, freq_levels, "% of time"),
    ]

    nrows, ncols = len(WEEKS), 3
    fig = plt.figure(figsize=(3.5 * ncols, 4.2 * nrows))
    ims = [None, None, None]

    for r, (wname, m) in enumerate(metrics.items()):
        for c in range(3):
            ax = fig.add_subplot(nrows, ncols, r * ncols + c + 1, projection=proj)
            title, cmap, levels, _ = cols[c]
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
        cb = fig.colorbar(ims[c],
                          cax=cax,
                          orientation="horizontal"
                          )
        cb.set_label(cols[c][3], fontsize=11)

    fig.suptitle(f"Weekly {horizon[var]} Heat Stress Summary",
                 fontsize=17,
                 fontweight="bold",
                 y=1.0
                 )
    fig.subplots_adjust(left=0.06,
                        right=0.98,
                        top=0.93,
                        bottom=0.10,
                        wspace=0.12,
                        hspace=0.15
                        )

    out = f"plots/{horizon[var]}_weekly_forecast_summary_{today.strftime('%Y-%m-%d')}.png"
    fig.savefig(out,
                dpi=200,
                bbox_inches="tight"
                )
    plt.close(fig)
    print(f"Saved: {out}")


def plot_categories(var, da, metrics):
    """Operational category map of the weekly peak WBGT, one panel per week."""
    fig = plt.figure(figsize=(5.4 * len(WEEKS) + 1.2, 7.0))
    map_axes = []
    for i, (wname, m) in enumerate(metrics.items()):
        ax = fig.add_subplot(1, len(WEEKS), i + 1, projection=proj)
        map_axes.append(ax)
        pm = ax.contourf(da.lon,
                         da.lat,
                         m[1],
                         levels=CAT_BOUNDS,
                         colors=CAT_COLORS,
                         extend="neither",
                         transform=proj
                         )
        clip_to_region(ax, pm)
        add_map_features(ax)
        ax.set_title(wname, fontsize=14, fontweight="bold")

    fig.subplots_adjust(left=0.04,
                        right=0.88,
                        top=0.90,
                        bottom=0.19,
                        wspace=0.10
                        )

    cat_cmap = ListedColormap(CAT_COLORS)
    cat_norm = BoundaryNorm(np.arange(len(CATEGORIES) + 1), cat_cmap.N)

    fig.canvas.draw()
    panel = map_axes[-1].get_position()

    cax = fig.add_axes([panel.x1 + 0.025, panel.y0, 0.016, panel.height])
    cbar = fig.colorbar(plt.cm.ScalarMappable(cmap=cat_cmap, norm=cat_norm),
                        cax=cax,
                        orientation="vertical",
                        ticks=np.arange(len(CATEGORIES)) + 0.5
                        )
    cbar.ax.set_yticklabels(CAT_NAMES, fontsize=11)
    cbar.ax.tick_params(length=0)
    cbar.outline.set_linewidth(0.6)
    cbar.set_label("Heat stress category",
                   fontsize=11,
                   fontweight="bold",
                   labelpad=10
                   )

    legend_handles = [
        Patch(facecolor=color,
              edgecolor="black",
              linewidth=0.6,
              label=f"{rng}  :  {advice}")
        for (_, rng, advice, color) in CATEGORIES
    ]

    LEGEND_ROWS = 2

    ncol = int(np.ceil(len(legend_handles) / LEGEND_ROWS))
    legend_handles = [legend_handles[r * ncol + c]
                      for c in range(ncol)
                      for r in range(LEGEND_ROWS)
                      if r * ncol + c < len(legend_handles)]

    leg = fig.legend(handles=legend_handles,
                     loc="lower center",
                     bbox_to_anchor=(0.46, 0.01),
                     ncol=ncol,
                     fontsize=10,
                     frameon=True,
                     handlelength=1.4,
                     handleheight=1.0,
                     columnspacing=1.6,
                     title="WBGT range and recommended action",
                     title_fontsize=11,
                     )
    leg.get_title().set_fontweight("bold")

    fig.suptitle(f"Weekly {horizon[var]} WBGT Heat Stress Outlooks",
                 fontsize=16, fontweight="bold")

    out = f"plots/{horizon[var]}_weekly_forecast_categories_{today.strftime('%Y-%m-%d')}.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}")


os.makedirs('plots', exist_ok=True)

for VAR in VARIABLES:
    if VAR not in ds:
        print(f"  {VAR} not in {NC_FILE} - skipped")
        continue

    print(f"\n{horizon[VAR]} ({VAR}):")
    da = ds[VAR]
    week_data = {name: da.sel(time=slice(start, end))
                 for name, (start, end) in WEEKS.items()}
    metrics = {name: weekly_metrics(sub) for name, sub in week_data.items()}

    plot_summary(VAR, da, metrics)
    plot_categories(VAR, da, metrics)
