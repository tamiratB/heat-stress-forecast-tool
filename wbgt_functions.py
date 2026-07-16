"""
wbgt_functions.py
=================
Robust computation of the outdoor Wet-Bulb Globe Temperature (WBGT)
from gridded daily meteorological data.

WBGT formulation
---------------------------------------------
    WBGT = 0.7 · Tw  +  0.2 · Tg  +  0.1 · Ta

where
    Tw  – natural wet-bulb temperature   [K]
    Tg  – black globe temperature        [K]
    Ta  – 2m air temperature   [K]
---------------------------------------------

Developed by: Tamirat B. Jimma
"""

import logging
import os
from typing import Optional

import numpy as np
import pandas as pd
import xarray as xr
import warnings
import dask

logger = logging.getLogger(__name__)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

# physical constants
R_AIR = 287.058          # J kg-1 K-1  specific gas constant for dry air
CP_AIR = 1005.0          # J kg-1 K-1  specific heat at constant pressure
LV = 2.501e6             # J kg-1      latent heat of vaporisation (0 °C)
R_V = 461.5              # J kg-1 K-1  specific gas constant for water vapour
SIGMA = 5.670374419e-8   # W m-2 K-4   Stefan-Boltzmann constant
I_SOLAR = 1367.0         # W m-2       solar constant
PR_AIR = 0.71            # –           Prandtl number
G = 9.80665              # m s-2       standard acceleration due to gravity
M_H2O = 0.018015         # [kg mol-1]  molar mass of water
M_AIR = 0.028964         # [kg mol-1]  molar mass of dry air
SC_AIR = 0.60            # [-]         Schmidt number of air


def nc_write_encoding(ds, compress=False, least_significant_digit=None):
    enc = {}
    for name, da in ds.data_vars.items():
        e = {"dtype": "float32"}
        if da.chunks is not None:
            e["chunksizes"] = tuple(c[0] for c in da.chunks)
        if compress:
            e.update(zlib=True, complevel=1, shuffle=True)
            if least_significant_digit is not None:
                e["least_significant_digit"] = least_significant_digit
        enc[name] = e
    return enc


def _safe_divide(
    num: xr.DataArray,
    den: xr.DataArray,
    floor: float = 1e-6,
) -> xr.DataArray:

    safe_den = xr.where(np.abs(den) < floor, floor, den)
    return num / safe_den


def es_hpa_from_Tk(Tk: xr.DataArray) -> xr.DataArray:

    Tc = Tk - 273.15
    return 6.112 * np.exp(17.67 * Tc / (Tc + 243.5))


def prata_emissivity(
    ea_hpa: xr.DataArray,
    Ta_K: xr.DataArray,
) -> xr.DataArray:

    w = 46.5 * ea_hpa / Ta_K

    epsilon_a = 1.0 - (1.0 + w) * np.exp(
        -np.sqrt(1.2 + 3.0 * w)
    )

    return epsilon_a.clip(0.0, 1.0)


def psychrometric_stull_wbt(
    Ta: xr.DataArray,
    RH: xr.DataArray,
) -> xr.DataArray:

    Tc = Ta - 273.15
    RH = RH.clip(5.0, 99.0)
    Tw_c = (
        Tc * np.arctan(0.151977 * (RH + 8.313659) ** 0.5)
        + np.arctan(Tc + RH)
        - np.arctan(RH - 1.676331)
        + 0.00391838 * (RH ** 1.5) * np.arctan(0.023101 * RH)
        - 4.686035
    )
    return Tw_c + 273.15


def _accumulation_seconds(ssrd: xr.DataArray) -> float:

    if ssrd.time.size < 2:
        return 3600.0

    dt_ns = np.diff(ssrd.time.values).astype("timedelta64[s]").astype(np.float64)
    median_s = float(np.median(dt_ns))

    if median_s <= 0:
        raise ValueError(
            f"Time axis appears non-monotonic (median change in  time is = {median_s} s). "
            "Check your input dataset."
        )

    return median_s


def compute_instantaneous_solar(
    ssrd: xr.DataArray,
    lat_name: str = "lat",
    lon_name: str = "lon",
    debug: bool = False,
) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:

    if lat_name not in ssrd.coords:
        raise ValueError(
            f"Latitude coordinate '{lat_name}' not found in ssrd. "
            f"Available coords: {list(ssrd.coords)}"
        )
    if lon_name not in ssrd.coords:
        raise ValueError(
            f"Longitude coordinate '{lon_name}' not found in ssrd. "
            f"Available coords: {list(ssrd.coords)}"
        )
    if "time" not in ssrd.dims:
        raise ValueError(
            "ssrd must have a 'time' dimension. "
            f"Found dimensions: {list(ssrd.dims)}"
        )

    step_s = _accumulation_seconds(ssrd)

    if step_s <= 0:
        raise ValueError(f"Invalid accumulation period: {step_s} s")
    GHI = (ssrd / step_s).clip(min=0.0)

    t = pd.DatetimeIndex(ssrd.time.values)
    days_in_year = np.where(t.is_leap_year, 366.0, 365.0)

    doy = xr.DataArray(
        t.dayofyear.astype(np.float64),
        dims="time",
        coords={"time": ssrd.time},
    )

    days_in_year_da = xr.DataArray(
        days_in_year,
        dims="time",
        coords={"time": ssrd.time},
    )

    utc_hour = xr.DataArray(
        t.hour + t.minute / 60.0 + t.second / 3600.0,
        dims="time",
        coords={"time": ssrd.time},
    )

    gamma = 2.0 * np.pi * (doy - 1.0) / days_in_year_da

    decl = (
        0.006918
        - 0.399912 * np.cos(gamma)
        + 0.070257 * np.sin(gamma)
        - 0.006758 * np.cos(2 * gamma)
        + 0.000907 * np.sin(2 * gamma)
        - 0.002697 * np.cos(3 * gamma)
        + 0.001480 * np.sin(3 * gamma)
    )

    eot = 229.18 * (
        0.000075
        + 0.001868 * np.cos(gamma)
        - 0.032077 * np.sin(gamma)
        - 0.014615 * np.cos(2 * gamma)
        - 0.040849 * np.sin(2 * gamma)
    )

    E0 = (
        1.000110
        + 0.034221 * np.cos(gamma)
        + 0.001280 * np.sin(gamma)
        + 0.000719 * np.cos(2 * gamma)
        + 0.000077 * np.sin(2 * gamma)
    )

    I0_E0 = (I_SOLAR * E0)

    lat_rad = np.deg2rad(ssrd[lat_name])
    lon_deg = ssrd[lon_name]

    lat3 = lat_rad.broadcast_like(ssrd)
    lon3 = lon_deg.broadcast_like(ssrd)
    decl3 = decl.broadcast_like(ssrd)
    eot3 = eot.broadcast_like(ssrd)
    utc3 = utc_hour.broadcast_like(ssrd)
    I0_E03 = I0_E0.broadcast_like(ssrd)

    solar_time = utc3 + lon3 / 15.0 + eot3 / 60.0
    H = np.deg2rad(15.0 * (solar_time - 12.0))

    cos_theta = (
        np.sin(lat3) * np.sin(decl3)
        + np.cos(lat3) * np.cos(decl3) * np.cos(H)
    )

    cos_theta_safe = cos_theta.clip(0.0, 1.0)

    theta = xr.apply_ufunc(
        np.arccos,
        cos_theta_safe,
        dask="allowed",
    )

    S_max = (I0_E03 * cos_theta_safe).clip(min=0.0)
    S_star = _safe_divide(GHI, S_max).clip(0.0, 1.0)
    zenith_limit = np.deg2rad(89.5)

    S_star_safe = S_star.clip(min=1e-4)
    fdir_raw = xr.where(
        (cos_theta > np.cos(zenith_limit)) & (S_star > 1e-4),
        np.exp(3.0 - 1.34 * S_star - 1.65 / S_star_safe),
        0.0,
    ).clip(0.0, 1.0)

    fdir = fdir_raw
    daylight = cos_theta > 0.01

    GHI = GHI.where(daylight, 0.0)
    fdir = fdir.where(daylight, 0.0)
    theta = theta.where(daylight, np.pi / 2.0)

    if debug:
        print(f"[solar] accumulation period  : {step_s:.0f} s  "
              f"({step_s/3600:.2f} h)")
        print(f"[solar] GHI  max             : {float(GHI.max()):.2f} W m-2")
        print(f"[solar] S_max max            : {float(S_max.max()):.2f} W m-2")
        print(f"[solar] S*   range (daytime) : "
              f"{float(S_star.where(daylight, np.nan).min()):.4f} – "
              f"{float(S_star.where(daylight, np.nan).max()):.4f}")
        print(f"[solar] fdir max             : {float(fdir.max()):.4f}")
        print(f"[solar] theta range          : "
              f"{float(theta.min()):.4f} – {float(theta.max()):.4f} rad  "
              f"({np.rad2deg(float(theta.min())):.1f}° – "
              f"{np.rad2deg(float(theta.max())):.1f}°)")
        n_day = int(daylight.sum())
        n_total = int(daylight.size)
        print(f"[solar] daylight cells       : {n_day} / {n_total} "
              f"({100*n_day/n_total:.1f}%)")

    return (
        GHI.astype(np.float32),
        fdir.astype(np.float32),
        theta.astype(np.float32),
    )


def _mu_air(Tf: xr.DataArray) -> xr.DataArray:

    return 1.458e-6 * Tf**1.5 / (Tf + 110.4)


def _k_air(Tf: xr.DataArray) -> xr.DataArray:

    return (
        1.5207e-11 * Tf*Tf*Tf
        - 4.8574e-8 * Tf*Tf
        + 1.0184e-4 * Tf
        - 3.9333e-4
    )


def compute_Tw(
    Td: xr.DataArray,
    Ta: xr.DataArray,
    RH: xr.DataArray,
    P: xr.DataArray,
    U: xr.DataArray,
    S: xr.DataArray,
    theta: xr.DataArray,
    fdir: xr.DataArray,
    Dw: float = 0.007,
    Lw: float = 0.0254,
    alpha_w: float = 0.4,
    epsilon_w: float = 0.95,
    alpha_sfc: float = 0.45,
    max_iter: int = 50,
    tol: float = 1e-3,
    dask_iters: int = 20,
    debug: bool = False,
) -> xr.DataArray:

    D_4L = Dw / (4.0 * Lw)

    RH = RH.clip(0.0, 100.0)
    U = U.clip(min=0.05)          # avoid Re → 0 in calm conditions

    valid = (
        np.isfinite(Td) & np.isfinite(Ta) & np.isfinite(RH)
        & np.isfinite(P) & np.isfinite(U)
        & np.isfinite(S) & np.isfinite(theta) & np.isfinite(fdir)
        & (Ta > 200.0) & (Ta < 340.0)
        & (Td <= Ta)
        & (P > 100.0)
    )

    es_Ta = es_hpa_from_Tk(Ta)
    ea = (RH / 100.0) * es_Ta

    epsilon_a = prata_emissivity(ea, Ta)

    Ta2 = Ta * Ta
    L_in_coeff = 0.5 * (1.0 + epsilon_a) * SIGMA * (Ta2 * Ta2)

    S_eff = xr.where(theta < np.pi / 2.0, S, 0.0)

    cos_t = np.cos(theta).clip(min=1e-4)
    sin_t = np.sin(theta).clip(min=0.0)
    tan_t = (sin_t / cos_t).clip(0.0, 20.0)

    fdir_safe = xr.where(theta < np.deg2rad(89.5), fdir, 0.0)

    S_abs = (1.0 - alpha_w) * S_eff * (
          (1.0 - fdir_safe) * (1.0 + D_4L)
          + fdir_safe * (tan_t / np.pi + D_4L)
          + alpha_sfc
         )
    S_abs = S_abs.clip(min=0.0)

    psy_ratio = (1.0 / (CP_AIR * M_AIR)) * (PR_AIR / SC_AIR) ** 0.56

    M_psy = M_H2O * psy_ratio

    Tpsy = psychrometric_stull_wbt(Ta, RH)
    Tw = Tpsy.clip(min=Td, max=Ta)

    _churchill_chu_denom = (1.0 + (0.469 / PR_AIR)**(9.0 / 16.0))**(4.0 / 9.0)

    def _convection(Tw_state):

        Tf = 0.5 * (Tw_state + Ta)
        mu = _mu_air(Tf)
        k = _k_air(Tf)
        rho = (P * 100.0) / (R_AIR * Tf)

        Re = (rho * U * Dw / mu).clip(min=1.0)
        Nu_forced = 0.281 * Re**0.6 * PR_AIR**0.44

        nu_kin = mu / rho
        dT = np.abs(Tw_state - Ta).clip(min=0.1)
        Gr = (G * (1.0 / Tf) * dT * Dw**3 / (nu_kin * nu_kin)).clip(min=1.0)
        Ra = Gr * PR_AIR
        Nu_nat = 2.0 + 0.589 * np.sqrt(np.sqrt(Ra)) / _churchill_chu_denom

        Nu = np.cbrt(Nu_forced * Nu_forced * Nu_forced
                     + Nu_nat * Nu_nat * Nu_nat)
        h_c = Nu * k / Dw

        return rho, Re, Nu_forced, Gr, Ra, Nu_nat, Nu, h_c

    def _lv(T):
        return LV - 2.36e3 * (T - 273.15)

    max_delta = np.inf

    lazy = Tw.chunks is not None
    n_iters = dask_iters if lazy else max_iter

    for i in range(n_iters):

        rho, Re, Nu_forced, Gr, Ra, Nu_nat, Nu, h_c = _convection(Tw)

        es_Tw = es_hpa_from_Tk(Tw)

        Tw2 = Tw * Tw
        Tw3 = Tw2 * Tw
        Tw4 = Tw2 * Tw2

        L_out = epsilon_w * SIGMA * Tw4
        Fnet = epsilon_w * L_in_coeff - L_out + S_abs

        lv = _lv(Tw)

        P_minus_esTw = P - es_Tw
        lv_M_psy = lv * M_psy
        lat_term = lv_M_psy * (es_Tw - ea) / P_minus_esTw

        f = (Tw - Ta) + lat_term - Fnet / h_c

        des_dTw = es_Tw * lv / (R_V * Tw2)
        dlat_dTw = lv_M_psy * des_dTw * P / (P_minus_esTw * P_minus_esTw)

        dFnet_dTw = -4.0 * epsilon_w * SIGMA * Tw3

        fp = 1.0 + dlat_dTw - dFnet_dTw / h_c
        fp = xr.where(np.abs(fp) > 1e-8, fp, np.nan)

        delta = xr.where(np.isfinite(fp), f / fp, 0.0)
        delta = delta.clip(-2.0, 2.0)

        if i >= n_iters - 2:
            delta = 0.5 * delta

        Tw_new = (Tw - delta).clip(min=Td, max=Ta + 15.0)

        if not lazy:
            abs_delta = np.abs(Tw_new - Tw).where(valid, 0.0)
            max_delta = float(abs_delta.max())

        Tw = Tw_new

        if not lazy:
            if debug:
                warnings.warn(
                    f"compute_Tw iter {i+1:2d}: max |change in Tw| = {max_delta:.5f} K",
                    RuntimeWarning, stacklevel=2,
                )

            if max_delta < tol:
                if debug:
                    warnings.warn(
                        f"compute_Tw converged in {i+1} iteration(s).",
                        RuntimeWarning, stacklevel=2,
                    )
                break

        else:
            if not lazy and max_delta >= tol:
                warnings.warn(
                    f"compute_Tw: Newton solver did not converge after {max_iter} "
                    f"iterations; max |Change in Tw| = {max_delta:.4f} K. "
                    "Results in non-converged cells may be inaccurate. "
                    "Consider increasing max_iter or checking input ranges.",
                    RuntimeWarning,
                    stacklevel=2,
                )

    Tw = xr.where(valid & np.isfinite(Tw), Tw, np.nan)

    return Tw.astype(np.float32)


def compute_Tg(
    Ta: xr.DataArray,
    U: xr.DataArray,
    S: xr.DataArray,
    theta: xr.DataArray,
    fdir: xr.DataArray,
    RH: xr.DataArray,
    P: xr.DataArray,
    Dg: float = 0.0508,
    alpha_g: float = 0.05,
    epsilon_g: float = 0.95,
    alpha_sfc: float = 0.45,
    max_iter: int = 50,
    tol: float = 1e-3,
    dask_iters: int = 20,
    debug: bool = False,
) -> xr.DataArray:

    RH = RH.clip(0.0, 100.0)
    U  = U.clip(min=0.05)

    valid = (
        np.isfinite(Ta) & np.isfinite(U) & np.isfinite(S)
        & np.isfinite(theta) & np.isfinite(fdir)
        & np.isfinite(RH) & np.isfinite(P)
        & (Ta > 200.0) & (Ta < 340.0)
        & (P > 100.0)
    )

    es_Ta = es_hpa_from_Tk(Ta)
    ea = (RH / 100.0) * es_Ta
    epsilon_a = prata_emissivity(ea, Ta)

    Ta2 = Ta * Ta
    R_long_in = 0.5 * epsilon_g * (1.0 + epsilon_a) * SIGMA * (Ta2 * Ta2)

    mu_sun = np.cos(theta).clip(min=1e-3)
    S_eff = xr.where(theta < np.pi / 2, S, 0.0)

    inv_2mu = (1.0 / (2.0 * mu_sun)).clip(max=3.0)

    R_short_in = (
        (1.0 - alpha_g)
        * S_eff
        / 2.0
        * (1.0 + (inv_2mu - 1.0) * fdir + alpha_sfc)
    ).clip(min=0.0)

    a = epsilon_g * SIGMA

    Tg = Ta.copy()

    _churchill_chu_denom = (1.0 + (0.469 / PR_AIR)**(9.0 / 16.0))**(4.0 / 9.0)

    if debug:
        logger.debug(
            "Tg solver: Ta_mean=%.2f K, R_short_in_mean=%.2f W m-2, "
            "R_long_in_mean=%.2f W m-2",
            float(Ta.mean()), float(R_short_in.mean()),
            float(R_long_in.mean()),
        )

    lazy = Tg.chunks is not None
    n_iters = dask_iters if lazy else max_iter

    for i in range(n_iters):
        Tf = 0.5 * (Tg + Ta)
        MU_AIR = 1.458e-6 * Tf**1.5 / (Tf + 110.4)
        K_AIR = _k_air(Tf)

        rho = (P * 100.0) / (R_AIR * Tf)
        Re = (rho * U * Dg / MU_AIR).clip(min=1.0)
        Nu_forced = 2.0 + 0.6 * Re**0.5 * PR_AIR**(1.0 / 3.0)

        beta = 1.0 / Tf

        nu = MU_AIR / rho

        dT = np.abs(Tg - Ta).clip(min=0.1)
        Gr = (
            G
            * beta
            * dT
            * Dg**3
            / (nu * nu)
        ).clip(min=1.0)

        Ra = Gr * PR_AIR
        Nu_nat = 2.0 + 0.589 * np.sqrt(np.sqrt(Ra)) / _churchill_chu_denom

        Nu_forced2 = Nu_forced * Nu_forced
        Nu_nat2 = Nu_nat * Nu_nat
        Nu = np.cbrt(Nu_forced2 * Nu_forced + Nu_nat2 * Nu_nat)

        h = Nu * K_AIR / Dg
        c = h * Ta + R_long_in + R_short_in

        Tg2 = Tg * Tg
        Tg3 = Tg2 * Tg
        Tg4 = Tg2 * Tg2
        f = a * Tg4 + h * Tg - c
        fp = 4.0 * a * Tg3 + h

        delta = f / fp
        delta = delta.clip(-5.0, 5.0)

        Tg_new = Tg - delta
        Tg_new = Tg_new.clip(min=200.0, max=370.0)

        if not lazy:
            max_delta = float(np.abs(Tg_new - Tg).where(valid, 0.0).max())

        Tg = Tg_new

        if not lazy:
            if debug:
                logger.debug(
                    "Tg iter %2d: max. difference is = %.5f K, h_mean = %.3f W m-2 K-2",
                    i + 1, max_delta, float(h.mean()),
                )

            if max_delta < tol:
                logger.debug("Tg converged in %d iterations.", i + 1)
                break

    Tg = xr.where(valid & np.isfinite(Tg), Tg, np.nan)

    return Tg.astype(np.float32)


def compute_WBGT(
    Td: xr.DataArray,
    Ta: xr.DataArray,
    RH: xr.DataArray,
    P_hPa: xr.DataArray,
    U: xr.DataArray,
    S_instant: xr.DataArray,
    fdir_instant: xr.DataArray,
    theta_instant: xr.DataArray,
    Dg: float = 0.0508,
    alpha_g: float = 0.05,
    epsilon_g: float = 0.95,
    alpha_sfc: float = 0.45,
    Dw: float = 0.007,
    Lw: float = 0.0254,
    alpha_w: float = 0.4,
    epsilon_w: float = 0.95,
    max_iter: int = 10,
    tol: float = 0.001,
    test_year: int = None,
    debug: bool = False,
    save_diagnostics: bool = True,
    diag_locations: Optional[dict] = None,
) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:

    Td, Ta, RH, P_hPa, U, S_instant, fdir_instant, theta_instant = dask.persist(
        Td, Ta, RH, P_hPa, U, S_instant, fdir_instant, theta_instant
    )

    logger.info("Computing Tw...")
    Tw = compute_Tw(
        Td, Ta, RH, P_hPa, U,
        S_instant, theta_instant, fdir_instant,
        Dw=Dw, Lw=Lw, alpha_w=alpha_w, epsilon_w=epsilon_w, alpha_sfc=alpha_sfc,
        max_iter=max_iter, tol=tol, debug=debug,
    )

    logger.info("Computing Tg...")
    Tg = compute_Tg(
        Ta, U,
        S_instant, theta_instant, fdir_instant,
        RH, P_hPa,
        Dg=Dg, alpha_g=alpha_g, epsilon_g=epsilon_g, alpha_sfc=alpha_sfc,
        max_iter=max_iter, tol=tol, debug=debug,
    )

    valid = np.isfinite(Tw) & np.isfinite(Tg) & np.isfinite(Ta)

    logger.info("Computing WBGT...")
    WBGT = (0.7 * Tw + 0.2 * Tg + 0.1 * Ta).where(valid)
    WBGT = WBGT.persist()

    return Tw, Tg, WBGT.astype(np.float32)
