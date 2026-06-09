#!/usr/bin/env python3
"""
modify_wrf_mid_lowcloud.py

Modify WRF restart (wrfrst) and boundary (wrfbdy) files to increase low cloud
coverage by enhancing moisture and cloud water in the lower troposphere.

This is the intermediate version between the original (no modification) and
modify_wrf_more_lowcloud.py (QVAPOR×1.20, QCLOUD×2.00, RH_seed=95%).

Modification parameters (halfway between original and more_lowcloud):
  QVAPOR_SCALE      : ×1.10  (was ×1.20)
  QCLOUD_SCALE      : ×1.50  (was ×2.00)
  RH_SEED_THRESHOLD : 0.97   (was 0.95 — requires higher RH to seed new cloud)
  QCLOUD_SEED       : 5.0e-6 (was 1.0e-5)

Usage:
  /global/homes/h/heroplr/.conda/envs/esmac_diags/bin/python modify_wrf_mid_lowcloud.py
"""

import xarray as xr
import numpy as np

# ============================================================
# Paths
# ============================================================

#DATA_DIR   = "/pscratch/sd/h/heroplr/Data_To_Yan/20170714case/"
DATA_DIR   = "/pscratch/sd/h/heroplr/RGMA_homecopy/WRF_dm_Stam3/test/20170714_mid_3aer/"
OUT_DIR    = DATA_DIR

WRFRST_IN  = DATA_DIR + "wrfrst_d01_2017-07-15_07:00:00"
WRFBDY_IN  = DATA_DIR + "wrfbdy_d01"

WRFRST_OUT = OUT_DIR  + "wrfrst_d01_2017-07-15_07:00:00_midlowcloud"
WRFBDY_OUT = OUT_DIR  + "wrfbdy_d01_midlowcloud"

# ============================================================
# Modification parameters
# ============================================================

K_LOW         = 28     # levels 0–K_LOW-1 are in the low-cloud target layer
K_TAPER_START = 23     # full scaling below this; taper to 1.0 over [K_TAPER_START, K_LOW)

QVAPOR_SCALE  = 1.10   # halfway between 1.0 (original) and 1.20 (more_lowcloud)
QCLOUD_SCALE  = 1.50   # halfway between 1.0 (original) and 2.00 (more_lowcloud)

# RH-based cloud seeding (creates new QCLOUD where air is nearly saturated)
RH_SEED_THRESHOLD = 0.97    # more conservative than more_lowcloud (0.95); less than original (no seeding)
QCLOUD_SEED       = 5.0e-6  # kg kg-1 – half of more_lowcloud (1.0e-5)
QCLOUD_THRESH     = 1.0e-7  # kg kg-1 – threshold for CLDFRA = 1 diagnosis (unchanged)

# WRF thermodynamic constants
P00   = 1.0e5          # reference pressure [Pa]
Rd_Cp = 0.28571428     # R_d / C_p  (used for Exner function)

# ============================================================
# Helper functions
# ============================================================

def build_scale_profile(n_levels, k_low, k_taper_start, full_scale):
    """
    1-D scale profile of length n_levels:
      k < k_taper_start          : full_scale
      k_taper_start <= k < k_low : linear taper from full_scale → 1.0
      k >= k_low                 : 1.0
    """
    p = np.ones(n_levels, dtype=np.float32)
    for k in range(k_taper_start):
        p[k] = full_scale
    for k in range(k_taper_start, k_low):
        t = (k - k_taper_start) / float(k_low - k_taper_start)
        p[k] = float(full_scale) * (1.0 - t) + 1.0 * t
    return p


def scale_3d_var(data_np, scale_profile, clip_min=None, clip_max=None):
    """
    In-place scale of a (Time, bottom_top, SN, WE) numpy array.
    scale_profile : 1-D array of length bottom_top.
    """
    for k in range(data_np.shape[1]):
        s = float(scale_profile[k])
        if s != 1.0:
            data_np[:, k, :, :] *= s
    if clip_min is not None:
        np.clip(data_np, clip_min, np.inf, out=data_np)
    if clip_max is not None:
        np.clip(data_np, -np.inf, clip_max, out=data_np)
    return data_np


def scale_bdy_var(data_np, scale_profile, clip_min=None, clip_max=None):
    """
    In-place scale of a (Time, bdy_width, bottom_top, bdy_len) numpy array.
    scale_profile : 1-D array of length bottom_top.
    """
    for k in range(data_np.shape[2]):
        s = float(scale_profile[k])
        if s != 1.0:
            data_np[:, :, k, :] *= s
    if clip_min is not None:
        np.clip(data_np, clip_min, np.inf, out=data_np)
    if clip_max is not None:
        np.clip(data_np, -np.inf, clip_max, out=data_np)
    return data_np


def compute_rh(qvapor, T_pert, P_pert, P_base):
    """
    Compute relative humidity from WRF prognostic fields.

    Parameters
    ----------
    qvapor  : (1, nz, ny, nx) float32  – water vapour mixing ratio [kg kg-1]
    T_pert  : (1, nz, ny, nx) float32  – perturbation potential temperature [K]
                                          (T in wrfrst = theta - 300)
    P_pert  : (1, nz, ny, nx) float32  – perturbation pressure [Pa]
    P_base  : (1, nz, ny, nx) float32  – base-state pressure [Pa]

    Returns
    -------
    rh : (1, nz, ny, nx) float32  – relative humidity (0–1; can exceed 1)
    """
    P_total = P_pert + P_base
    theta   = T_pert.astype(np.float64) + 300.0
    T_act   = theta * (P_total.astype(np.float64) / P00) ** Rd_Cp

    e_s   = 610.78 * np.exp(17.2693882 * (T_act - 273.16) / (T_act - 35.86))
    q_sat = 0.622 * e_s / (P_total.astype(np.float64) - 0.378 * e_s)

    rh = qvapor.astype(np.float64) / np.maximum(q_sat, 1e-12)
    return rh.astype(np.float32)


# ============================================================
# Print configuration
# ============================================================

print("=" * 60)
print("WRF Low-Cloud Enhancement Script  [MID version]")
print("=" * 60)
print(f"  Target levels  : 0–{K_LOW-1}  (~surface to ~3 km)")
print(f"  Taper zone     : levels {K_TAPER_START}–{K_LOW-1}")
print(f"  QVAPOR scale   : ×{QVAPOR_SCALE}  (more_lowcloud: ×1.20)")
print(f"  QCLOUD scale   : ×{QCLOUD_SCALE}  (more_lowcloud: ×2.00)")
print(f"  RH seed thresh : {RH_SEED_THRESHOLD*100:.0f}%  →  QCLOUD_seed = {QCLOUD_SEED} kg/kg  (more_lowcloud: 95%, 1e-5)")
print()

# ============================================================
# 1.  Modify wrfrst (restart file)
# ============================================================

print(f"Opening wrfrst: {WRFRST_IN}")
rst = xr.open_dataset(WRFRST_IN)
n_lev = rst.sizes["bottom_top"]
print(f"  Dimensions: {rst.sizes['west_east']} × {rst.sizes['south_north']} × {n_lev}")
print(f"  Total variables: {len(rst.variables)}")
print()

qv_prof = build_scale_profile(n_lev, K_LOW, K_TAPER_START, QVAPOR_SCALE)
qc_prof = build_scale_profile(n_lev, K_LOW, K_TAPER_START, QCLOUD_SCALE)

# ── QVAPOR ──────────────────────────────────────────────────
print("  [1/6] Scaling QVAPOR ...")
qv = rst["QVAPOR"].values.copy()
before = qv[0, :K_LOW].mean()
qv = scale_3d_var(qv, qv_prof, clip_min=0.0)
after  = qv[0, :K_LOW].mean()
print(f"        low-level mean: {before*1000:.2f} → {after*1000:.2f} g/kg  (×{after/before:.3f})")

# ── QCLOUD ──────────────────────────────────────────────────
print("  [2/6] Scaling QCLOUD ×1.5 in existing cloud regions ...")
qc = rst["QCLOUD"].values.copy()
before_qc = qc[0, :K_LOW].mean()
qc = scale_3d_var(qc, qc_prof, clip_min=0.0)
after_qc  = qc[0, :K_LOW].mean()
print(f"        low-level mean: {before_qc*1e6:.4f} → {after_qc*1e6:.4f} mg/kg")

# ── RH-based QCLOUD seeding ─────────────────────────────────
print(f"  [3/6] RH-based QCLOUD seeding (RH > {RH_SEED_THRESHOLD*100:.0f}%, levels 0–{K_LOW-1}) ...")
T_pert = rst["T"].values.copy()
P_pert = rst["P"].values.copy()
P_base = rst["PB"].values.copy()
rh     = compute_rh(qv, T_pert, P_pert, P_base)

n_seeded = 0
for k in range(K_LOW):
    mask = (rh[0, k] > RH_SEED_THRESHOLD) & (qc[0, k] < QCLOUD_THRESH)
    if mask.any():
        qc[0, k, mask] = QCLOUD_SEED
        n_seeded += int(mask.sum())

print(f"        Grid points seeded: {n_seeded:,}")
print(f"        QCLOUD low mean after seeding: {qc[0,:K_LOW].mean()*1e6:.4f} mg/kg")
del T_pert, P_pert, P_base, rh

# ── CLDFRA ──────────────────────────────────────────────────
print("  [4/6] Re-diagnosing CLDFRA from QCLOUD_new ...")
cf_orig = rst["CLDFRA"].values.copy()
cf_new  = np.zeros_like(cf_orig)
cf_new[0, :K_LOW][qc[0, :K_LOW] > QCLOUD_THRESH] = 1.0
cf_new[0, K_LOW:] = cf_orig[0, K_LOW:]
orig_cloud_frac = float(cf_orig[0, :K_LOW].mean())
new_cloud_frac  = float(cf_new[0, :K_LOW].mean())
print(f"        Low-level cloud fraction: {orig_cloud_frac:.4f} → {new_cloud_frac:.4f}")

# ── CLDFRA_OLD ──────────────────────────────────────────────
print("  [5/6] Setting CLDFRA_OLD = CLDFRA_new ...")
cf_old_new = cf_new.copy()
if "CLDFRA_OLD" in rst:
    cf_old_new[0, K_LOW:] = rst["CLDFRA_OLD"].values[0, K_LOW:]

# ── Q2 ──────────────────────────────────────────────────────
print("  [6/6] Scaling Q2 (surface 2-m moisture) ...")
q2 = rst["Q2"].values.copy()
q2 *= QVAPOR_SCALE
q2 = np.clip(q2, 0.0, None)
print(f"        Q2 scaled by ×{QVAPOR_SCALE}")

# ── Assign modified variables back to dataset ────────────────
print()
print("  Assigning modified variables to dataset ...")

def assign_back(ds, name, data_np):
    ds[name] = xr.DataArray(data_np, dims=ds[name].dims,
                             coords=ds[name].coords, attrs=ds[name].attrs)

assign_back(rst, "QVAPOR",    qv)
assign_back(rst, "QCLOUD",    qc)
assign_back(rst, "CLDFRA",    cf_new)
if "CLDFRA_OLD" in rst:
    assign_back(rst, "CLDFRA_OLD", cf_old_new)
assign_back(rst, "Q2",        q2)

print(f"\n  Writing wrfrst output: {WRFRST_OUT}")
rst.to_netcdf(WRFRST_OUT)
rst.close()
print("  wrfrst done.\n")

# ============================================================
# 2.  Modify wrfbdy (boundary file)
# ============================================================

print(f"Opening wrfbdy: {WRFBDY_IN}")
bdy = xr.open_dataset(WRFBDY_IN)
n_lev_bdy = bdy.sizes["bottom_top"]
print(f"  Time steps: {bdy.sizes['Time']}, bdy_width: {bdy.sizes['bdy_width']}, "
      f"bottom_top: {n_lev_bdy}")
print(f"  Total variables: {len(bdy.variables)}")
print()

qv_bdy_prof = build_scale_profile(n_lev_bdy, K_LOW, K_TAPER_START, QVAPOR_SCALE)
qc_bdy_prof = build_scale_profile(n_lev_bdy, K_LOW, K_TAPER_START, QCLOUD_SCALE)

SIDES = [("BXS", "BTXS"), ("BXE", "BTXE"), ("BYS", "BTYS"), ("BYE", "BTYE")]

print("  Scaling QVAPOR boundaries and tendencies ...")
for val_sfx, tend_sfx in SIDES:
    val_name  = f"QVAPOR_{val_sfx}"
    tend_name = f"QVAPOR_{tend_sfx}"

    val = bdy[val_name].values.copy()
    before = val[:, :, :K_LOW, :].mean()
    val = scale_bdy_var(val, qv_bdy_prof, clip_min=0.0)
    after  = val[:, :, :K_LOW, :].mean()
    print(f"    {val_name}: low mean {before:.1f} → {after:.1f}  (ratio {after/before:.3f})")

    bdy[val_name] = xr.DataArray(val, dims=bdy[val_name].dims,
                                  coords=bdy[val_name].coords, attrs=bdy[val_name].attrs)

    if tend_name in bdy:
        bt = bdy[tend_name].values.copy()
        bt = scale_bdy_var(bt, qv_bdy_prof)
        bdy[tend_name] = xr.DataArray(bt, dims=bdy[tend_name].dims,
                                       coords=bdy[tend_name].coords, attrs=bdy[tend_name].attrs)

print()
print("  Scaling QCLOUD boundaries and tendencies ...")
for val_sfx, tend_sfx in SIDES:
    val_name  = f"QCLOUD_{val_sfx}"
    tend_name = f"QCLOUD_{tend_sfx}"

    val = bdy[val_name].values.copy()
    before = val[:, :, :K_LOW, :].mean()
    val = scale_bdy_var(val, qc_bdy_prof, clip_min=0.0)
    after  = val[:, :, :K_LOW, :].mean()
    print(f"    {val_name}: low mean {before:.4f} → {after:.4f}  (ratio x{QCLOUD_SCALE})")

    bdy[val_name] = xr.DataArray(val, dims=bdy[val_name].dims,
                                  coords=bdy[val_name].coords, attrs=bdy[val_name].attrs)

    if tend_name in bdy:
        bt = bdy[tend_name].values.copy()
        bt = scale_bdy_var(bt, qc_bdy_prof)
        bdy[tend_name] = xr.DataArray(bt, dims=bdy[tend_name].dims,
                                       coords=bdy[tend_name].coords, attrs=bdy[tend_name].attrs)

print(f"\n  Writing wrfbdy output: {WRFBDY_OUT}")
bdy.to_netcdf(WRFBDY_OUT)
bdy.close()
print("  wrfbdy done.\n")

# ============================================================
# Summary
# ============================================================
print("=" * 60)
print("Summary of modifications  [MID version]")
print("=" * 60)
print(f"  wrfrst output : {WRFRST_OUT}")
print(f"    QVAPOR      : ×{QVAPOR_SCALE} in levels 0–{K_LOW-1} (taper {K_TAPER_START}–{K_LOW-1})")
print(f"    QCLOUD      : ×{QCLOUD_SCALE} in existing clouds + RH-seed at RH>{RH_SEED_THRESHOLD*100:.0f}%")
print(f"    CLDFRA      : re-diagnosed (=1 where QCLOUD_new > {QCLOUD_THRESH})")
print(f"    CLDFRA_OLD  : = CLDFRA_new")
print(f"    Q2          : ×{QVAPOR_SCALE}")
print()
print(f"  wrfbdy output : {WRFBDY_OUT}")
print(f"    QVAPOR_B{{XS,XE,YS,YE}} : ×{QVAPOR_SCALE} in levels 0–{K_LOW-1}")
print(f"    QVAPOR_BT{{XS,XE,YS,YE}}: tendencies scaled consistently")
print(f"    QCLOUD_B{{XS,XE,YS,YE}} : ×{QCLOUD_SCALE} in levels 0–{K_LOW-1}")
print(f"    QCLOUD_BT{{XS,XE,YS,YE}}: tendencies scaled consistently")
print()
print("Done.")
