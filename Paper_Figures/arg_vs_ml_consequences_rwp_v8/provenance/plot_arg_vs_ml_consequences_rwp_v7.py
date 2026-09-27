#!/usr/bin/env python
"""
Consequences of the ARG activation bias for cloud and aerosol properties.

Section~\\ref{subsubsec:full_record_validation} establishes that the diagnostic
ARG activated fraction sits systematically below parcel-model truth (pooled
median relative error -0.503 in the accumulation mode).  This script traces
what that bias does once ARG is the ACTIVE scheme, by comparing the four
ML-ACT runs against the four ARG-ACT runs along the causal chain

    activated fraction -> droplet number -> droplet size -> precipitation
                       -> wet removal -> aerosol / CCN -> cloud water -> albedo

Panels (5 rows x 2 columns)
    a  ac0        accumulation-mode number            [cm-3]
    b  nu0        Aitken-mode number                  [cm-3]
    c  CCN(S)     CCN at the selected supersaturation [cm-3]
    d  QNDROP     cloud droplet number                [cm-3]
    e  r_v        droplet mean-volume radius          [um]   (derived)
    f  LWP_TOT    liquid water path                   [g m-2]
    g  precip     surface precipitation rate          [mm h-1] (derived)
    h  q_c        in-cloud cloud water                [g kg-1]
    i  RWP        rain water path                     [g m-2] (derived)
    j  CLDCOL     cloudy-column fraction              [-]

PANEL (i): RAIN WATER PATH instead of TOA shortwave CRE.  There is no rain-water
    -path diagnostic in wrfout -- LWP/LWP_TOT come from the solar_diag package
    but have no rain counterpart -- so RWP is integrated here the same way the
    model integrates LWP:
        RWP = sum_k ( q_rain * rho * dz )        [kg m-2] -> x1000 -> [g m-2]
    over the FULL column, then area-averaged over cloudy columns, matching the
    LWP_TOT convention in panel (f).  Note the full-column integral: rain falls
    below cloud base, so restricting it to Z_TOP as the 3-D in-cloud scalars are
    would discard most of the mass.
    The panel is drawn on a LOG axis.  RWP spans ~0.2 g m-2 in the ARG org case
    to ~220 g m-2 in the mid cases; on a linear axis all four org curves collapse
    onto zero.
    SW CRE is still computed and still appears in the summary table; it is only
    no longer plotted.

WHERE THE FIELDS COME FROM
    The ML runs (cosine_more_fix) carry CCN1..CCN6 in wrfout.  The ARG runs
    (WRF_dm_v2) do NOT -- for those cases CCN is read from the matching wrfrst,
    which does carry CCN1..CCN6 (restart and history are written on the same
    2-minute cadence, one restart per history time).  LWP_TOT and TAU_QC_TOT
    are wrfout-only in both trees.  TAU_QC_TOT is written by both but is
    IDENTICALLY ZERO everywhere in the old WRF_dm (ARG) build -- the
    optical-depth diagnostic is not computed there.  WRF_dm_v2 does populate it,
    but it is still not plotted, so that this figure stays directly comparable
    with v3; in-cloud cloud water is shown instead.

AVERAGING (three conventions, deliberately)
    * 3-D fields (ac0, nu0, CCN, QNDROP, QCLOUD): air-mass-weighted mean below
      Z_TOP, over CLOUDY COLUMNS only -- the aerosol actually feeding the cloud.
    * 2-D cloud fields (LWP_TOT, TAU_QC_TOT): area mean over cloudy columns.
    * Domain-mean fields (precip rate, SW CRE): ALL columns.  Precipitation and
      the radiative effect are domain-integrated quantities; masking them to
      instantaneously cloudy columns would make the time series depend on how
      the cloud mask moves rather than on the physics.

DERIVED QUANTITIES
    r_v  = (3 q_c / (4 pi rho_w N))^(1/3), from the ratio of the averaged cloud
           water and averaged droplet number (ratio of means, not mean of
           ratios -- robust where N is small).  Air density cancels.
    rain = d(RAINNC + RAINC)/dt between consecutive output times.  RAINNC is
           carried through restarts, so differencing is safe across them.
    CRE  = -(SWUPT - SWUPTC), the usual sign convention (negative = cooling).

WHAT v5 CHANGES FROM v3
    Same eight curves and the same layout as v3 -- four ML-ACT runs (solid) and
    four ARG-ACT runs (dashed), one colour per case -- but the ARG side is now
    read from /scratch/.../WRF_dm_v2/test instead of /scratch/.../WRF_dm/test.
    Those are the SAME four 20170714 ARG-ACT cases restarted from the same
    09:02 restart; the tree differs from WRF_dm only in Registry.EM_COMMON
    (re_cloud/re_ice/re_snow promoted to "rh", plus a morr_two_moment package
    line), so activation, microphysics and the aerosol treatment are unchanged.
    The v4 figure carried both trees at once, WRF_dm dashed and WRF_dm_v2
    dotted, to show the reproduction; that check is done -- over the frames
    the two trees share, every field agrees to within 1.2% (all but mid_3aer
    to within 0.5%), and only SW CRE departs, by 1.9-5.5% -- so the duplicate
    lines are dropped here.

    ONE CAVEAT, and it is only about radiation.  Both trees run mp_physics=10,
    ra_sw_physics=4 with use_mp_re=1, but in WRF_dm re_cloud never reached
    radiation (it stayed identically zero), while in WRF_dm_v2 it is populated.
    Radiation therefore sees different droplet radii between the trees, worth
    1.9-5.5% in SW CRE (org -1.9%, org_3aer -3.1%, mid_3aer -3.7%, mid -5.5%,
    ARG v2 relative to ARG over their shared frames).  SW CRE is not plotted -- panel (i) is RWP -- but it is
    still in the summary table, so that column is NOT interchangeable with v3's.
    Every other field is.

    The four ARG runs do not all end at the same time: the campaign gave the
    two org cases a 13:00 target and the two mid cases 14:02, so ARG org stops
    at 13:04 and ARG org_3aer at 13:20 while the mid pair runs to 14:04.  (The
    mid_3aer run also wrote a 14:06 wrfout, but the campaign cancelled the job
    mid-write and left it truncated at 9.6 of 14.7 GB; it fails to read and is
    skipped, which is why that case ends at 14:04 like the rest.)
    The figure is drawn WITHOUT --tend-common by default so the two short
    curves simply end early rather than clipping the other six back with them.

    The cache is incremental per case and is seeded from the v4 cache: the v4
    run already read these four WRF_dm_v2 cases (under its V2_* keys), so
    nothing is re-read off Lustre unless a case is missing or --refresh names
    it.  One wrfout here is ~14 GB, so a full re-read is hours.

USAGE
  python plot_arg_vs_ml_consequences_rwp_v5.py --workers 8 --font-scale 1.6
  python plot_arg_vs_ml_consequences_rwp_v5.py --ccn CCN4      # S=0.2% instead
  # re-read a case whose run has grown since the cache was written:
  python plot_arg_vs_ml_consequences_rwp_v5.py --refresh ARG_mid,ARG_mid_3aer
"""

import os
import re
import glob
import argparse
import textwrap
import datetime as dt
from multiprocessing import Pool

import numpy as np
import netCDF4 as nc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# ----------------------------------------------------------------------------
ML_DIR = ("/scratch/07088/tg863871/Perlm_Backup/"
          "WRF_stam3_ml_mixout_bce_lr1e-4_cosine_more_fix/test")
# The ARG-ACT record.  v3 and earlier read this from WRF_dm; from v5 it is the
# WRF_dm_v2 re-run of the same four cases from the same 09:02 restart (see the
# module docstring -- identical physics, Registry-only difference, except that
# radiation now sees non-zero re_cloud).
ARG_DIR = "/scratch/07088/tg863871/WRF_dm_v2/test"

case_dirs = {
    "org":          os.path.join(ML_DIR,  "20170714_1aer_org"),
    "org_3aer":     os.path.join(ML_DIR,  "20170714_3aer_org"),
    "mid":          os.path.join(ML_DIR,  "20170714_1aer_mid"),
    "mid_3aer":     os.path.join(ML_DIR,  "20170714_3aer_mid"),
    "ARG_org":      os.path.join(ARG_DIR, "20170714_org_1aer_ARG"),
    "ARG_org_3aer": os.path.join(ARG_DIR, "20170714_org_3aer_ARG"),
    "ARG_mid":      os.path.join(ARG_DIR, "20170714_mid_1aer_ARG"),
    "ARG_mid_3aer": os.path.join(ARG_DIR, "20170714_mid_3aer_ARG"),
}
# The v4 cache stored these same four directories under V2_* keys.  Seeding
# translates them onto the ARG_* keys used here; the v4 cache's own ARG_* keys
# (the WRF_dm tree) are deliberately NOT inherited.
SEED_RENAME = {"ARG_org":      "V2_org",
               "ARG_org_3aer": "V2_org_3aer",
               "ARG_mid":      "V2_mid",
               "ARG_mid_3aer": "V2_mid_3aer"}

# ARG-ACT runs, i.e. the ones whose wrfout carries no CCN and for which CCN has
# to be read from the matching wrfrst.  True for BOTH ARG trees.
IS_ARG = {c: c.startswith("ARG_") for c in case_dirs}

CASE_STYLE = {
    "org":          dict(color="#1f77b4", ls="-",  label="org (ML)"),
    "ARG_org":      dict(color="#1f77b4", ls="--", label="org (ARG)"),
    "org_3aer":     dict(color="#d62728", ls="-",  label="org_3aer (ML)"),
    "ARG_org_3aer": dict(color="#d62728", ls="--", label="org_3aer (ARG)"),
    "mid":          dict(color="#2ca02c", ls="-",  label="mid (ML)"),
    "ARG_mid":      dict(color="#2ca02c", ls="--", label="mid (ARG)"),
    "mid_3aer":     dict(color="#9467bd", ls="-",  label="mid_3aer (ML)"),
    "ARG_mid_3aer": dict(color="#9467bd", ls="--", label="mid_3aer (ARG)"),
}
CASE_LW = {c: 2.0 for c in case_dirs}

Z_TOP = 2000.0
QC_THRESH = 1.0e-6
RHO_W = 1000.0
G, RD, CP, P0 = 9.81, 287.06, 1004.6, 1.0e5
MIN_AGE = 60
CCN_VAR = "CCN2"
CCN_SS = {"CCN1": "0.02%", "CCN2": "0.05%", "CCN3": "0.1%",
          "CCN4": "0.2%",  "CCN5": "0.5%",  "CCN6": "1.0%"}

# ---------------------------------------------------------------------------
# Known-bad output times, excluded from every panel and from the summary.
#
# ARG_mid_3aer 09:52-10:00 : one restart segment of that run, written in a
# single burst on 2026-07-24, sits on a flattened trajectory that catches up in
# a single step at 10:02.  With those five frames included the ARG/ML
# accumulation-number ratio balloons 1.006 -> 1.018 -> 1.030 -> 1.041 -> 1.052
# -> 1.063 across them and then falls straight back to 1.025, which no physical
# divergence does.  A re-run from wrfrst_09:50 was attempted and hung at
# 09:51, so the frames cannot currently be regenerated; they are dropped
# instead.  Same list as plot_aerosol_number_timeseries_clean.py.
RUN_DATE = "2017-07-15"
# Those five frames belong to the OLD WRF_dm run.  WRF_dm_v2 re-integrated
# that stretch normally -- its ac0 ratio walks smoothly through 09:52-10:00
# with no catch-up step -- so nothing is excluded here.  The entry is kept,
# empty, as the place to add times if a v2 segment ever needs the same
# treatment.
EXCLUDE_TIMES = {}

WRFOUT_GLOB = "wrfout_d01_*"
TIME_RE = re.compile(r"wrfout_d\d\d_(\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:\d{2})$")

# raw per-file means that get cached
SERIES = ["ac0", "nu0", "ccn", "QNDROP", "qc_kg", "nd_kg",
          "LWP_TOT", "TAU_QC_TOT", "RAIN_ACC", "SWCRE", "CLDCOL", "RWP"]

MOIST = ("mid", "mid_3aer", "ARG_mid", "ARG_mid_3aer")
DRY   = ("org", "org_3aer", "ARG_org", "ARG_org_3aer")

# (variable, title, ylabel, may-be-log, case subset or None for all)
#
# v6 drops the surface precipitation rate entirely and gives rain water path
# TWO panels instead of one, split by background.  On a single shared axis the
# moist RWP (~170-270 g m-2) and the default-humidity RWP (~0.06-5 g m-2) are
# four decades apart; even on a log axis that compresses the org curves --
# whose collapse and slow recovery is the whole point of Section 4 -- into the
# bottom decade.  Splitting them lets each pair set its own scale.
PANELS = [
    ("ac0",        "(a) accumulation-mode number",      r"$N_{acc}$ [cm$^{-3}$]",     True,  None),
    ("nu0",        "(b) Aitken-mode number",            r"$N_{ait}$ [cm$^{-3}$]",     True,  None),
    ("ccn",        "(c) CCN",                           r"CCN [cm$^{-3}$]",           True,  None),
    ("QNDROP",     "(d) cloud droplet number",          r"$N_d$ [cm$^{-3}$]",         True,  None),
    ("r_v",        "(e) droplet mean-volume radius",    r"$r_v$ [$\mu$m]",            False, None),
    ("LWP_TOT",    "(f) liquid water path",             r"LWP [g m$^{-2}$]",          False, None),
    ("RWP",        "(g) rain water path, mid",          r"RWP [g m$^{-2}$]",          False, MOIST),
    ("rain_rate",  "(h) precipitation rate, mid",       r"[mm h$^{-1}$]",             False, MOIST),
    ("RWP",        "(i) rain water path, org",          r"RWP [g m$^{-2}$]",          True,  DRY),
    ("rain_rate",  "(j) precipitation rate, org",       r"[mm h$^{-1}$]",             True,  DRY),
]



def parse_time(s, default_date="2017-07-15"):
    if s is None:
        return None
    s = s.strip().replace("_", " ")
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%H:%M:%S", "%H:%M"):
        try:
            t = dt.datetime.strptime(s, fmt)
        except ValueError:
            continue
        if fmt.startswith("%H"):
            dd = dt.datetime.strptime(default_date, "%Y-%m-%d")
            t = t.replace(year=dd.year, month=dd.month, day=dd.day)
        return t
    raise argparse.ArgumentTypeError("cannot parse time %r" % s)


def _read_ccn_from_restart(path, ccn_var):
    """ARG wrfout has no CCN; the matching wrfrst does.  Returns the 3-D field
    or None."""
    rst = path.replace("wrfout_d01_", "wrfrst_d01_")
    if not os.path.exists(rst):
        return None
    try:
        with nc.Dataset(rst) as r:
            if ccn_var not in r.variables:
                return None
            return np.asarray(r.variables[ccn_var][0], dtype=np.float32)
    except Exception:
        return None


def process_file(args):
    path, is_arg, ccn_var = args
    out = {k: np.nan for k in SERIES}
    try:
        with nc.Dataset(path) as d:
            zf = (np.asarray(d.variables["PH"][0], dtype=np.float32) +
                  np.asarray(d.variables["PHB"][0], dtype=np.float32)) / G
            zf -= np.asarray(d.variables["HGT"][0], dtype=np.float32)[None, :, :]
            zh = 0.5 * (zf[:-1] + zf[1:])
            dz = np.diff(zf, axis=0)
            del zf

            p = (np.asarray(d.variables["P"][0], dtype=np.float32) +
                 np.asarray(d.variables["PB"][0], dtype=np.float32))
            theta = np.asarray(d.variables["T"][0], dtype=np.float32) + 300.0
            tk = theta * (p / P0) ** (RD / CP)
            qv = np.asarray(d.variables["QVAPOR"][0], dtype=np.float32)
            rho = p / (RD * tk * (1.0 + 0.61 * qv))
            del p, theta, tk, qv

            # air mass per column-layer.  w_full keeps the WHOLE column for the
            # rain-water path; w is then clipped to Z_TOP for the in-cloud
            # scalars, as in every earlier version.
            w_full = np.asarray(dz * rho, dtype=np.float32)
            w = w_full.copy()
            low = zh <= Z_TOP
            w[~low] = 0.0
            del dz, zh

            qc = np.asarray(d.variables["QCLOUD"][0], dtype=np.float32)
            cloudy = (qc > QC_THRESH) & low
            del low
            cloud_col = cloudy.any(axis=0)
            out["CLDCOL"] = float(cloud_col.mean())
            ncol = float(cloud_col.sum())

            # ---- domain-mean 2-D fields (ALL columns) ----
            rain = np.zeros(cloud_col.shape, dtype=np.float64)
            for v in ("RAINNC", "RAINC"):
                if v in d.variables:
                    rain += np.asarray(d.variables[v][0], dtype=np.float64)
            out["RAIN_ACC"] = float(rain.mean())
            if "SWUPT" in d.variables and "SWUPTC" in d.variables:
                su = np.asarray(d.variables["SWUPT"][0], dtype=np.float64)
                sc = np.asarray(d.variables["SWUPTC"][0], dtype=np.float64)
                out["SWCRE"] = float(-(su - sc).mean())

            # ---- 2-D cloud fields (cloudy columns) ----
            if ncol > 0:
                for v in ("LWP_TOT", "TAU_QC_TOT"):
                    if v in d.variables:
                        x = np.asarray(d.variables[v][0], dtype=np.float64)
                        out[v] = float(x[cloud_col].mean())

                # rain water path: full-column integral of q_rain, then the same
                # cloudy-column area mean LWP_TOT gets.  kg m-2 -> g m-2.
                if "QRAIN" in d.variables:
                    qr = np.asarray(d.variables["QRAIN"][0], dtype=np.float32)
                    rwp2d = (qr * w_full).sum(axis=0, dtype=np.float64)
                    out["RWP"] = float(rwp2d[cloud_col].mean()) * 1000.0
                    del qr, rwp2d
            del w_full

            # ---- 3-D fields, air-mass weighted over cloudy columns ----
            w[:, ~cloud_col] = 0.0
            del cloudy, cloud_col
            wsum = w.sum(dtype=np.float64)
            if wsum > 0:
                for v in ("ac0", "nu0", "QNDROP"):
                    if v not in d.variables:
                        continue
                    x = np.asarray(d.variables[v][0], dtype=np.float32)
                    out[v] = float(((x * rho * 1.0e-6) * w).sum(dtype=np.float64) / wsum)
                    if v == "QNDROP":                       # keep per-kg for r_v
                        out["nd_kg"] = float((x * w).sum(dtype=np.float64) / wsum)
                    del x
                out["qc_kg"] = float((qc * w).sum(dtype=np.float64) / wsum)

                if ccn_var in d.variables:                  # ML runs
                    x = np.asarray(d.variables[ccn_var][0], dtype=np.float32)
                    out["ccn"] = float((x * w).sum(dtype=np.float64) / wsum)
                elif is_arg:                                # ARG runs -> restart
                    x = _read_ccn_from_restart(path, ccn_var)
                    if x is not None and x.shape == w.shape:
                        out["ccn"] = float((x * w).sum(dtype=np.float64) / wsum)
                    del x
            del w, rho, qc

            tstr = d.variables["Times"][0].tobytes().decode().strip("\x00").strip()
    except Exception as exc:                                # noqa: BLE001
        print("  !! failed %s: %s" % (os.path.basename(path), exc), flush=True)
        return None

    if not tstr:
        m = TIME_RE.search(path)
        tstr = m.group(1) if m else None
    return dt.datetime.strptime(tstr, "%Y-%m-%d_%H:%M:%S"), out


def compute_case(case, cdir, pool, ccn_var, stride=1, min_age=MIN_AGE):
    files = sorted(glob.glob(os.path.join(cdir, WRFOUT_GLOB)))
    files = [f for f in files if not f.endswith(".anomalous")]
    if min_age > 0:
        now = __import__("time").time()
        fresh = [f for f in files if now - os.path.getmtime(f) < min_age]
        if fresh:
            print("[%s] skipping %d file(s) still being written" % (case, len(fresh)),
                  flush=True)
        files = [f for f in files if now - os.path.getmtime(f) >= min_age]
    files = files[::stride]
    if not files:
        print("[%s] no usable wrfout in %s" % (case, cdir), flush=True)
        return None
    print("[%s] %d files" % (case, len(files)), flush=True)

    tasks = [(f, IS_ARG[case], ccn_var) for f in files]
    res = [r for r in pool.imap_unordered(process_file, tasks, chunksize=1)
           if r is not None]
    if not res:
        return None
    res.sort(key=lambda r: r[0])

    times = np.array([r[0] for r in res])
    s = {k: np.array([r[1].get(k, np.nan) for r in res]) for k in SERIES}

    # ---- derived ----
    with np.errstate(divide="ignore", invalid="ignore"):
        vol = s["qc_kg"] / (s["nd_kg"] * RHO_W)             # m3 per droplet
        # v7 BUG FIX: through v6 this was cbrt(vol) -- the cube root of the
        # volume per droplet with NO sphere factor -- so every plotted radius
        # was (4 pi/3)^(1/3) = 1.612x the mean-volume radius the docstring
        # specifies.  (4/3) pi r^3 = vol  =>  r = (3 vol / 4 pi)^(1/3).
        s["r_v"] = np.where(s["nd_kg"] > 0,
                            np.cbrt(np.maximum(3.0 * vol / (4.0 * np.pi), 0.0))
                            * 1.0e6, np.nan)
    dt_h = np.diff(times).astype("timedelta64[s]").astype(float) / 3600.0
    rate = np.full(times.shape, np.nan)
    with np.errstate(invalid="ignore"):
        rate[1:] = np.diff(s["RAIN_ACC"]) / np.where(dt_h > 0, dt_h, np.nan)
    rate[rate < 0] = np.nan                                  # restart resets, if any
    s["rain_rate"] = rate
    s["LWP_TOT"] = s["LWP_TOT"] * 1000.0                     # kg m-2 -> g m-2
    s["qc_gkg"]  = s["qc_kg"] * 1000.0                       # kg/kg -> g/kg
    return times, s


def _vislen(s):
    """Rendered length of a mathtext label: $, braces, ^ and backslashes cost
    no width on the page, so len() badly overestimates e.g. r'$N_{acc}$ [cm$^{-3}$]'."""
    s = s.replace(r"\mu", "u")
    return len(re.sub(r"[${}^\\]", "", s))


def make_plot(data, outstem, ccn_var, dpi=200, logy=True,
              tstart=None, tend=None, cases=None, fscale=1.0,
              hide_title_i=False):
    """fscale multiplies every font size.

    At fscale > 1 the figure WIDTH stays 15 in -- that is what fixes the
    apparent text size once the figure is scaled to a column -- and everything
    else is laid out in inches around the enlarged text: titles and long
    y-labels are wrapped, the margins/legend band are reserved explicitly
    (tight_layout stops filling its rect once the text is that large), and the
    figure HEIGHT grows so the panels keep a usable aspect instead of being
    squeezed to ~1 in.  fscale == 1 reproduces the original layout exactly.
    """
    figw = 15.0
    fs_title, fs_lab = 16.0 * fscale, 15.0 * fscale
    fs_tick, fs_leg = 13.0 * fscale, 15.0 * fscale
    fs_sup = 17.0 * fscale
    scaled = abs(fscale - 1.0) > 1.0e-9          # 1x keeps the original layout
    pt = 1.0 / 72.0

    plotted = [c for c in (cases or case_dirs) if c in data]
    nleg = len(plotted)

    sup = ""          # v7: figure title removed (requested)
    # Column count is driven by the curve count: 8 curves stack 4-and-4 in two
    # columns, 12 need a third or the legend block grows taller than a panel.
    # A scaled column is ~4.6 in wide at fscale 1.6, so three still fit in 15 in
    # but four do not.
    ncol_leg = 4 if not scaled else (3 if nleg > 8 else 2)
    titles = [("(c) CCN at S=%s" % CCN_SS.get(ccn_var, "?")) if v == "ccn"
              else t for v, t, _, _, _ in PANELS]

    if not scaled:
        figh = 20.0
        ylabs = [y for _, _, y, _, _ in PANELS]
    else:
        # --- widths: fixed 15 in, split between two panels and their gutters
        gap_w = 2 * 1.4 * fs_lab * pt + 5 * 0.6 * fs_tick * pt
        left, right = gap_w + 0.15, 0.15
        axw = (figw - left - right - gap_w) / 2.0
        titles = ["\n".join(textwrap.wrap(t, max(12, int(axw / (0.52 * fs_title * pt)))))
                  for t in titles]
        ntitle = max(t.count("\n") + 1 for t in titles)

        # --- heights: derived, so the panels keep roughly their 1x aspect and
        #     a rotated y-label is never taller than the panel it labels
        ylab_cols = 14
        ylabs = ["\n".join(textwrap.wrap(y, ylab_cols)) if _vislen(y) > ylab_cols
                 else y for _, _, y, _, _ in PANELS]
        axh = max(0.465 * axw, ylab_cols * 0.52 * fs_lab * pt)
        sup_lines = []
        sup_h = 0.10
        row_h = ntitle * 1.3 * fs_title * pt + 0.12           # per-panel title
        bot_h = (-(-nleg // ncol_leg) * 1.7 * fs_leg * pt     # legend block
                 + 2 * 1.35 * fs_tick * pt                    # 2-line x labels
                 + 1.9 * fs_lab * pt + 0.25)                  # "time (UTC)"
        figh = sup_h + bot_h + 5 * (row_h + axh)

    fig, axes = plt.subplots(5, 2, figsize=(figw, figh), sharex=True)
    axes = axes.ravel()

    for ax, title, ylab, (v, _, _, canlog, subset) in zip(axes, titles, ylabs, PANELS):
        for case in (cases or case_dirs):
            if case not in data:
                continue
            if subset is not None and case not in subset:
                continue
            times, s = data[case]
            if v not in s:
                continue
            y = np.asarray(s[v], dtype=float)
            t = np.asarray(times)
            keep = np.isfinite(y)
            if tstart is not None:
                keep &= t >= tstart
            if tend is not None:
                keep &= t <= tend
            if not keep.any():
                continue
            st = CASE_STYLE[case]
            ax.plot(t[keep], y[keep], color=st["color"], ls=st["ls"],
                    lw=CASE_LW.get(case, 2.0), label=st["label"])
        if hide_title_i and title.lstrip().startswith("(i)"):
            ax.set_title("")
        else:
            ax.set_title(title, fontsize=fs_title)
        ax.set_ylabel(ylab, fontsize=fs_lab)
        ax.tick_params(labelsize=fs_tick)
        if logy and canlog:
            ax.set_yscale("log")
            ax.grid(alpha=0.3, lw=0.5, which="major")
            ax.grid(alpha=0.15, lw=0.4, which="minor")
        else:
            ax.grid(alpha=0.3, lw=0.5)

    for ax in axes[8:]:
        ax.set_xlabel("time (UTC)", fontsize=fs_lab)
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d\n%H:%M"))
        if scaled:
            # An explicit hour locator rather than AutoDateLocator: at 2x a tick
            # label is ~1 in wide, the tick budget drops to 4-5, and
            # AutoDateLocator answers a budget that small by giving up and
            # falling back to a 30-minute default -- denser than what it
            # replaced.
            span_h = (ax.get_xlim()[1] - ax.get_xlim()[0]) * 24.0
            nmax = max(3, int(round(axw / (0.52 * 5.5 * fs_tick * pt))))
            step = next((h for h in (1, 2, 3, 6, 12, 24)
                         if span_h / h <= nmax), 24)
            ax.xaxis.set_major_locator(mdates.HourLocator(interval=step))
        else:
            ax.xaxis.set_major_locator(mdates.AutoDateLocator(maxticks=6))
        ax.tick_params(labelbottom=True)

    handles, labels = axes[0].get_legend_handles_labels()
    if not scaled:
        fig.legend(handles, labels, loc="lower center", ncol=ncol_leg,
                   frameon=False, fontsize=fs_leg, bbox_to_anchor=(0.5, -0.004))
        fig.tight_layout(rect=(0, 0.035, 1, 0.965))
    else:
        fig.subplots_adjust(left=left / figw, right=1.0 - right / figw,
                            top=1.0 - (sup_h + row_h) / figh,
                            bottom=bot_h / figh,
                            wspace=gap_w / axw, hspace=row_h / axh)
        # legend sits inside the reserved bottom band, below the x labels
        fig.legend(handles, labels, loc="lower center", ncol=ncol_leg,
                   frameon=False, fontsize=fs_leg,
                   bbox_to_anchor=(0.5, 0.15 / figh))
    for ext, kw in (("png", dict(dpi=dpi)), ("eps", dict(format="eps"))):
        path = "%s.%s" % (outstem, ext)
        fig.savefig(path, bbox_inches="tight", **kw)
        print("wrote %s" % path, flush=True)
    plt.close(fig)


def main():
    global CCN_VAR
    ap = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    ap.add_argument("--workers", type=int,
                    default=int(os.environ.get("SLURM_CPUS_ON_NODE",
                                               os.cpu_count() or 8)))
    ap.add_argument("--out", default="arg_vs_ml_consequences_rwp_v5")
    ap.add_argument("--ccn", default="CCN2", choices=sorted(CCN_SS))
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--tstart", default=None)
    ap.add_argument("--tend", default=None)
    ap.add_argument("--tend-common", action="store_true")
    ap.add_argument("--hide-title-i", action="store_true",
                    help="v6 variant B: suppress the panel-(i) title, so the "
                         "rain-water-path heading appears once, on panel (g)")
    ap.add_argument("--font-scale", type=float, default=1.0,
                    help="multiply every font size (2 = twice as big)")
    ap.add_argument("--cache", default=None)
    ap.add_argument("--seed-cache",
                    default="arg_vs_ml_consequences_rwp_v4_series_%s.npz",
                    help="cache from an earlier run to inherit cases from when "
                         "--cache does not exist yet ('%%s' -> the CCN name)")
    ap.add_argument("--refresh", default="", metavar="CASES",
                    help="comma-separated cases to re-read even if cached, or "
                         "'all'")
    args = ap.parse_args()
    CCN_VAR = args.ccn

    tstart, tend = parse_time(args.tstart), parse_time(args.tend)
    cache = args.cache or ("%s_series_%s.npz" % (args.out, args.ccn))
    schema = dict(series=tuple(SERIES), z_top=Z_TOP, qc=QC_THRESH, ccn=args.ccn)

    def _load(path):
        if not os.path.exists(path):
            return None
        print("loading cache %s" % path, flush=True)
        z = np.load(path, allow_pickle=True)
        if z["schema"].item() != schema:
            print("  schema mismatch, ignoring", flush=True)
            return None
        return z["data"].item()

    # Incremental by case.  Adding the WRF_dm_v2 runs must not force a re-read
    # of the eight cases that are already cached: one wrfout here is ~14 GB, so
    # a full re-read is ~17 TB off Lustre and takes hours.  Cases absent from
    # the cache -- and only those -- are computed and folded back in.
    seed = args.seed_cache
    if seed and "%s" in seed:
        seed = seed % args.ccn
    data = _load(cache)
    if data is None and seed and os.path.abspath(seed) != os.path.abspath(cache):
        raw = _load(seed)
        if raw is not None:
            # Take the ML cases under their own names and the four WRF_dm_v2
            # cases from the V2_* keys the v4 script wrote them under.  Nothing
            # else is inherited -- in particular not the v4 cache's ARG_* keys,
            # which hold the OLD WRF_dm tree and would silently reinstate the
            # very data this version replaces.
            data = {c: raw[c] for c in case_dirs
                    if c not in SEED_RENAME and c in raw}
            data.update({new: raw[old] for new, old in SEED_RENAME.items()
                         if old in raw})
            print("  seeded %d case(s) from %s (%s remapped from %s)"
                  % (len(data), seed,
                     ", ".join(sorted(k for k in SEED_RENAME if k in data)),
                     ", ".join(sorted(SEED_RENAME[k] for k in SEED_RENAME
                                      if k in data))), flush=True)
    data = dict(data or {})

    stale = ([c for c in case_dirs] if args.refresh.strip() == "all"
             else [c.strip() for c in args.refresh.split(",") if c.strip()])
    for c in stale:
        if c not in case_dirs:
            raise SystemExit("--refresh: unknown case %r" % c)
        data.pop(c, None)

    missing = [c for c in case_dirs if c not in data]
    if missing:
        print("computing %d case(s): %s" % (len(missing), ", ".join(missing)),
              flush=True)
        print("using %d workers" % args.workers, flush=True)
        with Pool(args.workers) as pool:
            for case in missing:
                r = compute_case(case, case_dirs[case], pool, args.ccn,
                                 stride=args.stride)
                if r is not None:
                    data[case] = r
        np.savez(cache, data=np.array(data, dtype=object),
                 schema=np.array(schema, dtype=object))
        print("cached -> %s" % cache, flush=True)
    else:
        print("all %d case(s) already cached" % len(data), flush=True)

    if not data:
        raise SystemExit("no data")

    # Drop known-bad times.  Applied to the loaded series, not at read time, so
    # the cache stays valid and the exclusion list can change without a re-read.
    for case, tstrs in EXCLUDE_TIMES.items():
        if case not in data:
            continue
        drop = {dt.datetime.strptime(f"{RUN_DATE}_{t}:00", "%Y-%m-%d_%H:%M:%S")
                for t in tstrs}
        times, ser = data[case]
        keep = np.array([t not in drop for t in times], dtype=bool)
        n_out = int((~keep).sum())
        if n_out:
            print("excluding %d known-bad time(s) from %s: %s"
                  % (n_out, case, ", ".join(sorted(tstrs))), flush=True)
            data[case] = (np.asarray(times)[keep],
                          {k: np.asarray(v)[keep] for k, v in ser.items()})

    cases = [c for c in case_dirs if c in data]
    print("plotting: %s" % ", ".join(cases), flush=True)
    if args.tend_common:
        # OFF by default here.  The two ARG org runs were only ever targeted at
        # 13:00, so a common end clips the mid pair and the ML runs back an hour
        # to match them.  Curves are left to end where their runs end instead.
        common = min(data[c][0][-1] for c in cases)
        tend = common if tend is None else min(tend, common)
        print("common end: %s (over %s)" % (common, ", ".join(cases)), flush=True)

    make_plot(data, args.out, args.ccn, logy=True,
              tstart=tstart, tend=tend, cases=cases, fscale=args.font_scale,
              hide_title_i=args.hide_title_i)

    _summary(data, "ML vs ARG (WRF_dm_v2)", "", "ARG_", "ML", "ARG")


def _summary(data, heading, pre_a, pre_b, lab_a, lab_b):
    """Mean of each series over the times the two runs share."""
    print("\n=== %s, mean over the overlapping window ===" % heading)
    print(f"{'pair':12} {'var':12} {lab_a:>12} {lab_b:>12} "
          f"{lab_b + '/' + lab_a:>12}  (matched times)")
    for pair in ("org", "org_3aer", "mid", "mid_3aer"):
        a, b = pre_a + pair, pre_b + pair
        if a not in data or b not in data:
            continue
        ta, sa = data[a]; tb, sb = data[b]
        # Compare on MATCHED times only.  Masking each series independently to
        # a common end time is wrong when one of them has a gap: ARG_mid_3aer is
        # missing 09:52-10:00, so an unmatched mean skipped five declining
        # points on the ARG side that the ML side still had, and reported a
        # spurious 5% aerosol divergence where the matched ratio is 1.004.
        idx = {t: j for j, t in enumerate(tb)}
        pairs = [(i, idx[t]) for i, t in enumerate(ta) if t in idx]
        if not pairs:
            continue
        ma = [q[0] for q in pairs]; mb = [q[1] for q in pairs]
        print(f"-- {pair} ({len(pairs)} matched frames, "
              f"{ta[ma[0]]:%H:%M} - {ta[ma[-1]]:%H:%M})")
        for v in ("ac0", "ccn", "QNDROP", "r_v", "rain_rate", "LWP_TOT",
                  "RWP", "qc_gkg", "SWCRE"):
            with np.errstate(invalid="ignore"):
                xa = np.nanmean(np.asarray(sa[v], float)[ma])
                xb = np.nanmean(np.asarray(sb[v], float)[mb])
            if np.isfinite(xa) and np.isfinite(xb) and xa != 0:
                print(f"{pair:12} {v:12} {xa:>12.4g} {xb:>12.4g} {xb/xa:>12.3f}")


if __name__ == "__main__":
    main()
