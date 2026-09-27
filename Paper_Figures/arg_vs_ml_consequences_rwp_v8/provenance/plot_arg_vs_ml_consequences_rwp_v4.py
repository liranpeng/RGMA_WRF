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
    (WRF_dm) do NOT -- for those cases CCN is read from the matching wrfrst,
    which does carry CCN1..CCN6 (restart and history are written on the same
    2-minute cadence, one restart per history time).  LWP_TOT and TAU_QC_TOT
    are wrfout-only in both trees.  TAU_QC_TOT is written by both but is
    IDENTICALLY ZERO everywhere in the WRF_dm (ARG) build -- the optical-depth
    diagnostic is not computed there -- so it is not plotted; in-cloud cloud
    water is shown instead.

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

WHAT v4 ADDS OVER v3
    The four ARG cases were re-run in /scratch/.../WRF_dm_v2/test from the same
    09:02 restart.  That tree differs from WRF_dm only in Registry.EM_COMMON
    (re_cloud/re_ice/re_snow promoted to "rh" so the effective radii reach
    wrfout); activation, microphysics and the aerosol treatment are unchanged.
    They are plotted as a THIRD line style -- DOTTED, same colour as the case
    they repeat -- so the reproduction can be read off the figure directly: a
    dotted curve should lie on top of the dashed one of its colour.  A second
    summary table reports the ARG/ARG-v2 ratios over the frames both trees
    have.

    Two consequences of that campaign still being mid-flight:
      * the dotted curves are SHORT (they end wherever the runs have got to);
      * --tend-common ignores them when it picks the common end, otherwise a
        30-minute re-run would clip every other curve back to 09:30.

    The cache is now INCREMENTAL PER CASE.  Adding four cases no longer forces
    a re-read of the eight already cached -- one wrfout here is ~14 GB, so a
    full re-read is several hours of Lustre traffic.  Cases missing from the
    cache are computed and folded back in; --refresh forces named cases (or
    'all') to be re-read.

USAGE
  python plot_arg_vs_ml_consequences_rwp_v4.py --workers 12 --tend-common
  python plot_arg_vs_ml_consequences_rwp_v4.py --ccn CCN4      # S=0.2% instead
  # pick up frames the v2 campaign has written since the last figure:
  python plot_arg_vs_ml_consequences_rwp_v4.py --refresh V2_org,V2_org_3aer,V2_mid,V2_mid_3aer
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
ARG_DIR = "/scratch/07088/tg863871/WRF_dm/test"
# WRF_dm_v2: the SAME four ARG-ACT cases re-run from the same 09:02 restart.
# The tree differs from WRF_dm only in Registry.EM_COMMON (re_cloud/re_ice/
# re_snow promoted to "rh" so the effective radii reach wrfout); the physics
# and the activation scheme are identical, so these curves are a reproduction
# check on the ARG record, not a new experiment.  The campaign is still
# running -- as of writing it covers only the first ~30 min of the window --
# so the dotted curves are short by design and will lengthen on a re-run.
ARG_V2_DIR = "/scratch/07088/tg863871/WRF_dm_v2/test"

case_dirs = {
    "org":          os.path.join(ML_DIR,  "20170714_1aer_org"),
    "org_3aer":     os.path.join(ML_DIR,  "20170714_3aer_org"),
    "mid":          os.path.join(ML_DIR,  "20170714_1aer_mid"),
    "mid_3aer":     os.path.join(ML_DIR,  "20170714_3aer_mid"),
    "ARG_org":      os.path.join(ARG_DIR, "20170714_org_1aer_ARG"),
    "ARG_org_3aer": os.path.join(ARG_DIR, "20170714_org_3aer_ARG"),
    "ARG_mid":      os.path.join(ARG_DIR, "20170714_mid_1aer_ARG"),
    "ARG_mid_3aer": os.path.join(ARG_DIR, "20170714_mid_3aer_ARG"),
    "V2_org":       os.path.join(ARG_V2_DIR, "20170714_org_1aer_ARG"),
    "V2_org_3aer":  os.path.join(ARG_V2_DIR, "20170714_org_3aer_ARG"),
    "V2_mid":       os.path.join(ARG_V2_DIR, "20170714_mid_1aer_ARG"),
    "V2_mid_3aer":  os.path.join(ARG_V2_DIR, "20170714_mid_3aer_ARG"),
}
# ARG-ACT runs, i.e. the ones whose wrfout carries no CCN and for which CCN has
# to be read from the matching wrfrst.  True for BOTH ARG trees.
IS_ARG = {c: (c.startswith("ARG_") or c.startswith("V2_")) for c in case_dirs}
V2_CASES = [c for c in case_dirs if c.startswith("V2_")]

CASE_STYLE = {
    "org":          dict(color="#1f77b4", ls="-",  label="org (ML)"),
    "ARG_org":      dict(color="#1f77b4", ls="--", label="org (ARG)"),
    "V2_org":       dict(color="#1f77b4", ls=":",  label="org (ARG v2)"),
    "org_3aer":     dict(color="#d62728", ls="-",  label="org_3aer (ML)"),
    "ARG_org_3aer": dict(color="#d62728", ls="--", label="org_3aer (ARG)"),
    "V2_org_3aer":  dict(color="#d62728", ls=":",  label="org_3aer (ARG v2)"),
    "mid":          dict(color="#2ca02c", ls="-",  label="mid (ML)"),
    "ARG_mid":      dict(color="#2ca02c", ls="--", label="mid (ARG)"),
    "V2_mid":       dict(color="#2ca02c", ls=":",  label="mid (ARG v2)"),
    "mid_3aer":     dict(color="#9467bd", ls="-",  label="mid_3aer (ML)"),
    "ARG_mid_3aer": dict(color="#9467bd", ls="--", label="mid_3aer (ARG)"),
    "V2_mid_3aer":  dict(color="#9467bd", ls=":",  label="mid_3aer (ARG v2)"),
}
# The v2 runs reproduce their WRF_dm counterparts to ~0.1%, so a plain dotted
# line lands exactly under the dashed one and is invisible -- in panel (a) the
# ML and ARG curves already coincide, so three lines share the same pixels.
# Open circles every few frames survive that overlap: they read as markers
# ON the curve rather than as a fourth line, and the 16-frame v2 records give
# ~5 of them.  Drawn last (V2_* is last in case_dirs) so they sit on top.
CASE_LW = {c: (2.4 if c.startswith("V2_") else 2.0) for c in case_dirs}
V2_MARKER = dict(marker="o", markevery=3, markersize=8, markerfacecolor="none",
                 markeredgewidth=1.8)

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
EXCLUDE_TIMES = {
    "ARG_mid_3aer": ["09:52", "09:54", "09:56", "09:58", "10:00"],
}

WRFOUT_GLOB = "wrfout_d01_*"
TIME_RE = re.compile(r"wrfout_d\d\d_(\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:\d{2})$")

# raw per-file means that get cached
SERIES = ["ac0", "nu0", "ccn", "QNDROP", "qc_kg", "nd_kg",
          "LWP_TOT", "TAU_QC_TOT", "RAIN_ACC", "SWCRE", "CLDCOL", "RWP"]

PANELS = [
    ("ac0",        "(a) accumulation-mode number",      r"$N_{acc}$ [cm$^{-3}$]",     True),
    ("nu0",        "(b) Aitken-mode number",            r"$N_{ait}$ [cm$^{-3}$]",     True),
    ("ccn",        "(c) CCN",                           r"CCN [cm$^{-3}$]",           True),
    ("QNDROP",     "(d) cloud droplet number",          r"$N_d$ [cm$^{-3}$]",         True),
    ("r_v",        "(e) droplet mean-volume radius",    r"$r_v$ [$\mu$m]",            False),
    ("LWP_TOT",    "(f) liquid water path",             r"LWP [g m$^{-2}$]",          False),
    ("rain_rate",  "(g) surface precipitation rate",    r"[mm h$^{-1}$]",             True),
    ("qc_gkg",     "(h) in-cloud cloud water",          r"$q_c$ [g kg$^{-1}$]",       False),
    ("RWP",        "(i) rain water path",               r"RWP [g m$^{-2}$]",          True),
    ("CLDCOL",     "(j) cloudy-column fraction",        "fraction of domain [-]",     False),
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
        s["r_v"] = np.where(s["nd_kg"] > 0,
                            np.cbrt(np.maximum(vol, 0.0)) * 1.0e6, np.nan)
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
              tstart=None, tend=None, cases=None, fscale=1.0):
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
    has_v2 = any(c.startswith("V2_") for c in plotted)

    sup = ("Consequences of the ARG activation bias: ML-ACT (solid) vs "
           "ARG-ACT (dashed)%s|"
           "aerosol / cloud fields averaged over cloudy columns below "
           "%d m AGL; LWP and RWP are full-column integrals over cloudy "
           "columns; precipitation is a domain mean"
           % (" vs the WRF_dm_v2 ARG re-runs (dotted, open circles)"
              if has_v2 else "",
              int(Z_TOP)))
    # Column count is driven by the curve count: 8 curves stack 4-and-4 in two
    # columns, 12 need a third or the legend block grows taller than a panel.
    # A scaled column is ~4.6 in wide at fscale 1.6, so three still fit in 15 in
    # but four do not.
    ncol_leg = 4 if not scaled else (3 if nleg > 8 else 2)
    titles = [("(c) CCN at S=%s" % CCN_SS.get(ccn_var, "?")) if v == "ccn"
              else t for v, t, _, _ in PANELS]

    if not scaled:
        figh = 20.0
        ylabs = [y for _, _, y, _ in PANELS]
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
                 else y for _, _, y, _ in PANELS]
        axh = max(0.465 * axw, ylab_cols * 0.52 * fs_lab * pt)
        sup_lines = [w for part in sup.split("|") for w in
                     textwrap.wrap(part, max(30, int(figw / (0.52 * fs_sup * pt))))]
        sup_h = len(sup_lines) * 1.25 * fs_sup * pt + 0.20
        row_h = ntitle * 1.3 * fs_title * pt + 0.12           # per-panel title
        bot_h = (-(-nleg // ncol_leg) * 1.7 * fs_leg * pt     # legend block
                 + 2 * 1.35 * fs_tick * pt                    # 2-line x labels
                 + 1.9 * fs_lab * pt + 0.25)                  # "time (UTC)"
        figh = sup_h + bot_h + 5 * (row_h + axh)

    fig, axes = plt.subplots(5, 2, figsize=(figw, figh), sharex=True)
    axes = axes.ravel()

    for ax, title, ylab, (v, _, _, canlog) in zip(axes, titles, ylabs, PANELS):
        for case in (cases or case_dirs):
            if case not in data:
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
            kw = dict(V2_MARKER) if case.startswith("V2_") else {}
            ax.plot(t[keep], y[keep], color=st["color"], ls=st["ls"],
                    lw=CASE_LW.get(case, 2.0), label=st["label"], **kw)
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
        fig.suptitle(sup.replace("|", "\n"), fontsize=fs_sup)
        fig.tight_layout(rect=(0, 0.035, 1, 0.965))
    else:
        fig.subplots_adjust(left=left / figw, right=1.0 - right / figw,
                            top=1.0 - (sup_h + row_h) / figh,
                            bottom=bot_h / figh,
                            wspace=gap_w / axw, hspace=row_h / axh)
        fig.suptitle("\n".join(sup_lines), fontsize=fs_sup,
                     y=1.0 - 0.10 / figh, verticalalignment="top")
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
    ap.add_argument("--out", default="arg_vs_ml_consequences_rwp_v4")
    ap.add_argument("--ccn", default="CCN2", choices=sorted(CCN_SS))
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--tstart", default=None)
    ap.add_argument("--tend", default=None)
    ap.add_argument("--tend-common", action="store_true")
    ap.add_argument("--font-scale", type=float, default=1.0,
                    help="multiply every font size (2 = twice as big)")
    ap.add_argument("--cache", default=None)
    ap.add_argument("--seed-cache", default="arg_vs_ml_consequences_rwp_series_%s.npz",
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
        data = _load(seed)
        if data is not None:
            print("  seeded %d case(s) from %s" % (len(data), seed), flush=True)
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
        # The WRF_dm_v2 campaign is mid-flight and only a few frames deep, so
        # it is deliberately NOT part of the common end: letting it vote would
        # clip every other curve back to ~09:30 and throw away the whole
        # figure.  The dotted curves simply stop where the runs currently are.
        vote = [c for c in cases if c not in V2_CASES] or cases
        common = min(data[c][0][-1] for c in vote)
        tend = common if tend is None else min(tend, common)
        print("common end: %s (over %s)" % (common, ", ".join(vote)), flush=True)

    make_plot(data, args.out, args.ccn, logy=True,
              tstart=tstart, tend=tend, cases=cases, fscale=args.font_scale)

    _summary(data, "ML vs ARG", "", "ARG_", "ML", "ARG")
    # The v2 tree re-runs the same cases with the same physics, so this second
    # table is a reproduction check, not a comparison: over the frames both
    # trees have, every ratio should sit at 1.  Precipitation RATE is the one
    # exception worth ignoring at the first frame -- it is a time difference,
    # so it is NaN there and the mean rests on very few points.
    if any(c in data for c in V2_CASES):
        _summary(data, "ARG (WRF_dm) vs ARG v2 (WRF_dm_v2)", "ARG_", "V2_",
                 "ARG", "ARG v2")


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
