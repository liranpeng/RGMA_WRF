#!/usr/bin/env python3
"""
trace_east_boundary_subdomains_lwp_v22.py

Changes from v21:
  25. All output goes to ./dcre_figures_v22 (override with WRF_OUTDIR), so the
      v21 output in ./dcre_figures_more_fix* is preserved.  Everything v21
      wrote is written again here from the current model record.

  26. NEW Figure 10: the ML-vs-ARG domain comparison, drawn once PER HOUR plus
      once over the whole shared window.  Layout follows the two-panel
      pixel-wise style of replot_dcre_ml_arg_comparison.py rather than the 2x2
      of v21's Figure 8: (a) org background, (b) mid background, solid = ML-ACT,
      hatched = ARG-ACT, and the four ARG sub-terms greyed because TAU_QC_TOT
      and RE_QC are identically zero in the WRF_dm build, so they are the
      analytic 1/3:2/3 fallback and not a measurement.  Files:
        dCRE_ML_ARG_domain_comparison_HH<hh>.png   (one per hour)
        dCRE_ML_ARG_domain_comparison_ALL.png      (all shared times)
      An hour is kept only if BOTH the ML and the ARG pair of that background
      have at least MIN_TIMES_PER_HOUR samples in it; error bars on a two- or
      three-sample hour are not worth drawing.

  27. NEW Figure 11: half-hourly stacked-bar TIME SERIES of the domain-mean
      dCRE decomposition, one figure per pair.  One stacked bar per 30-minute
      bin, x-axis is time; components stack up/down about zero by sign and the
      black diamonds are DCRE_total +/- SE within the bin.  This is the
      met-binned Figure 6 stacking applied to a time axis instead of a met
      regime.  Files: dCRE_timeseries_stacked_<pair>.png
      Bars are drawn from the DOMAIN decomposition (dcre_pixelwise_<pair>
      _domain.csv), so every output time contributes once.

  28. The per-timestep map figures (Figures 1 and 2) are OFF by default -- they
      are ~2.8 MB each, two per output time, and dominate the runtime.  Set
      WRF_MAP_FIGURES=1 to restore them.

  29. The analysis window now ends at 2017-07-15 12:58 (TIME_END_DEFAULT,
      override with WRF_TEND; WRF_TEND="" keeps the whole record).  The cap is
      applied to the shared-time set before anything is traced, so the stacked
      time series, the time-mean bars and the totals all stop at the same
      instant instead of each ending wherever its own case record does.

Changes from v20:
  24. [FIX] mean-based decomposition: all box/domain averages are now taken
      over ONE joint valid mask -- points inside the box that are sunlit and
      non-NaN in BOTH the control and perturbed run.  Previously CF, CRE and
      the LW/SFC differences were each reduced with np.nanmean independently,
      so the control and perturbed means could be drawn from different pixel
      sets and DCRE_total = CRE_p - CRE_c mixed two samples.  In-cloud scalars
      remain per-run (the cloudy sets legitimately differ) but are now also
      restricted to finite points within that mask.  The pixel-wise routine
      already used joint masks (`any_cld`, `both`) and is unchanged.

  22. ML (non-ARG) cases now read the cosine_more_FIX tree:
        BASE_DIR = .../WRF_stam3_ml_mixout_bce_lr1e-4_cosine_more_fix/test
      v20 pointed at .../cosine (no "_more"), which predates both the
      cosine_more runs and the fn13/fn23 fix.  The _more_fix tree is the one
      in which module_mixactivate.F writes the diagnostic ARG activated
      fraction at BOTH activate_ml call sites (4 BUGFIX markers, 2 fn13
      writes in the GROW_SHRINK block), so the diagnostic field advances on
      the same inputs as the emulator.  Case subdirectory names are unchanged
      (20170714_1aer_org, 20170714_3aer_org, 20170714_1aer_mid,
      20170714_3aer_mid) and all four exist in the new tree.
      The four ARG cases in WRF_dm are UNCHANGED.

  23. All figures and CSVs now go to ./dcre_figures_more_fix (was
      ./east_boundary_traced_subdomains_lwp_v31), so the v20 output is
      preserved for comparison.  Override with WRF_OUTDIR as before.

  NOTE ON THE TIME WINDOW.  As in v18-v20 every case is restricted to the
  set intersection of timestamps across all eight cases.  The _more_fix ML
  runs were still integrating when this version was written and are much
  shorter than the ARG runs, so that intersection is currently bounded by
  the shortest ML case rather than by the ARG record.  Check the "common
  timestamps" count printed at startup before interpreting the dCRE bars:
  a short window means large sampling error on the decomposition terms.

Changes from v18:
  21. ML (non-ARG) cases now read the freshly re-run "*_new" wrfout directories
      (under /scratch/07088/tg863871/Perlm_Backup):
        org      = 20170714_org_new
        org_3aer = 20170714_org_3aer_new
        mid      = 20170714_mid_new
        mid_3aer = 20170714_mid_3aer_new
      These carry wrfout + wrfrst at every step (CCN1-6 read from wrfrst as
      before).  The four ARG cases are UNCHANGED (no "*_new" ARG runs exist).
      As in v18, all cases are restricted to the time steps shared by every
      ML and ARG case (set intersection of timestamps).
      Default output dir bumped to ./east_boundary_traced_subdomains_lwp_v26.

Changes from v17:
  18. ML (non-ARG) cases now read the same wrfout directories used by
      evaluate_activation_ML_vs_ARG_actcalled.py (the "*_out" emulator runs):
        org      = 20170714_org_09_out
        org_3aer = 20170714_org_3aer_out
        mid      = 20170714_mid_09_out
        mid_3aer = 20170714_mid_3aer_out

  19. All cases are restricted to the time steps SHARED by every ML and ARG
      case (set intersection of timestamps across all 8 cases).  Tracing,
      time series, dCRE decompositions and figures therefore sample exactly
      the same instants for every case.

  20. Bar-plot error bars now show the STANDARD ERROR of the mean
      (std / sqrt(n)) instead of the standard deviation.

Changes from v16:
  11. Four ARG cases added alongside the existing four ML cases:
        ARG_org      = 20170714_ARG           (ARG control, org background)
        ARG_org_3aer = 20170714_3aer_ARG      (ARG 3× aer,  org background)
        ARG_mid      = 20170714_ARG_mid       (ARG control, mid background)
        ARG_mid_3aer = 20170714_3aer_ARG_mid  (ARG 3× aer,  mid background)
      Two new DECOMP_PAIRS added (scheme="ARG").

  12. RAINNC accumulated precipitation read every timestep; instantaneous
      rate (mm hr⁻¹) computed by forward/backward finite-differencing and
      added as the 8th met-binned variable ("precip_rate").

  13. Vertical velocity stored and displayed in cm s⁻¹ (×100 from WRF m s⁻¹).

  14. X-axis tick labels in Figure 6 use ≤2 significant digits (:.2g).

  15. Figure 8 (NEW): domain-mean dCRE decomposition, ML vs ARG side by side
      for org and mid backgrounds, pixel-wise and mean-based columns.
      Saved as dCRE_ML_ARG_domain_comparison.png/.pdf

  16. Per-timestep Figures 1 & 2 continue to show only the 4 ML cases
      (PLOT_MAP_CASES) to keep map panels manageable.

  17. Figure 9 (NEW): dispersion (β–Nc) consistency check.  Recovers the
      gamma-PSD shape factor β = r_e / r_v (r_v = (3·LWP/4π ρ_w N_c)^{1/3})
      from RE_QC, LWP_TOT and column N_c, and regresses ln β on ln N_c per
      case (and perturbed−control per pair).  The slope b = dlnβ/dlnNc gives
      the dispersion offset to the Twomey number term, −3b·100%, compared
      against Wang et al. (2025, |b|≈0.024 ≈7%).  This is the same factor that
      DCRE_A_re isolates (Liu & Daum 2002).  Outputs:
        dispersion_beta_Nc_scatter.png/.pdf,
        dispersion_twomey_offset.png/.pdf,
        dispersion_beta_Nc_by_case.csv, dispersion_beta_Nc_by_pair.csv

All v16 features retained.
"""
import os
import sys
import glob
import warnings
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from datetime import datetime
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import matplotlib.dates as mdates

sys.stdout.reconfigure(line_buffering=True)

# Benign NumPy reductions over cloud-free / empty-in-cloud columns return NaN
# (handled downstream with has_cloud / nanmean masks).  Silence only these two
# specific messages so genuinely new RuntimeWarnings stay visible.
warnings.filterwarnings("ignore", message="All-NaN slice encountered")
warnings.filterwarnings("ignore", message="Mean of empty slice")


# ══════════════════════════════════════════════════════════════════════════════
# Case directories
# ══════════════════════════════════════════════════════════════════════════════
BASE_DIR = "/scratch/07088/tg863871/Perlm_Backup/WRF_stam3_ml_mixout_bce_lr1e-4_cosine_more_fix/test"

# v23: the ARG cases now come from WRF_dm_v2, not WRF_dm.
#
# WRF_dm_v2 re-runs the SAME four ARG-ACT cases from the same 09:02 restart; the
# trees differ only in Registry.EM_COMMON, where re_cloud/re_ice/re_snow are
# promoted to "rh" AND added to the morr_two_moment package.  In WRF_dm that
# package line omitted them, so re_cloud was never filled and TAU_QC_TOT / RE_QC
# came out identically zero -- which is why every earlier version of this
# analysis had to fall back to the analytic 1/3:2/3 Twomey split for the ARG
# sub-terms.  In WRF_dm_v2 both fields are populated, so R_tau and R_re are real
# measurements and the ARG apportionment is no longer a fallback.
#
# Override with WRF_ARG_DIR to point back at WRF_dm and reproduce the old run.
ARG_DIR = os.environ.get("WRF_ARG_DIR",
                         "/scratch/07088/tg863871/WRF_dm_v2/test")
case_dirs = {
    # ML (non-ARG) cases — v19: freshly re-run "*_new" emulator wrfout dirs
    "org":          os.path.join(BASE_DIR, "20170714_1aer_org"),
    "org_3aer":     os.path.join(BASE_DIR, "20170714_3aer_org"),
    "mid":          os.path.join(BASE_DIR, "20170714_1aer_mid"),
    "mid_3aer":     os.path.join(BASE_DIR, "20170714_3aer_mid"),
    # ARG cases
    "ARG_org":      os.path.join(ARG_DIR, "20170714_org_1aer_ARG"),
    "ARG_org_3aer": os.path.join(ARG_DIR, "20170714_org_3aer_ARG"),
    "ARG_mid":      os.path.join(ARG_DIR, "20170714_mid_1aer_ARG"),
    "ARG_mid_3aer": os.path.join(ARG_DIR, "20170714_mid_3aer_ARG"),
    # v24: the mid2 background -- mid with a moisture-only edit placing it
    # halfway between org and mid (see mid2_inputs_comparison.png).  Aerosol
    # is bit-identical to the mid pair, so the 1aer/3aer contrast is the same.
    "mid2":           os.path.join(BASE_DIR, "20170714_1aer_mid2"),
    "mid2_3aer":      os.path.join(BASE_DIR, "20170714_3aer_mid2"),
    "ARG_mid2":       os.path.join(ARG_DIR,  "20170714_mid2_1aer_ARG"),
    "ARG_mid2_3aer":  os.path.join(ARG_DIR,  "20170714_mid2_3aer_ARG"),
}

# Only ML cases shown in per-timestep LWP map panels (Figures 1 & 2)
ML_CASES       = ["org", "org_3aer", "mid", "mid_3aer", "mid2", "mid2_3aer"]
ARG_CASES      = ["ARG_org", "ARG_org_3aer", "ARG_mid", "ARG_mid_3aer",
                  "ARG_mid2", "ARG_mid2_3aer"]
# The LWP map grid is 2x2, so it still shows only the original four.
PLOT_MAP_CASES = ["org", "org_3aer", "mid", "mid_3aer"]

OUTDIR = "./dcre_figures_v24_mid2_1300"

# ── user options ──────────────────────────────────────────────────────────────
BOX_SIZE = 100
X_MARGIN = 2
Y_MARGIN = 10

TRACE_TOP_M_AGL       = 1000.0
WIND_VECTOR_TOP_M_AGL = 1000.0
MAX_SUBSTEP_SEC       = 60

LWP_SCALE   = 1000.0
Q_THRESHOLD = 1.0e-6

USE_FIXED_LWP_COLORBAR = True
FIXED_LWP_VMIN = 0.0
FIXED_LWP_VMAX = 200.0

WIND_SKIP    = 25
QUIVER_SCALE = 250
WIND_COLOR   = "white"

BOX_COLORS = ["tab:blue", "tab:orange", "tab:green", "tab:red"]
BOX_LABELS = ["Box 1", "Box 2", "Box 3", "Box 4"]

CASE_LINESTYLES = {
    "org":          "-",
    "org_3aer":     "--",
    "mid":          "-.",
    "mid_3aer":     ":",
    "ARG_org":      "-",
    "ARG_org_3aer": "--",
    "ARG_mid":      "-.",
    "ARG_mid_3aer": ":",
}

DECOMP_PAIRS = {
    "org_3aer-org":         {"perturbed": "org_3aer",     "control": "org",
                             "linestyle": "--", "scheme": "ML",  "background": "org"},
    "mid_3aer-mid":         {"perturbed": "mid_3aer",     "control": "mid",
                             "linestyle": ":",  "scheme": "ML",  "background": "mid"},
    "ARG_org_3aer-ARG_org": {"perturbed": "ARG_org_3aer", "control": "ARG_org",
                             "linestyle": "--", "scheme": "ARG", "background": "org"},
    "ARG_mid_3aer-ARG_mid": {"perturbed": "ARG_mid_3aer", "control": "ARG_mid",
                             "linestyle": ":",  "scheme": "ARG", "background": "mid"},
    "mid2_3aer-mid2":       {"perturbed": "mid2_3aer",    "control": "mid2",
                             "linestyle": "-.", "scheme": "ML",  "background": "mid2"},
    "ARG_mid2_3aer-ARG_mid2": {"perturbed": "ARG_mid2_3aer", "control": "ARG_mid2",
                             "linestyle": "-.", "scheme": "ARG", "background": "mid2"},
}

ML_PAIRS  = ["org_3aer-org", "mid_3aer-mid", "mid2_3aer-mid2"]
ARG_PAIRS = ["ARG_org_3aer-ARG_org", "ARG_mid_3aer-ARG_mid",
             "ARG_mid2_3aer-ARG_mid2"]

# Background grouping for Figure 8
BACKGROUND_GROUPS = {
    "org": {"ML": "org_3aer-org",         "ARG": "ARG_org_3aer-ARG_org",
            "label": "org background"},
    "mid": {"ML": "mid_3aer-mid",         "ARG": "ARG_mid_3aer-ARG_mid",
            "label": "mid background"},
    "mid2": {"ML": "mid2_3aer-mid2",      "ARG": "ARG_mid2_3aer-ARG_mid2",
             "label": "mid2 background"},
}

LAUNCHED_TRACE_CASE = "org"

# ── analysis window ───────────────────────────────────────────────────────────
# Hard end of the dCRE window.  Applied to the shared-time set in main(), so it
# caps EVERYTHING built downstream at once: the half-hourly stacked time series,
# the time-mean bars and the totals printed with them.  Set WRF_TEND="" to keep
# the whole shared record.
#
# 12:58 is the last output time mid_3aer -- the shortest of the eight records --
# had written when this window was fixed.  Pinning it here rather than letting
# the intersection float keeps the window identical from run to run while the
# longer cases keep integrating.
TIME_END_DEFAULT = "2017-07-15_13:00:00"


def _parse_window_time(s):
    """'2017-07-15_13:00:00', '2017-07-15 13:00' or '13:00' -> datetime/None."""
    s = (s or "").strip().replace("_", " ")
    if not s:
        return None
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%H:%M:%S", "%H:%M"):
        try:
            t = datetime.strptime(s, fmt)
        except ValueError:
            continue
        if fmt.startswith("%H"):                      # time only -> run date
            d = datetime.strptime(TIME_END_DEFAULT[:10], "%Y-%m-%d")
            t = t.replace(year=d.year, month=d.month, day=d.day)
        return t
    raise SystemExit(f"cannot parse WRF_TEND={s!r}")


TIME_END = _parse_window_time(os.environ.get("WRF_TEND", TIME_END_DEFAULT))

# ── environment overrides ─────────────────────────────────────────────────────
BOX_SIZE = int(os.environ.get("WRF_BOX_SIZE", str(BOX_SIZE)))
OUTDIR   = os.environ.get("WRF_OUTDIR", OUTDIR)
os.makedirs(OUTDIR, exist_ok=True)

# v22: the per-timestep map figures (1 & 2) are ~2.8 MB each, two per output
# time, and take the bulk of the runtime.  Off unless asked for.
MAKE_MAP_FIGURES = os.environ.get("WRF_MAP_FIGURES", "0") not in ("0", "", "no")

# v22 Figure 10: hourly ML-vs-ARG comparison.  An hour with only a couple of
# output times gives an SE that is noise, so require a floor.
HOURLY_BIN_MINUTES   = 60
MIN_TIMES_PER_HOUR   = 5

# v22 Figure 11: half-hourly stacked time series.
STACK_BIN_MINUTES    = 30
MIN_TIMES_PER_STACK  = 2

N_MET_BINS = 10

DOMAIN_LINE_COLOR = "#2d2d2d"
DOMAIN_LINE_WIDTH = 2.8

CF_MIN = 0.05
EPS    = 1e-10
GRAV   = 9.81
RD     = 287.04
RHO_W  = 1000.0    # liquid water density (kg m⁻³), for r_v in the β–Nc check

QCLOUD_THRESH_3D = 1.0e-5
FDIAG_TOP_M      = 1200.0
QNDROP_SCALE     = 1.0e-6

P0    = 1.0e5
RD_CP = RD / 1004.5

# colours
_METHOD_COLORS = {"pixel-wise": "#4393c3", "mean-based": "#d6604d"}

# pair markers / colours (ARG pairs get '^'/'D' markers, darker shades)
_PAIR_MARKERS = {
    "org_3aer-org":         "o",
    "mid_3aer-mid":         "s",
    "ARG_org_3aer-ARG_org": "^",
    "ARG_mid_3aer-ARG_mid": "D",
    "mid2_3aer-mid2":         "v",
    "ARG_mid2_3aer-ARG_mid2": "P",
}
_PAIR_BAR_COLORS = {
    "org_3aer-org":         "#4393c3",
    "mid_3aer-mid":         "#d6604d",
    "ARG_org_3aer-ARG_org": "#4393c3",   # same hue, distinguished by hatch
    "ARG_mid_3aer-ARG_mid": "#d6604d",
    "mid2_3aer-mid2":         "#4a3aa7",   # violet, validated against the other two
    "ARG_mid2_3aer-ARG_mid2": "#4a3aa7",
}
_PAIR_HATCH = {
    "org_3aer-org":         "",
    "mid_3aer-mid":         "",
    "ARG_org_3aer-ARG_org": "///",
    "ARG_mid_3aer-ARG_mid": "///",
    "mid2_3aer-mid2":         "",
    "ARG_mid2_3aer-ARG_mid2": "///",
}
_PAIR_ALPHA = {
    "org_3aer-org":         0.85,
    "mid_3aer-mid":         0.85,
    "ARG_org_3aer-ARG_org": 0.55,
    "ARG_mid_3aer-ARG_mid": 0.55,
    "mid2_3aer-mid2":         0.85,
    "ARG_mid2_3aer-ARG_mid2": 0.55,
}

_LWP_GROUP_COLORS = {"low": "#74add1", "high": "#f46d43"}

_DCRE_BAR_COLS   = ["DCRE_total", "DCRE_CF", "DCRE_A_Nc",
                    "DCRE_A_LWP", "DCRE_A_re", "DCRE_A_cov"]
_DCRE_BAR_LABELS = [r"$\Delta$CRE$_{\rm total}$",
                    r"$\Delta$CRE$_{\rm CF}$",
                    r"$\Delta$CRE$_{A,Nc}$",
                    r"$\Delta$CRE$_{A,LWP}$",
                    r"$\Delta$CRE$_{A,re}$",
                    r"$\Delta$CRE$_{A,cov}$"]

_DCRE_COMP_SPECS = [
    ("DCRE_total",  "black",    2.0, r"$\Delta$CRE$_{\rm total}$"),
    ("DCRE_CF",     "#1f77b4",  1.5, r"$\Delta$CRE$_{\rm CF}$"),
    ("DCRE_A_Nc",   "#ff7f0e",  1.3, r"$\Delta$CRE$_{A,Nc}$ (Twomey)"),
    ("DCRE_A_LWP",  "#2ca02c",  1.3, r"$\Delta$CRE$_{A,LWP}$ (2nd indirect)"),
    ("DCRE_A_re",   "#d62728",  1.1, r"$\Delta$CRE$_{A,re}$ (beyond-Twomey $r_e$)"),
    ("DCRE_A_cov",  "#9467bd",  0.8, r"$\Delta$CRE$_{A,cov}$ (RRTMG resid.)"),
]

# ── CCN variables read from wrfrst files (not in wrfout for ARG cases) ────────
CCN_VARS_RST = ("ccn1", "ccn2", "ccn3")   # internal keys; WRF names are CCN1/CCN2/CCN3
CCN_RST_MAP  = {"ccn1": "CCN1", "ccn2": "CCN2", "ccn3": "CCN3"}   # key → WRF varname

# ── met-bin variable definitions (11 variables; w in cm s⁻¹; CCN from wrfrst) ─
_MET_VAR_DEFS = [
    ("cldfrac",     "Cloud fraction",          ""),
    ("lwp",         "LWP",                     r"(g m$^{-2}$)"),
    ("qndrop",      "QNDROP",                  r"($\times 10^6$ kg$^{-1}$)"),
    ("w",           "Vertical velocity",        r"(cm s$^{-1}$)"),
    ("wspd",        "Horiz. wind speed",        r"(m s$^{-1}$)"),
    ("cld_thick",   "Cloud thickness",          "(m)"),
    ("rh",          "Relative humidity",        "(%)"),
    ("precip_rate", "Precip. rate (RAINNC)",    r"(mm hr$^{-1}$)"),
    ("ccn1",        "CCN1",                    r"(cm$^{-3}$)"),
    ("ccn2",        "CCN2",                    r"(cm$^{-3}$)"),
    ("ccn3",        "CCN3",                    r"(cm$^{-3}$)"),
]

_STACK_SPECS = [
    ("DCRE_CF",    "#1f77b4", r"$\Delta$CRE$_{\rm CF}$"),
    ("DCRE_A_Nc",  "#ff7f0e", r"$\Delta$CRE$_{A,Nc}$"),
    ("DCRE_A_LWP", "#2ca02c", r"$\Delta$CRE$_{A,LWP}$"),
    ("DCRE_A_re",  "#d62728", r"$\Delta$CRE$_{A,re}$"),
    ("DCRE_A_cov", "#9467bd", r"$\Delta$CRE$_{A,cov}$"),
]


# ══════════════════════════════════════════════════════════════════════════════
# Error-bar helpers  (v18: standard error of the mean instead of std)
# ══════════════════════════════════════════════════════════════════════════════

def sem_series(series):
    """Standard error of the mean of a pandas Series (NaNs skipped)."""
    n = series.count()
    return series.std() / np.sqrt(n) if n > 1 else np.nan


def sem_list(vals):
    """Standard error of the mean of a list/array of values (NaNs skipped)."""
    a = np.asarray(vals, dtype=float)
    a = a[np.isfinite(a)]
    n = a.size
    return np.std(a, ddof=1) / np.sqrt(n) if n > 1 else np.nan


# ══════════════════════════════════════════════════════════════════════════════
# WRF helpers  (unchanged from v16)
# ══════════════════════════════════════════════════════════════════════════════

def parse_time_from_ds(ds):
    t = ds["Times"].values[0]
    if isinstance(t, np.ndarray):
        s = b"".join(t).decode("utf-8")
    elif isinstance(t, bytes):
        s = t.decode("utf-8")
    else:
        s = str(t)
    return datetime.strptime(s, "%Y-%m-%d_%H:%M:%S")


def get_case_files(case_dir):
    files = sorted(glob.glob(os.path.join(case_dir, "wrfout*")))
    out = []
    for f in files:
        try:
            with xr.open_dataset(f, decode_times=False) as ds:
                t = parse_time_from_ds(ds)
            out.append((t, f))
        except Exception as e:
            print(f"  Skip {f}: {e}")
    return sorted(out, key=lambda x: x[0])


def get_case_rst_files(case_dir):
    """Return {datetime: filepath} for all wrfrst files in case_dir.

    wrfrst files carry CCN1-CCN6 which are absent from wrfout in ARG runs.
    The timestamp is parsed from the filename (wrfrst_d01_YYYY-MM-DD_HH:MM:SS).
    """
    import re as _re
    rst_by_time = {}
    for f in sorted(glob.glob(os.path.join(case_dir, "wrfrst*"))):
        m = _re.search(r"wrfrst_d\d+_(\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:\d{2})", f)
        if m:
            try:
                t = datetime.strptime(m.group(1), "%Y-%m-%d_%H:%M:%S")
                rst_by_time[t] = f
            except ValueError:
                pass
    return rst_by_time


def read_ccn_from_rst(rst_file, z_agl, ccn_keys=CCN_VARS_RST, top_m=1000.0):
    """Read CCN1/CCN2/CCN3 from a wrfrst file and return 0-1 km layer means.

    Uses netCDF4 directly (faster than xarray for a single variable read).
    z_agl : numpy array (nz, ny, nx) already computed from the matching wrfout.
    Returns dict {ccn_key: 2D array (ny, nx)}, values NaN where variable absent.
    """
    import netCDF4 as _nc4
    out = {}
    try:
        with _nc4.Dataset(rst_file, "r") as ds:
            for key in ccn_keys:
                wrf_name = CCN_RST_MAP[key]
                if wrf_name not in ds.variables:
                    out[key] = None
                    continue
                # shape (Time, nz, ny, nx) or (nz, ny, nx)
                arr = np.array(ds.variables[wrf_name][:], dtype=np.float32)
                while arr.ndim > 3 and arr.shape[0] == 1:
                    arr = arr[0]
                # 0-1 km AGL layer mean (all columns, not cloud-masked)
                out[key] = layer_mean_0_1km(arr, z_agl, top_m=top_m)
    except Exception as e:
        print(f"    [CCN wrfrst] {rst_file}: {e}")
        for key in ccn_keys:
            out.setdefault(key, None)
    return out


def destagger_x(u_stag):
    return 0.5 * (u_stag[:, :, :-1] + u_stag[:, :, 1:])


def destagger_y(v_stag):
    return 0.5 * (v_stag[:, :-1, :] + v_stag[:, 1:, :])


def get_height_agl_numpy(ds):
    ph  = ds["PH"].isel(Time=0).values
    phb = ds["PHB"].isel(Time=0).values
    z_w    = (ph + phb) / GRAV
    z_mass = 0.5 * (z_w[:-1, :, :] + z_w[1:, :, :])
    hgt    = ds["HGT"].isel(Time=0).values
    return z_mass - hgt[None, :, :]


def layer_mean_0_1km(field3d, z_agl, top_m=1000.0):
    mask         = (z_agl >= 0.0) & (z_agl <= top_m)
    field_masked = np.where(mask, field3d, np.nan)
    mean2d       = np.nanmean(field_masked, axis=0)
    missing      = ~np.isfinite(mean2d)
    if np.any(missing):
        mean2d[missing] = field3d[0, :, :][missing]
    return mean2d


def layer_mean_incloud(field3d, z_agl, qcloud_3d, top_m=1200.0, qc_thresh=1.0e-5):
    mask         = (z_agl >= 0.0) & (z_agl <= top_m) & (qcloud_3d > qc_thresh)
    field_masked = np.where(mask, field3d, np.nan)
    return np.nanmean(field_masked, axis=0)


def _extract_extra_met(ds0, z_agl, ny, nx):
    try:
        W_stag = ds0["W"].values
        w_mass = 0.5 * (W_stag[:-1] + W_stag[1:])
        # convert m/s → cm/s here
        w_2d   = layer_mean_0_1km(w_mass, z_agl) * 100.0
    except KeyError:
        w_2d = np.full((ny, nx), np.nan)

    try:
        u3d     = destagger_x(ds0["U"].values)
        v3d     = destagger_y(ds0["V"].values)
        u2d     = layer_mean_0_1km(u3d, z_agl)
        v2d     = layer_mean_0_1km(v3d, z_agl)
        wspd_2d = np.sqrt(u2d**2 + v2d**2)
    except KeyError:
        wspd_2d = np.full((ny, nx), np.nan)

    try:
        qc3d         = np.clip(ds0["QCLOUD"].values, 0, None)
        cloud_mask   = qc3d > QCLOUD_THRESH_3D
        has_cloud    = cloud_mask.any(axis=0)
        z_cld        = np.where(cloud_mask, z_agl, np.nan)
        cld_top      = np.nanmax(z_cld, axis=0)
        cld_base     = np.nanmin(z_cld, axis=0)
        cld_thick_2d = np.where(has_cloud, cld_top - cld_base, np.nan)
    except KeyError:
        cld_thick_2d = np.full((ny, nx), np.nan)

    try:
        P_arr = ds0["P"].values + ds0["PB"].values
        T_K   = (ds0["T"].values + 300.0) * (P_arr / 1e5)**RD_CP
        T_C   = T_K - 273.15
        e_s   = 611.2 * np.exp(17.67 * T_C / (T_C + 243.5))
        q     = np.clip(ds0["QVAPOR"].values, 0, None)
        e     = q * P_arr / (0.622 + q)
        rh_3d = np.clip(e / (e_s + 1e-10) * 100.0, 0.0, 100.0)
        rh_2d = layer_mean_0_1km(rh_3d, z_agl)
    except KeyError:
        rh_2d = np.full((ny, nx), np.nan)

    return w_2d, wspd_2d, cld_thick_2d, rh_2d


def bilinear_2d(field, y, x):
    ny, nx = field.shape
    x  = np.clip(x, 0, nx - 1.001)
    y  = np.clip(y, 0, ny - 1.001)
    x0 = np.floor(x).astype(int)
    y0 = np.floor(y).astype(int)
    x1 = np.clip(x0 + 1, 0, nx - 1)
    y1 = np.clip(y0 + 1, 0, ny - 1)
    wx = x - x0; wy = y - y0
    f00 = field[y0, x0]; f10 = field[y0, x1]
    f01 = field[y1, x0]; f11 = field[y1, x1]
    return (f00*(1-wx)*(1-wy) + f10*wx*(1-wy)
          + f01*(1-wx)*wy    + f11*wx*wy)


def interp_latlon_from_index(xlat, xlon, y, x):
    return bilinear_2d(xlat, y, x), bilinear_2d(xlon, y, x)


def squeeze_time(ds):
    if "Time" in ds.dims:
        return ds.isel(Time=0)
    return ds


def get_lwp_var(ds):
    if "LWP_TOT" in ds:
        return ds["LWP_TOT"]
    if "LWP" in ds:
        return ds["LWP"]
    raise KeyError("Neither LWP_TOT nor LWP found in dataset.")


def normalize_cldfrac(cf):
    vmax = float(cf.max(skipna=True))
    if vmax > 2.0:
        return (cf / 100.0).clip(0.0, 1.0)
    return cf.clip(0.0, 1.0)


# ══════════════════════════════════════════════════════════════════════════════
# Box masks  (unchanged)
# ══════════════════════════════════════════════════════════════════════════════

def polygon_mask_from_box_numpy(pos_box, ny, nx):
    y_poly = pos_box[:, 0]; x_poly = pos_box[:, 1]
    ymin = max(0,      int(np.floor(np.nanmin(y_poly))))
    ymax = min(ny - 1, int(np.ceil( np.nanmax(y_poly))))
    xmin = max(0,      int(np.floor(np.nanmin(x_poly))))
    xmax = min(nx - 1, int(np.ceil( np.nanmax(x_poly))))
    mask = np.zeros((ny, nx), dtype=bool)
    mask[ymin:ymax+1, xmin:xmax+1] = True
    return mask


# ══════════════════════════════════════════════════════════════════════════════
# dCRE field extraction  (unchanged from v16)
# ══════════════════════════════════════════════════════════════════════════════

def extract_fields_for_dcre(ds0):
    P_arr = ds0["P"].values + ds0["PB"].values
    T_arr = (ds0["T"].values + 300.0) * (P_arr / 1e5) ** RD_CP
    rho   = P_arr / (RD * T_arr * (1.0 + 0.608 * ds0["QVAPOR"].values))
    dz    = np.diff((ds0["PH"].values + ds0["PHB"].values) / GRAV, axis=0)

    LWP = np.sum(np.clip(ds0["QCLOUD"].values, 0, None) * rho * dz, axis=0)
    Nc  = np.sum(np.clip(ds0["QNDROP"].values, 0, None) * rho * dz, axis=0)

    cf  = ds0["CLDFRAC2D"].values
    tau = ds0["TAU_QC_TOT"].values

    re_raw = ds0["RE_QC"].values
    if re_raw.ndim == 3:
        qc3d = np.clip(ds0["QCLOUD"].values, 0, None)
        wsum = qc3d.sum(axis=0)
        re   = np.where(wsum > EPS,
                        (re_raw * qc3d).sum(axis=0) / wsum,
                        np.nan) * 1e6
    else:
        re = re_raw * 1e6

    swdnt      = ds0["SWDNT"].values
    cre_sw_toa = -(ds0["SWUPT"].values  - ds0["SWUPTC"].values)
    cre_sw_sfc = ((ds0["SWDNB"].values  - ds0["SWUPB"].values) -
                  (ds0["SWDNBC"].values - ds0["SWUPBC"].values))
    cre_lw_toa = -(ds0["LWUPT"].values  - ds0["LWUPTC"].values)
    cre_lw_sfc = ((ds0["LWDNB"].values  - ds0["LWUPB"].values) -
                  (ds0["LWDNBC"].values - ds0["LWUPBC"].values))

    alpha = np.where(
        (cf > CF_MIN) & (swdnt > 5.0),
        np.clip(-cre_sw_toa / (cf * swdnt + EPS), 0.0, 1.0),
        np.nan)

    return dict(LWP=LWP, Nc=Nc, cf=cf, tau=tau, re=re, alpha=alpha,
                swdnt=swdnt,
                cre_sw_toa=cre_sw_toa, cre_sw_sfc=cre_sw_sfc,
                cre_lw_toa=cre_lw_toa, cre_lw_sfc=cre_lw_sfc)


# ══════════════════════════════════════════════════════════════════════════════
# dCRE pixel-wise decomposition  (unchanged from v16)
# ══════════════════════════════════════════════════════════════════════════════

def decompose_dcre_one_box(d1, d3, box_mask):
    nan = float("nan")
    f1, f3 = d1["cf"], d3["cf"]
    S0     = d1["swdnt"]
    cld1   = f1 > CF_MIN
    cld3   = f3 > CF_MIN
    sunlit = S0 > 5.0
    any_cld = (cld1 | cld3) & sunlit & box_mask
    both    = (cld1 & cld3 & sunlit & box_mask &
               np.isfinite(d1["alpha"]) & np.isfinite(d3["alpha"]))

    _keys = ["DCRE_total", "DCRE_CF", "DCRE_A",
             "DCRE_A_Nc", "DCRE_A_LWP", "DCRE_A_re", "DCRE_A_cov",
             "DCRE_sfc", "DCRE_lw_toa", "DCRE_lw_sfc",
             "frac_Nc", "frac_LWP", "frac_re", "frac_cov",
             "R_Nc", "R_LWP", "R_tau", "R_re",
             "n_any_pix", "n_both_pix"]

    if any_cld.sum() == 0:
        out = {k: nan for k in _keys}
        out["n_any_pix"] = 0; out["n_both_pix"] = 0
        return out

    a1_ext = np.where(cld1 & np.isfinite(d1["alpha"]), d1["alpha"], 0.0)
    a3_ext = np.where(cld3 & np.isfinite(d3["alpha"]), d3["alpha"], 0.0)
    a_bar = 0.5 * (a1_ext + a3_ext)
    da    = a3_ext - a1_ext
    f_bar = 0.5 * (f1 + f3)
    df    = f3 - f1

    DCRE_total_2d = np.where(any_cld, d3["cre_sw_toa"] - d1["cre_sw_toa"], np.nan)
    DCRE_CF_2d    = np.where(any_cld, -S0 * a_bar * df,                     np.nan)
    DCRE_A_2d     = np.where(any_cld, -S0 * f_bar * da,                     np.nan)

    mean_ln_Rtau = mean_ln_RNc = mean_ln_RLWP = mean_ln_Rre = 0.0
    if both.sum() > 0:
        Nc1_ic  = d1["Nc"][both]  / (f1[both] + EPS)
        Nc3_ic  = d3["Nc"][both]  / (f3[both] + EPS)
        LWP1_ic = d1["LWP"][both] / (f1[both] + EPS)
        LWP3_ic = d3["LWP"][both] / (f3[both] + EPS)
        re1_ic  = d1["re"][both]
        re3_ic  = d3["re"][both]
        ok = ((Nc1_ic  > EPS) & (Nc3_ic  > EPS) &
              (LWP1_ic > EPS) & (LWP3_ic > EPS) &
              (d1["tau"][both] > EPS) & (d3["tau"][both] > EPS) &
              np.isfinite(re1_ic) & np.isfinite(re3_ic) &
              (re1_ic > 0.1) & (re3_ic > 0.1))
        if ok.sum() > 0:
            mean_ln_Rtau = float(np.mean(np.log(d3["tau"][both][ok] / d1["tau"][both][ok])))
            mean_ln_RNc  = float(np.mean(np.log(Nc3_ic[ok]  / Nc1_ic[ok])))
            mean_ln_RLWP = float(np.mean(np.log(LWP3_ic[ok] / LWP1_ic[ok])))
            mean_ln_Rre  = float(np.mean(np.log(re3_ic[ok]  / re1_ic[ok])))

    mean_dln_Nc        = (1./3.) * mean_ln_RNc
    mean_dln_LWP       = (2./3.) * mean_ln_RLWP
    mean_ln_Rre_Twomey = (1./3.) * (mean_ln_RLWP - mean_ln_RNc)
    mean_ln_Rre_adj    = mean_ln_Rre - mean_ln_Rre_Twomey

    if abs(mean_ln_Rtau) > 0.01:
        frac_Nc  =  mean_dln_Nc     / mean_ln_Rtau
        frac_LWP =  mean_dln_LWP    / mean_ln_Rtau
        frac_re  = -mean_ln_Rre_adj / mean_ln_Rtau
        frac_cov =  1.0 - frac_Nc - frac_LWP - frac_re
    else:
        frac_Nc = 1./3.; frac_LWP = 2./3.; frac_re = 0.; frac_cov = 0.

    DCRE_A_Nc_2d  = DCRE_A_2d * frac_Nc
    DCRE_A_LWP_2d = DCRE_A_2d * frac_LWP
    DCRE_A_re_2d  = DCRE_A_2d * frac_re
    DCRE_A_cov_2d = DCRE_A_2d * frac_cov

    def nm(arr):
        v = float(np.nanmean(arr))
        return v if np.isfinite(v) else nan

    return dict(
        DCRE_total  = nm(DCRE_total_2d),
        DCRE_CF     = nm(DCRE_CF_2d),
        DCRE_A      = nm(DCRE_A_2d),
        DCRE_A_Nc   = nm(DCRE_A_Nc_2d),
        DCRE_A_LWP  = nm(DCRE_A_LWP_2d),
        DCRE_A_re   = nm(DCRE_A_re_2d),
        DCRE_A_cov  = nm(DCRE_A_cov_2d),
        DCRE_sfc    = nm(np.where(any_cld, d3["cre_sw_sfc"] - d1["cre_sw_sfc"], np.nan)),
        DCRE_lw_toa = nm(np.where(any_cld, d3["cre_lw_toa"] - d1["cre_lw_toa"], np.nan)),
        DCRE_lw_sfc = nm(np.where(any_cld, d3["cre_lw_sfc"] - d1["cre_lw_sfc"], np.nan)),
        frac_Nc  = frac_Nc,  frac_LWP = frac_LWP,
        frac_re  = frac_re,  frac_cov = frac_cov,
        R_Nc  = float(np.exp(mean_ln_RNc)),
        R_LWP = float(np.exp(mean_ln_RLWP)),
        R_tau = float(np.exp(mean_ln_Rtau)),
        R_re  = float(np.exp(mean_ln_Rre)),
        n_any_pix  = int(any_cld.sum()),
        n_both_pix = int(both.sum()),
    )


# ══════════════════════════════════════════════════════════════════════════════
# dCRE mean-based decomposition  (unchanged from v16)
# ══════════════════════════════════════════════════════════════════════════════

def decompose_dcre_mean_based_one_box(d1, d3, box_mask):
    nan = float("nan")
    cld1   = (d1["cf"] > CF_MIN) & box_mask
    cld3   = (d3["cf"] > CF_MIN) & box_mask
    sunlit = (d1["swdnt"] > 5.0) & box_mask
    # note: `valid` (built below) intersects these with joint finiteness
    n_cld1 = int(cld1.sum()); n_cld3 = int(cld3.sum())

    _keys = ["DCRE_total", "DCRE_CF", "DCRE_A",
             "DCRE_A_Nc", "DCRE_A_LWP", "DCRE_A_re", "DCRE_A_cov",
             "DCRE_sfc", "DCRE_lw_toa", "DCRE_lw_sfc",
             "frac_Nc", "frac_LWP", "frac_re", "frac_cov",
             "R_Nc", "R_LWP", "R_tau", "R_re",
             "n_any_pix", "n_both_pix"]
    out_nan = {k: nan for k in _keys}
    out_nan["n_any_pix"] = n_cld1; out_nan["n_both_pix"] = n_cld3

    if not np.any(sunlit) or not np.any(box_mask):
        return out_nan

    # ── JOINT VALID MASK ──────────────────────────────────────────────────
    # Every domain/box average below must be taken over the SAME set of
    # points, and only over points that are non-NaN in BOTH runs.  Averaging
    # each field independently with nanmean lets CF_c and CRE_c come from
    # different pixels than CF_p and CRE_p, so the difference CRE_p - CRE_c
    # mixes two different samples and the decomposition no longer closes.
    valid = (box_mask & sunlit
             & np.isfinite(d1["cf"])         & np.isfinite(d3["cf"])
             & np.isfinite(d1["cre_sw_toa"]) & np.isfinite(d3["cre_sw_toa"])
             & np.isfinite(d1["swdnt"]))
    if not np.any(valid):
        return out_nan

    S0 = float(np.mean(d1["swdnt"][valid]))
    if S0 < 5.0:
        return out_nan

    CF_c  = float(np.mean(d1["cf"][valid]))
    CF_p  = float(np.mean(d3["cf"][valid]))
    CRE_c = float(np.mean(d1["cre_sw_toa"][valid]))
    CRE_p = float(np.mean(d3["cre_sw_toa"][valid]))
    DCRE_total = CRE_p - CRE_c

    if CF_c < CF_MIN or CF_p < CF_MIN:
        out = out_nan.copy(); out["DCRE_total"] = DCRE_total
        return out

    alpha_c = np.clip(-CRE_c / (S0 * CF_c + EPS), 0.0, 1.0)
    alpha_p = np.clip(-CRE_p / (S0 * CF_p + EPS), 0.0, 1.0)
    alpha_bar = 0.5 * (alpha_c + alpha_p)
    CF_bar    = 0.5 * (CF_c   + CF_p)
    DCRE_CF = -S0 * alpha_bar * (CF_p - CF_c)
    DCRE_A  = -S0 * CF_bar   * (alpha_p - alpha_c)

    def _in_cloud_scalars(d, cld_mask):
        # In-cloud means are per-run by construction (the cloudy sets differ),
        # but each is still restricted to points where that field is finite.
        cld_mask = cld_mask & valid
        if not np.any(cld_mask):
            return nan, nan, nan, nan
        def _m(k):
            v = d[k][cld_mask]; f = np.isfinite(v)
            return float(np.mean(v[f])) if np.any(f) else nan
        LWP_ic = _m("LWP")
        Nc_ic  = _m("Nc")
        tau_ic = _m("tau")
        lwp_v = d["LWP"][cld_mask]; re_v = d["re"][cld_mask]
        ok = np.isfinite(re_v) & (lwp_v > EPS)
        re_ic = (float(np.sum(re_v[ok]*lwp_v[ok])/np.sum(lwp_v[ok]))
                 if np.any(ok) else nan)
        return LWP_ic, Nc_ic, tau_ic, re_ic

    LWP_c, Nc_c, tau_c, re_c = _in_cloud_scalars(d1, cld1)
    LWP_p, Nc_p, tau_p, re_p = _in_cloud_scalars(d3, cld3)

    frac_Nc = frac_LWP = frac_re = frac_cov = nan
    R_Nc = R_LWP = R_tau = R_re = nan
    DCRE_A_Nc = DCRE_A_LWP = DCRE_A_re = DCRE_A_cov = nan

    all_valid = all(np.isfinite(v) and v > EPS
                    for v in [LWP_c, LWP_p, Nc_c, Nc_p, tau_c, tau_p, re_c, re_p])
    if all_valid:
        dln_tau       = np.log(tau_p / tau_c)
        dln_Nc        = (1./3.) * np.log(Nc_p  / Nc_c)
        dln_LWP       = (2./3.) * np.log(LWP_p / LWP_c)
        dln_re_Twomey = (1./3.) * (np.log(LWP_p/LWP_c) - np.log(Nc_p/Nc_c))
        dln_re_adj    = np.log(re_p/re_c) - dln_re_Twomey
        R_Nc  = float(Nc_p/Nc_c);  R_LWP = float(LWP_p/LWP_c)
        R_tau = float(tau_p/tau_c); R_re  = float(re_p/re_c)
        if abs(dln_tau) > 0.01:
            frac_Nc  =  dln_Nc      / dln_tau
            frac_LWP =  dln_LWP     / dln_tau
            frac_re  = -dln_re_adj  / dln_tau
            frac_cov =  1.0 - frac_Nc - frac_LWP - frac_re
        else:
            frac_Nc = 1./3.; frac_LWP = 2./3.; frac_re = 0.; frac_cov = 0.
        DCRE_A_Nc  = DCRE_A * frac_Nc
        DCRE_A_LWP = DCRE_A * frac_LWP
        DCRE_A_re  = DCRE_A * frac_re
        DCRE_A_cov = DCRE_A * frac_cov

    def _bdiff(k):
        # Same rule: difference of means over one joint, all-finite sample.
        m = valid & np.isfinite(d1[k]) & np.isfinite(d3[k])
        if not np.any(m):
            return nan
        return float(np.mean(d3[k][m]) - np.mean(d1[k][m]))

    return dict(
        DCRE_total=DCRE_total, DCRE_CF=DCRE_CF, DCRE_A=DCRE_A,
        DCRE_A_Nc=DCRE_A_Nc, DCRE_A_LWP=DCRE_A_LWP,
        DCRE_A_re=DCRE_A_re, DCRE_A_cov=DCRE_A_cov,
        DCRE_sfc=_bdiff("cre_sw_sfc"),
        DCRE_lw_toa=_bdiff("cre_lw_toa"),
        DCRE_lw_sfc=_bdiff("cre_lw_sfc"),
        frac_Nc=frac_Nc, frac_LWP=frac_LWP,
        frac_re=frac_re, frac_cov=frac_cov,
        R_Nc=R_Nc, R_LWP=R_LWP, R_tau=R_tau, R_re=R_re,
        n_any_pix=n_cld1, n_both_pix=n_cld3,
    )


# ══════════════════════════════════════════════════════════════════════════════
# Box tracing  (unchanged from v16)
# ══════════════════════════════════════════════════════════════════════════════

def make_initial_boxes(ny, nx):
    x1 = nx - 1 - X_MARGIN; x0 = x1 - BOX_SIZE
    available_y = ny - 2 * Y_MARGIN; spacing = available_y / 4
    boxes = []
    for i in range(4):
        yc = Y_MARGIN + spacing * (i + 0.5)
        y0 = int(round(yc - BOX_SIZE / 2))
        y0 = max(0, min(y0, ny - BOX_SIZE - 1))
        corners = np.array([[y0, x0], [y0, x1], [y0+BOX_SIZE, x1],
                             [y0+BOX_SIZE, x0]], dtype=float)
        boxes.append(corners)
    return np.array(boxes)


def advance_rigid_boxes(pos, u2d, v2d, dt, dx, dy, ny, nx):
    for ibox in range(pos.shape[0]):
        mask = polygon_mask_from_box_numpy(pos[ibox], ny, nx)
        if mask.any():
            u_mean = float(np.nanmean(u2d[mask]))
            v_mean = float(np.nanmean(v2d[mask]))
        else:
            yc = pos[ibox, :, 0]; xc = pos[ibox, :, 1]
            u_mean = float(np.mean(bilinear_2d(u2d, yc, xc)))
            v_mean = float(np.mean(bilinear_2d(v2d, yc, xc)))
        pos[ibox, :, 1] += u_mean * dt / dx
        pos[ibox, :, 0] += v_mean * dt / dy
    pos[..., 1] = np.clip(pos[..., 1], 0, nx - 1)
    pos[..., 0] = np.clip(pos[..., 0], 0, ny - 1)
    return pos


def trace_one_case(case_name, files):
    file_by_time = dict(files)
    times = [t for t, _ in files]
    with xr.open_dataset(files[0][1], decode_times=False) as ds:
        u3d = destagger_x(ds["U"].isel(Time=0).values)
        _, ny, nx = u3d.shape
        dx   = float(ds.attrs["DX"]); dy = float(ds.attrs["DY"])
        xlat = ds["XLAT"].isel(Time=0).values
        xlon = ds["XLONG"].isel(Time=0).values
    pos = make_initial_boxes(ny, nx)
    states_by_time = {times[0]: pos.copy()}
    for it in range(len(files) - 1):
        t0, f0 = files[it]; t1, _ = files[it+1]
        dt_total = (t1 - t0).total_seconds()
        nsub = max(1, int(np.ceil(dt_total / MAX_SUBSTEP_SEC)))
        dt   = dt_total / nsub
        with xr.open_dataset(f0, decode_times=False) as ds:
            u3d   = destagger_x(ds["U"].isel(Time=0).values)
            v3d   = destagger_y(ds["V"].isel(Time=0).values)
            z_agl = get_height_agl_numpy(ds)
            u2d   = layer_mean_0_1km(u3d, z_agl, top_m=TRACE_TOP_M_AGL)
            v2d   = layer_mean_0_1km(v3d, z_agl, top_m=TRACE_TOP_M_AGL)
        for _ in range(nsub):
            advance_rigid_boxes(pos, u2d, v2d, dt, dx, dy, ny, nx)
        states_by_time[t1] = pos.copy()
    grid = {"ny": ny, "nx": nx, "dx": dx, "dy": dy, "XLAT": xlat, "XLONG": xlon}
    print(f"  {case_name}: traced {len(times)} times")
    return states_by_time, file_by_time, grid


def make_initial_boxes_nonoverlap(ny, nx):
    x1 = nx - 1 - X_MARGIN; x0 = x1 - BOX_SIZE
    available_y = ny - 2 * Y_MARGIN
    gap = max(2, (available_y - 4 * BOX_SIZE) // 5)
    boxes = []
    for i in range(4):
        y0 = Y_MARGIN + gap*(i+1) + BOX_SIZE*i
        y0 = max(0, min(y0, ny - BOX_SIZE - 1))
        corners = np.array([[y0, x0], [y0, x1], [y0+BOX_SIZE, x1],
                             [y0+BOX_SIZE, x0]], dtype=float)
        boxes.append(corners)
    return np.array(boxes)


def trace_all_launched_boxes(files, ny, nx, dx, dy):
    n_times = len(files); active = {}; all_states = {}
    for it, (t, f) in enumerate(files):
        pos_init = make_initial_boxes_nonoverlap(ny, nx)
        active[t] = pos_init.copy(); all_states[t] = {t: pos_init.copy()}
        if it < n_times - 1:
            t1, _ = files[it+1]
            dt_total = (t1 - t).total_seconds()
            nsub = max(1, int(np.ceil(dt_total / MAX_SUBSTEP_SEC))); dt = dt_total / nsub
            with xr.open_dataset(f, decode_times=False) as ds:
                u3d   = destagger_x(ds["U"].isel(Time=0).values)
                v3d   = destagger_y(ds["V"].isel(Time=0).values)
                z_agl = get_height_agl_numpy(ds)
                u2d   = layer_mean_0_1km(u3d, z_agl, top_m=TRACE_TOP_M_AGL)
                v2d   = layer_mean_0_1km(v3d, z_agl, top_m=TRACE_TOP_M_AGL)
            new_active = {}
            for t_launch, pos in active.items():
                p = pos.copy()
                for _ in range(nsub):
                    advance_rigid_boxes(p, u2d, v2d, dt, dx, dy, ny, nx)
                new_active[t_launch] = p; all_states[t_launch][t1] = p.copy()
            active = new_active
    print(f"  Launched boxes: {len(all_states)} trajectories")
    return all_states


def compute_launched_traj_data(launched_states, all_files, file_by_case_time, ny, nx):
    traj_lwp = {tl: {} for tl in launched_states}
    org_files = all_files[LAUNCHED_TRACE_CASE]
    for t, f in org_files:
        with xr.open_dataset(f, decode_times=False) as ds:
            lwp_field = get_lwp_var(squeeze_time(ds)).values * LWP_SCALE
        for t_launch, states in launched_states.items():
            if t not in states:
                continue
            pos = states[t]; box_lwp = np.full(4, np.nan)
            for ibox in range(4):
                mask = polygon_mask_from_box_numpy(pos[ibox], ny, nx)
                if np.any(mask):
                    box_lwp[ibox] = np.nanmean(lwp_field[mask])
            traj_lwp[t_launch][t] = box_lwp

    traj_dcre = {pname: {tl: [] for tl in launched_states} for pname in DECOMP_PAIRS}
    all_traj_times = set()
    for states in launched_states.values():
        all_traj_times.update(states.keys())

    for pname, pair in DECOMP_PAIRS.items():
        pcase = pair["perturbed"]; ccase = pair["control"]
        p_files = file_by_case_time[pcase]; c_files = file_by_case_time[ccase]
        print(f"  dCRE for launched trajectories — {pname} ...")
        for t in sorted(all_traj_times):
            if t not in p_files or t not in c_files:
                continue
            try:
                with xr.open_dataset(p_files[t], decode_times=False) as dsp, \
                     xr.open_dataset(c_files[t], decode_times=False) as dsc:
                    d3 = extract_fields_for_dcre(squeeze_time(dsp))
                    d1 = extract_fields_for_dcre(squeeze_time(dsc))
                for t_launch, states in launched_states.items():
                    if t not in states:
                        continue
                    pos = states[t]
                    for ibox in range(4):
                        mask = polygon_mask_from_box_numpy(pos[ibox], ny, nx)
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore", RuntimeWarning)
                            row = decompose_dcre_one_box(d1, d3, mask)
                        row["t_launch"] = t_launch; row["t_future"] = t; row["box"] = ibox
                        traj_dcre[pname][t_launch].append(row)
            except Exception as e:
                print(f"    Warning: dCRE launched traj at {t}: {e}")
    return traj_lwp, traj_dcre


# ══════════════════════════════════════════════════════════════════════════════
# Box/domain time series  (v17: adds RAINNC → precip_rate; w in cm/s)
# ══════════════════════════════════════════════════════════════════════════════

def compute_box_and_domain_time_series(case_results):
    box_ts = {}; domain_ts = {}

    # "rainnc" stored as accumulated mm; "precip_rate" derived after loop
    all_vars = ("lwp", "cldfrac", "qndrop", "fdiag", "swcre",
                "w", "wspd", "cld_thick", "rh", "rainnc")

    for case_name in case_dirs:
        states, file_by_time, grid = case_results[case_name]
        ny, nx = grid["ny"], grid["nx"]
        times  = sorted(states.keys())
        n = len(times)

        box_arrs = {v: np.full((n, 4), np.nan) for v in all_vars}
        box_arrs["precip_rate"] = np.full((n, 4), np.nan)
        dom_arrs = {v: np.full(n, np.nan) for v in all_vars}
        dom_arrs["precip_rate"] = np.full(n, np.nan)

        for it, t in enumerate(times):
            f = file_by_time[t]
            with xr.open_dataset(f, decode_times=False) as ds:
                ds0 = squeeze_time(ds)
                lwp = get_lwp_var(ds0).values * LWP_SCALE
                cf  = normalize_cldfrac(ds0["CLDFRAC2D"]).values
                z_agl     = get_height_agl_numpy(ds)
                qcloud_3d = np.clip(ds0["QCLOUD"].values, 0.0, 1.0)
                col_has_cloud = np.any(qcloud_3d > QCLOUD_THRESH_3D, axis=0)

                try:
                    qndrop_3d = ds0["QNDROP"].values
                    qndrop_2d = (layer_mean_incloud(qndrop_3d, z_agl, qcloud_3d,
                                                    top_m=FDIAG_TOP_M,
                                                    qc_thresh=QCLOUD_THRESH_3D)
                                 * QNDROP_SCALE)
                except KeyError:
                    qndrop_2d = np.full((ny, nx), np.nan)

                try:
                    nu0   = np.clip(ds0["nu0"].values,   0.0, 1e12)
                    ac0   = np.clip(ds0["ac0"].values,   0.0, 1e12)
                    nu0cw = np.clip(ds0["nu0cw"].values, 0.0, 1e12)
                    ac0cw = np.clip(ds0["ac0cw"].values, 0.0, 1e12)
                    total = nu0 + ac0 + nu0cw + ac0cw
                    fdiag_3d = np.where(total > 0, (nu0cw + ac0cw) / total, np.nan)
                    fdiag_2d = layer_mean_incloud(fdiag_3d, z_agl, qcloud_3d,
                                                  top_m=FDIAG_TOP_M,
                                                  qc_thresh=QCLOUD_THRESH_3D)
                except KeyError:
                    fdiag_2d = np.full((ny, nx), np.nan)

                try:
                    swcre_full = ds0["SWUPTC"].values - ds0["SWUPT"].values
                    swcre_2d   = np.where(col_has_cloud, swcre_full, np.nan)
                except KeyError:
                    swcre_2d = np.full((ny, nx), np.nan)

                # w already in cm/s from _extract_extra_met
                w_2d, wspd_2d, cld_thick_2d, rh_2d = \
                    _extract_extra_met(ds0, z_agl, ny, nx)

                try:
                    rainnc_2d = ds0["RAINNC"].values
                except (KeyError, AttributeError):
                    rainnc_2d = np.full((ny, nx), np.nan)

            fields = {
                "lwp":       lwp,
                "cldfrac":   cf,
                "qndrop":    qndrop_2d,
                "fdiag":     fdiag_2d,
                "swcre":     swcre_2d,
                "w":         w_2d,
                "wspd":      wspd_2d,
                "cld_thick": cld_thick_2d,
                "rh":        rh_2d,
                "rainnc":    rainnc_2d,
            }

            for v, fld in fields.items():
                dom_arrs[v][it] = np.nanmean(fld)

            pos = states[t]
            for ibox in range(4):
                mask = polygon_mask_from_box_numpy(pos[ibox], ny, nx)
                if np.any(mask):
                    for v, fld in fields.items():
                        box_arrs[v][it, ibox] = np.nanmean(fld[mask])

        # ── compute precip_rate from RAINNC differences ───────────────────
        rainnc_b = box_arrs["rainnc"]   # (n, 4)
        rainnc_d = dom_arrs["rainnc"]   # (n,)
        for it in range(n):
            if it == 0:
                if n >= 2:
                    dt_hr = (times[1] - times[0]).total_seconds() / 3600.0
                    if dt_hr > 0:
                        box_arrs["precip_rate"][0] = (
                            np.maximum(rainnc_b[1] - rainnc_b[0], 0.0) / dt_hr)
                        dom_arrs["precip_rate"][0] = (
                            max(float(rainnc_d[1] - rainnc_d[0]), 0.0) / dt_hr)
            else:
                dt_hr = (times[it] - times[it-1]).total_seconds() / 3600.0
                if dt_hr > 0:
                    box_arrs["precip_rate"][it] = (
                        np.maximum(rainnc_b[it] - rainnc_b[it-1], 0.0) / dt_hr)
                    dom_arrs["precip_rate"][it] = (
                        max(float(rainnc_d[it] - rainnc_d[it-1]), 0.0) / dt_hr)

        times_arr = np.array(times)
        box_ts[case_name]    = {"times": times_arr, **box_arrs}
        domain_ts[case_name] = {"times": times_arr, **dom_arrs}
        print(f"  {case_name}: box + domain TS done ({n} steps)")

    return box_ts, domain_ts


# ══════════════════════════════════════════════════════════════════════════════
# dCRE decomposition time series  (unchanged from v16, covers all 4 pairs)
# ══════════════════════════════════════════════════════════════════════════════

def compute_dcre_decomp_time_series(case_results, file_by_case_time, ny, nx):
    box_decomp_pw = {}; domain_decomp_pw = {}
    box_decomp_mb = {}; domain_decomp_mb = {}
    domain_mask = np.ones((ny, nx), dtype=bool)

    for pair_name, pair in DECOMP_PAIRS.items():
        pcase = pair["perturbed"]; ccase = pair["control"]
        p_states, _, _ = case_results[pcase]
        p_files = file_by_case_time[pcase]; c_files = file_by_case_time[ccase]
        pair_times = sorted(set(p_files.keys()) & set(c_files.keys()) & set(p_states.keys()))
        print(f"  dCRE (pw + mb) for {pair_name}: {len(pair_times)} times")

        box_rows_pw = []; box_rows_mb = []
        dom_rows_pw = []; dom_rows_mb = []

        for t in pair_times:
            fp = p_files[t]; fc = c_files[t]
            try:
                with xr.open_dataset(fp, decode_times=False) as dsp, \
                     xr.open_dataset(fc, decode_times=False) as dsc:
                    d3 = extract_fields_for_dcre(squeeze_time(dsp))
                    d1 = extract_fields_for_dcre(squeeze_time(dsc))
                    pos = p_states[t]
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", RuntimeWarning)
                        for ibox in range(4):
                            mask = polygon_mask_from_box_numpy(pos[ibox], ny, nx)
                            row_pw = decompose_dcre_one_box(d1, d3, mask)
                            row_pw.update({"time": t, "box": ibox,
                                           "box_label": BOX_LABELS[ibox], "pair": pair_name})
                            box_rows_pw.append(row_pw)
                            row_mb = decompose_dcre_mean_based_one_box(d1, d3, mask)
                            row_mb.update({"time": t, "box": ibox,
                                           "box_label": BOX_LABELS[ibox], "pair": pair_name})
                            box_rows_mb.append(row_mb)
                        dom_pw = decompose_dcre_one_box(d1, d3, domain_mask)
                        dom_pw.update({"time": t, "pair": pair_name}); dom_rows_pw.append(dom_pw)
                        dom_mb = decompose_dcre_mean_based_one_box(d1, d3, domain_mask)
                        dom_mb.update({"time": t, "pair": pair_name}); dom_rows_mb.append(dom_mb)
            except Exception as e:
                print(f"    Warning: dCRE decomp failed for {pair_name} at {t}: {e}")

        def _to_df(rows):
            return (pd.DataFrame(rows).set_index("time").sort_index()
                    if rows else pd.DataFrame())

        box_decomp_pw[pair_name]    = _to_df(box_rows_pw)
        domain_decomp_pw[pair_name] = _to_df(dom_rows_pw)
        box_decomp_mb[pair_name]    = _to_df(box_rows_mb)
        domain_decomp_mb[pair_name] = _to_df(dom_rows_mb)

    return box_decomp_pw, domain_decomp_pw, box_decomp_mb, domain_decomp_mb


def compute_global_lwp_range(file_by_case_time):
    if USE_FIXED_LWP_COLORBAR:
        print(f"  Fixed LWP colorbar: {FIXED_LWP_VMIN}–{FIXED_LWP_VMAX} g m⁻²")
        return FIXED_LWP_VMIN, FIXED_LWP_VMAX
    global_min = np.inf; global_max = -np.inf
    for case_name in ML_CASES:   # range from ML cases only
        for _, f in sorted(file_by_case_time[case_name].items()):
            with xr.open_dataset(f, decode_times=False) as ds:
                lwp = get_lwp_var(squeeze_time(ds)).values * LWP_SCALE
            global_min = min(global_min, np.nanmin(lwp))
            global_max = max(global_max, np.nanmax(lwp))
    return (global_min, global_max) if global_min < global_max else (0, 1)


# ══════════════════════════════════════════════════════════════════════════════
# Plotting helpers
# ══════════════════════════════════════════════════════════════════════════════

def add_wind_vectors(ax, ds, grid):
    xlat = grid["XLAT"]; xlon = grid["XLONG"]
    u3d  = destagger_x(ds["U"].isel(Time=0).values)
    v3d  = destagger_y(ds["V"].isel(Time=0).values)
    z_agl = get_height_agl_numpy(ds)
    u2d  = layer_mean_0_1km(u3d, z_agl, top_m=WIND_VECTOR_TOP_M_AGL)
    v2d  = layer_mean_0_1km(v3d, z_agl, top_m=WIND_VECTOR_TOP_M_AGL)
    sl   = (slice(None, None, WIND_SKIP), slice(None, None, WIND_SKIP))
    ax.quiver(xlon[sl], xlat[sl], u2d[sl], v2d[sl],
              scale=QUIVER_SCALE, width=0.002, alpha=0.9, color=WIND_COLOR)


def plot_timeseries_panel(ax, box_ts, variable, ylabel, current_time,
                          domain_ts=None, ts_cases=None):
    """Plot box-mean time series for a given variable.

    ts_cases: list of case names to plot (default: ML_CASES only).
    """
    if ts_cases is None:
        ts_cases = ML_CASES
    for case_name in ts_cases:
        if case_name not in box_ts:
            continue
        times = box_ts[case_name]["times"]
        vals  = box_ts[case_name].get(variable)
        if vals is None:
            continue
        for ibox in range(4):
            ax.plot(times, vals[:, ibox],
                    color=BOX_COLORS[ibox],
                    linestyle=CASE_LINESTYLES.get(case_name, "-"),
                    linewidth=1.4, alpha=0.85)
    if domain_ts is not None:
        for case_name in ts_cases:
            if case_name not in domain_ts:
                continue
            times = domain_ts[case_name]["times"]
            vals  = domain_ts[case_name].get(variable)
            if vals is None:
                continue
            ax.plot(times, vals,
                    color=DOMAIN_LINE_COLOR,
                    linestyle=CASE_LINESTYLES.get(case_name, "-"),
                    linewidth=DOMAIN_LINE_WIDTH, alpha=0.75, zorder=5)
    ax.axvline(current_time, color="k", linestyle="--", linewidth=1.3)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)


def plot_lwp_map_panels(fig, map_axes, t, case_results, file_by_case_time,
                        global_vmin, global_vmax, plot_cases=None):
    """Plot LWP maps for the given cases (default: PLOT_MAP_CASES)."""
    if plot_cases is None:
        plot_cases = PLOT_MAP_CASES
    mappable = None
    for ax, case_name in zip(map_axes, plot_cases):
        if t not in file_by_case_time.get(case_name, {}):
            ax.set_title(f"{case_name}\nNo data"); ax.axis("off"); continue
        states, _, grid = case_results[case_name]
        xlat = grid["XLAT"]; xlon = grid["XLONG"]
        f = file_by_case_time[case_name][t]
        with xr.open_dataset(f, decode_times=False) as ds:
            ds0 = squeeze_time(ds)
            lwp = get_lwp_var(ds0).values * LWP_SCALE
            mappable = ax.pcolormesh(xlon, xlat, lwp, shading="auto",
                                     cmap="viridis",
                                     vmin=global_vmin, vmax=global_vmax)
            add_wind_vectors(ax, ds, grid)
        pos = states[t]
        for ibox in range(4):
            y_poly = pos[ibox, :, 0]; x_poly = pos[ibox, :, 1]
            lat_poly, lon_poly = interp_latlon_from_index(xlat, xlon, y_poly, x_poly)
            ax.plot(np.r_[lon_poly, lon_poly[0]], np.r_[lat_poly, lat_poly[0]],
                    color=BOX_COLORS[ibox], linewidth=2.6)
        ax.set_title(case_name); ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
    return mappable


def plot_dcre_decomp_box_panel(ax, decomp_ts_pw, current_time, ibox,
                               decomp_ts_mb=None, pairs=None):
    """Time-series dCRE for one box.  pairs limits which pair names are shown."""
    if pairs is None:
        pairs = list(decomp_ts_pw.keys())
    for pair_name, df in decomp_ts_pw.items():
        if pair_name not in pairs or df.empty:
            continue
        sub = df[df["box"] == ibox].sort_index()
        if sub.empty:
            continue
        marker = _PAIR_MARKERS.get(pair_name, "o")
        for comp, color, lw, _ in _DCRE_COMP_SPECS:
            if comp not in sub.columns:
                continue
            ax.plot(sub.index, sub[comp], linestyle="-",
                    marker=marker, ms=2.5, color=color, linewidth=lw, alpha=0.90)
    if decomp_ts_mb is not None:
        for pair_name, df in decomp_ts_mb.items():
            if pair_name not in pairs or df.empty:
                continue
            sub = df[df["box"] == ibox].sort_index()
            if sub.empty:
                continue
            marker = _PAIR_MARKERS.get(pair_name, "o")
            for comp, color, lw, _ in _DCRE_COMP_SPECS:
                if comp not in sub.columns:
                    continue
                ax.plot(sub.index, sub[comp], linestyle="--",
                        marker=marker, ms=2.0, color=color,
                        linewidth=lw * 0.75, alpha=0.65)
    ax.axhline(0, lw=0.5, color="k")
    ax.axvline(current_time, color="k", linestyle="--", linewidth=1.2)
    ax.set_title(BOX_LABELS[ibox], color=BOX_COLORS[ibox], fontsize=10)
    ax.set_ylabel(r"W m$^{-2}$")
    ax.grid(True, alpha=0.3)


def plot_domain_dcre_timeseries_panel(ax, domain_decomp_ts_pw, current_time,
                                      pair_name, domain_decomp_ts_mb=None):
    df = domain_decomp_ts_pw.get(pair_name, pd.DataFrame())
    if not df.empty:
        for comp, color, lw, _ in _DCRE_COMP_SPECS:
            if comp in df.columns:
                ax.plot(df.sort_index().index, df.sort_index()[comp],
                        linestyle="-", color=color, linewidth=lw, alpha=0.9)
    if domain_decomp_ts_mb is not None:
        df_mb = domain_decomp_ts_mb.get(pair_name, pd.DataFrame())
        if not df_mb.empty:
            for comp, color, lw, _ in _DCRE_COMP_SPECS:
                if comp in df_mb.columns:
                    ax.plot(df_mb.sort_index().index, df_mb.sort_index()[comp],
                            linestyle="--", color=color,
                            linewidth=lw * 0.75, alpha=0.65)
    ax.axhline(0, lw=0.5, color="k")
    ax.axvline(current_time, color="k", linestyle="--", linewidth=1.2)
    ax.set_title(f"Domain dCRE — {pair_name}", fontsize=9)
    ax.set_ylabel(r"W m$^{-2}$"); ax.grid(True, alpha=0.3)


# ── legends ───────────────────────────────────────────────────────────────────

def add_figure_legends(fig, show_domain_line=False, show_method_legend=False,
                       ts_cases=None):
    if ts_cases is None:
        ts_cases = ML_CASES
    color_handles = [Line2D([0], [0], color=BOX_COLORS[i], lw=3, linestyle="-",
                            label=BOX_LABELS[i]) for i in range(4)]
    case_handles  = [Line2D([0], [0], color="k", lw=2.2,
                            linestyle=CASE_LINESTYLES.get(c, "-"), label=c)
                     for c in ts_cases]
    if show_domain_line:
        case_handles.append(
            Line2D([0], [0], color=DOMAIN_LINE_COLOR, lw=DOMAIN_LINE_WIDTH,
                   linestyle="-", label="domain mean (linestyle=case)"))
    fig.legend(handles=color_handles, loc="lower center",
               bbox_to_anchor=(0.27, -0.03), ncol=4, frameon=True,
               title="Color = subdomain", fontsize=9)
    fig.legend(handles=case_handles, loc="lower center",
               bbox_to_anchor=(0.73, -0.03),
               ncol=min(4, len(case_handles)), frameon=True,
               title="Linestyle = case", fontsize=9)
    if show_method_legend:
        mh = [Line2D([0], [0], color="k", lw=1.8, linestyle="-",
                     label="pixel-wise (original)"),
              Line2D([0], [0], color="k", lw=1.4, linestyle="--",
                     label="mean-based (new)")]
        fig.legend(handles=mh, loc="lower center",
                   bbox_to_anchor=(0.5, -0.07), ncol=2, frameon=True,
                   title="Decomposition method", fontsize=9)


def add_dcre_decomp_legends(fig, decomp_ts, pairs=None):
    if pairs is None:
        pairs = list(decomp_ts.keys())
    comp_handles = [Line2D([0], [0], color=s[1], lw=s[2], linestyle="-", label=s[3])
                    for s in _DCRE_COMP_SPECS]
    pair_handles = [Line2D([0], [0], color="k", lw=1.5,
                           marker=_PAIR_MARKERS.get(p, "o"),
                           linestyle="None", label=p)
                    for p in pairs if p in decomp_ts and not decomp_ts[p].empty]
    fig.legend(handles=comp_handles, loc="lower center",
               bbox_to_anchor=(0.36, -0.07), ncol=3, fontsize=7.5, frameon=True,
               title="Color = dCRE component")
    if pair_handles:
        fig.legend(handles=pair_handles, loc="lower center",
                   bbox_to_anchor=(0.80, -0.07),
                   ncol=max(1, len(pair_handles)), fontsize=8, frameon=True,
                   title="Marker = aerosol pair")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 1: LWP maps + time series + dCRE panels  (ML cases / ML pairs only)
# ══════════════════════════════════════════════════════════════════════════════

def plot_one_time(t, case_results, file_by_case_time, box_ts, domain_ts,
                  decomp_ts_pw, domain_decomp_ts_pw,
                  decomp_ts_mb, domain_decomp_ts_mb,
                  global_vmin, global_vmax):
    fig = plt.figure(figsize=(16, 26))
    gs  = GridSpec(7, 2,
                   height_ratios=[1, 1, 0.42, 0.42, 0.58, 0.58, 0.58],
                   hspace=0.46, wspace=0.18, figure=fig)

    map_axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]),
                fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
    ax_lwp_ts    = fig.add_subplot(gs[2, :])
    ax_cf_ts     = fig.add_subplot(gs[3, :])
    decomp_axes  = [fig.add_subplot(gs[4, 0]), fig.add_subplot(gs[4, 1]),
                    fig.add_subplot(gs[5, 0]), fig.add_subplot(gs[5, 1])]
    # v24: three ML pairs now, so row 6 is sub-divided by len(ML_PAIRS)
    # instead of being split in two.
    _gs_dom = gs[6, :].subgridspec(1, len(ML_PAIRS), wspace=0.22)
    domain_dcre_axes = [fig.add_subplot(_gs_dom[0, i])
                        for i in range(len(ML_PAIRS))]

    mappable = plot_lwp_map_panels(fig, map_axes, t, case_results,
                                   file_by_case_time, global_vmin, global_vmax,
                                   plot_cases=PLOT_MAP_CASES)
    if mappable is None:
        plt.close(fig); return

    cbar = fig.colorbar(mappable, ax=map_axes, shrink=0.88, pad=0.02)
    cbar.set_label("LWP_TOT (g m$^{-2}$)")

    plot_timeseries_panel(ax_lwp_ts, box_ts, "lwp",
                          "Box-mean LWP\n(g m$^{-2}$)", t,
                          domain_ts=domain_ts, ts_cases=ML_CASES)
    plot_timeseries_panel(ax_cf_ts, box_ts, "cldfrac",
                          "Box-mean CF", t,
                          domain_ts=domain_ts, ts_cases=ML_CASES)

    # dCRE panels — ML pairs only to avoid overcrowding
    for ibox, ax_d in enumerate(decomp_axes):
        plot_dcre_decomp_box_panel(ax_d, decomp_ts_pw, current_time=t,
                                   ibox=ibox, decomp_ts_mb=decomp_ts_mb,
                                   pairs=ML_PAIRS)

    for pair_name, ax_dom in zip(ML_PAIRS, domain_dcre_axes):
        plot_domain_dcre_timeseries_panel(
            ax_dom, domain_decomp_ts_pw, t, pair_name,
            domain_decomp_ts_mb=domain_decomp_ts_mb)
        ax_dom.set_xlabel("Time")

    decomp_axes[2].set_xlabel("Time"); decomp_axes[3].set_xlabel("Time")
    add_figure_legends(fig, show_domain_line=True, show_method_legend=True,
                       ts_cases=ML_CASES)
    add_dcre_decomp_legends(fig, decomp_ts_pw, pairs=ML_PAIRS)

    fig.suptitle(
        "LWP_TOT, 0–1 km wind, moving subdomains and dCRE decomposition (ML pairs)\n"
        r"Solid = pixel-wise  |  Dashed = mean-based  |  Thick dark = domain mean"
        f"\nCurrent time: {t:%Y-%m-%d %H:%M:%S}",
        fontsize=12, y=0.993)

    outname = os.path.join(OUTDIR, f"LWP_wind_boxes_dcre_decomp_{t:%Y%m%d_%H%M%S}.png")
    fig.savefig(outname, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outname}")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 2: LWP maps + 5 full-width met time series  (ML cases only)
# ══════════════════════════════════════════════════════════════════════════════

def plot_additional_figure(t, case_results, file_by_case_time, box_ts,
                           domain_ts, global_vmin, global_vmax):
    fig = plt.figure(figsize=(16, 26))
    gs  = GridSpec(7, 2,
                   height_ratios=[1, 1, 0.45, 0.45, 0.45, 0.45, 0.45],
                   hspace=0.48, wspace=0.18, figure=fig)

    map_axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]),
                fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
    ax_lwp_ts    = fig.add_subplot(gs[2, :])
    ax_cf_ts     = fig.add_subplot(gs[3, :])
    ax_qndrop_ts = fig.add_subplot(gs[4, :])
    ax_fdiag_ts  = fig.add_subplot(gs[5, :])
    ax_swcre_ts  = fig.add_subplot(gs[6, :])

    mappable = plot_lwp_map_panels(fig, map_axes, t, case_results,
                                   file_by_case_time, global_vmin, global_vmax,
                                   plot_cases=PLOT_MAP_CASES)
    if mappable is None:
        plt.close(fig); return

    cbar = fig.colorbar(mappable, ax=map_axes, shrink=0.88, pad=0.02)
    cbar.set_label("LWP_TOT (g m$^{-2}$)")

    kw = dict(domain_ts=domain_ts, ts_cases=ML_CASES)
    plot_timeseries_panel(ax_lwp_ts,    box_ts, "lwp",
                          "LWP (g m$^{-2}$)", t, **kw)
    plot_timeseries_panel(ax_cf_ts,     box_ts, "cldfrac", "CF", t, **kw)
    plot_timeseries_panel(ax_qndrop_ts, box_ts, "qndrop",
                          "QNDROP ≤1.2 km in-cloud\n"
                          r"($\times$10$^6$ kg$^{-1}$)", t, **kw)
    plot_timeseries_panel(ax_fdiag_ts,  box_ts, "fdiag",
                          r"$f_{\rm diag}$ ≤1.2 km in-cloud", t, **kw)
    ax_fdiag_ts.set_ylim(0.0, 1.0)
    plot_timeseries_panel(ax_swcre_ts,  box_ts, "swcre",
                          "SWCRE cloudy cols\n(W m$^{-2}$)", t, **kw)
    ax_swcre_ts.axhline(0, color="k", lw=0.6, ls=":")
    ax_swcre_ts.set_xlabel("Time")

    add_figure_legends(fig, show_domain_line=True, ts_cases=ML_CASES)
    fig.suptitle(
        r"LWP, QNDROP & $f_{\rm diag}$ (≤1.2 km AGL, in-cloud) and SWCRE"
        " — ML cases, moving subdomains\n"
        r"Thick dark = domain mean  |  "
        f"Current time: {t:%Y-%m-%d %H:%M:%S}",
        fontsize=14, y=0.993)

    outname = os.path.join(OUTDIR, f"QNDROP_fdiag_SWCRE_boxes_{t:%Y%m%d_%H%M%S}.png")
    fig.savefig(outname, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outname}")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 3: Time-averaged dCRE bar plot — all 4 pairs
# ══════════════════════════════════════════════════════════════════════════════

def plot_dCRE_timeavg_barplot(decomp_ts):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharey=False)
    axes = axes.flatten()
    all_pairs = list(decomp_ts.keys())
    n_pairs   = len(all_pairs)
    bar_width = 0.20
    x = np.arange(len(_DCRE_BAR_COLS))

    for ibox, ax in enumerate(axes):
        for ip, pname in enumerate(all_pairs):
            df = decomp_ts.get(pname, pd.DataFrame())
            if df.empty:
                continue
            sub = df[df["box"] == ibox]
            if sub.empty:
                continue
            means  = [sub[c].mean() if c in sub else np.nan for c in _DCRE_BAR_COLS]
            stds   = [sem_series(sub[c]) if c in sub else np.nan for c in _DCRE_BAR_COLS]
            offset = (ip - n_pairs / 2.0 + 0.5) * bar_width
            ax.bar(x + offset, means, bar_width,
                   yerr=stds, capsize=3, error_kw={"linewidth": 1.0},
                   color=_PAIR_BAR_COLORS.get(pname, f"C{ip}"),
                   hatch=_PAIR_HATCH.get(pname, ""),
                   alpha=_PAIR_ALPHA.get(pname, 0.85),
                   label=pname)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_title(BOX_LABELS[ibox], color=BOX_COLORS[ibox],
                     fontsize=12, fontweight="bold")
        ax.set_ylabel(r"W m$^{-2}$")
        ax.set_xticks(x)
        ax.set_xticklabels(_DCRE_BAR_LABELS, rotation=30, ha="right", fontsize=9)
        ax.grid(True, alpha=0.3, axis="y")
        if ibox == 0:
            ax.legend(fontsize=7, frameon=True, loc="best", ncol=2)

    fig.suptitle(
        "Time-averaged dCRE decomposition by box — pixel-wise  |  all 4 pairs\n"
        "Solid fill = ML, hatched = ARG  |  Error bars = ±1 SE",
        fontsize=13)
    fig.tight_layout()
    outname = os.path.join(OUTDIR, "dCRE_timeavg_barplot_all4pairs.png")
    fig.savefig(outname, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outname}")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 4: LWP-stratified dCRE — ML pairs only (launched traj from ML winds)
# ══════════════════════════════════════════════════════════════════════════════

def plot_lwp_stratified_dcre_barplot(traj_lwp, traj_dcre):
    # filter to ML pairs only
    ml_traj_dcre = {p: v for p, v in traj_dcre.items() if p in ML_PAIRS}

    all_launch_times = sorted(traj_lwp.keys())
    mean_lwp_by_traj = {ibox: {} for ibox in range(4)}
    for t_launch in all_launch_times:
        vals_by_t = traj_lwp[t_launch]
        if not vals_by_t:
            continue
        stacked = np.array(list(vals_by_t.values()))
        for ibox in range(4):
            col = stacked[:, ibox]
            finite = col[np.isfinite(col)]
            if len(finite):
                mean_lwp_by_traj[ibox][t_launch] = float(np.mean(finite))

    medians = {}
    for ibox in range(4):
        vals = list(mean_lwp_by_traj[ibox].values())
        medians[ibox] = float(np.median(vals)) if vals else np.nan

    def get_group(ibox, t_launch):
        v = mean_lwp_by_traj[ibox].get(t_launch, np.nan)
        if not np.isfinite(v) or not np.isfinite(medians[ibox]):
            return None
        return "high" if v >= medians[ibox] else "low"

    groups = ["low", "high"]
    data = {pname: {ibox: {g: {c: [] for c in _DCRE_BAR_COLS}
                           for g in groups}
                    for ibox in range(4)}
            for pname in ml_traj_dcre}

    for pname, tl_dict in ml_traj_dcre.items():
        for t_launch, rows in tl_dict.items():
            for row in rows:
                ibox = row["box"]; grp = get_group(ibox, t_launch)
                if grp is None:
                    continue
                for c in _DCRE_BAR_COLS:
                    v = row.get(c, np.nan)
                    if np.isfinite(v):
                        data[pname][ibox][grp][c].append(v)

    n_traj = {pname: {ibox: {g: sum(1 for tl in all_launch_times
                                     if get_group(ibox, tl) == g)
                              for g in groups}
                      for ibox in range(4)}
              for pname in ml_traj_dcre}

    n_pairs_ml = len(ml_traj_dcre)
    if n_pairs_ml == 0:
        print("  No ML traj dCRE data — skipping Figure 4.")
        return

    fig, axes = plt.subplots(4, n_pairs_ml, figsize=(7 * n_pairs_ml, 18),
                             sharey="row")
    if n_pairs_ml == 1:
        axes = axes[:, np.newaxis]
    bar_width = 0.35
    x = np.arange(len(_DCRE_BAR_COLS))

    for ibox in range(4):
        for ip, pname in enumerate(ml_traj_dcre):
            ax = axes[ibox, ip]
            for ig, grp in enumerate(groups):
                d = data[pname][ibox][grp]
                means = [np.mean(d[c]) if d[c] else np.nan for c in _DCRE_BAR_COLS]
                stds  = [sem_list(d[c])  if d[c] else np.nan for c in _DCRE_BAR_COLS]
                offset = (ig - 0.5) * bar_width
                n = n_traj[pname][ibox][grp]
                ax.bar(x + offset, means, bar_width, yerr=stds, capsize=3,
                       error_kw={"linewidth": 1.0},
                       color=_LWP_GROUP_COLORS[grp], alpha=0.85,
                       label=f"{grp} LWP (n={n})")
            ax.axhline(0, color="k", lw=0.8)
            if ibox == 0:
                ax.set_title(pname, fontsize=10)
            ax.set_ylabel(r"W m$^{-2}$" if ip == 0 else "")
            ax.set_xticks(x)
            ax.set_xticklabels(_DCRE_BAR_LABELS, rotation=30, ha="right", fontsize=8)
            ax.grid(True, alpha=0.3, axis="y")
            if ip == n_pairs_ml - 1:
                ax2 = ax.twinx()
                ax2.set_ylabel(BOX_LABELS[ibox], color=BOX_COLORS[ibox],
                               fontsize=11, fontweight="bold", rotation=270, labelpad=15)
                ax2.set_yticks([])
            if ibox == 0 and ip == 0:
                ax.legend(fontsize=8, frameon=False, loc="upper right")
            if ibox < 3:
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel("dCRE component")
        for ibox in range(4):
            med_val = medians[ibox]
            med_str = f"{med_val:.1f}" if np.isfinite(med_val) else "N/A"
            axes[ibox, 0].annotate(
                f"median LWP = {med_str} g m⁻²",
                xy=(0.02, 0.97), xycoords="axes fraction",
                ha="left", va="top", fontsize=7.5, color=BOX_COLORS[ibox])

    low_patch  = Patch(color=_LWP_GROUP_COLORS["low"],  alpha=0.85,
                       label="Low LWP (< median per box)")
    high_patch = Patch(color=_LWP_GROUP_COLORS["high"], alpha=0.85,
                       label="High LWP (≥ median per box)")
    fig.legend(handles=[low_patch, high_patch], loc="lower center",
               bbox_to_anchor=(0.5, -0.02), ncol=2, frameon=True, fontsize=10)
    fig.suptitle(
        "dCRE decomposition — low vs high LWP trajectories (pixel-wise, ML pairs)\n"
        f"Launched from east boundary ({LAUNCHED_TRACE_CASE} wind)\n"
        "Error bars = ±1 SE",
        fontsize=12)
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])
    outname = os.path.join(OUTDIR, "dCRE_LWP_stratified_barplot.png")
    fig.savefig(outname, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outname}")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 5: Domain-averaged dCRE bar plot — all 4 pairs
# ══════════════════════════════════════════════════════════════════════════════

def plot_domain_dcre_barplot(domain_decomp_ts):
    pairs_with_data = [p for p, df in domain_decomp_ts.items() if not df.empty]
    if not pairs_with_data:
        print("  No domain dCRE data — skipping Figure 5.")
        return

    fig, axes = plt.subplots(1, len(pairs_with_data),
                             figsize=(6.5 * len(pairs_with_data), 5), sharey=True)
    if len(pairs_with_data) == 1:
        axes = [axes]

    x = np.arange(len(_DCRE_BAR_COLS)); bar_width = 0.55
    comp_colors = [s[1] for s in _DCRE_COMP_SPECS]

    for ip, pname in enumerate(pairs_with_data):
        ax = axes[ip]
        df = domain_decomp_ts[pname]
        means = [df[c].mean() if c in df else np.nan for c in _DCRE_BAR_COLS]
        stds  = [sem_series(df[c]) if c in df else np.nan for c in _DCRE_BAR_COLS]
        for xi, (m, s, col) in enumerate(zip(means, stds, comp_colors)):
            ax.bar(xi, m, bar_width, yerr=s, capsize=5,
                   error_kw={"linewidth": 1.3},
                   color=col,
                   hatch=_PAIR_HATCH.get(pname, ""),
                   alpha=_PAIR_ALPHA.get(pname, 0.85))
        ax.axhline(0, color="k", lw=0.8)
        scheme = DECOMP_PAIRS[pname]["scheme"]
        ax.set_title(f"{pname}\n[{scheme}]", fontsize=11)
        ax.set_ylabel(r"W m$^{-2}$" if ip == 0 else "")
        ax.set_xticks(x)
        ax.set_xticklabels(_DCRE_BAR_LABELS, rotation=30, ha="right", fontsize=9)
        ax.grid(True, alpha=0.3, axis="y")

    comp_handles = [Patch(color=s[1], alpha=0.85, label=s[3]) for s in _DCRE_COMP_SPECS]
    fig.legend(handles=comp_handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.10), ncol=3, frameon=True,
               fontsize=9, title="dCRE component")
    fig.suptitle(
        "Domain-averaged dCRE decomposition — pixel-wise  |  all 4 pairs\n"
        "Solid = ML, hatched = ARG  |  Error bars = ±1 SE",
        fontsize=13)
    fig.tight_layout()
    outname = os.path.join(OUTDIR, "dCRE_domain_barplot_all4pairs.png")
    fig.savefig(outname, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {outname}")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 6: met-binned dCRE  (8 variables; w in cm/s; ≤2-digit x ticks)
# ══════════════════════════════════════════════════════════════════════════════

def build_sample_dataframe(decomp_ts, box_ts, pair_name):
    ctrl   = DECOMP_PAIRS[pair_name]["control"]
    df_dc  = decomp_ts[pair_name]
    if df_dc.empty:
        return pd.DataFrame()
    ctrl_met   = box_ts[ctrl]
    time_index = {t: i for i, t in enumerate(ctrl_met["times"])}
    dcre_cols  = ["DCRE_total", "DCRE_CF", "DCRE_A_Nc",
                  "DCRE_A_LWP", "DCRE_A_re", "DCRE_A_cov"]
    met_vars   = [v for v, *_ in _MET_VAR_DEFS]
    rows = []
    for t, row in df_dc.iterrows():
        it   = time_index.get(t)
        if it is None:
            continue
        ibox = int(row.get("box", -1))
        if not (0 <= ibox < 4):
            continue
        sample = {}
        for col in dcre_cols:
            v = row.get(col, np.nan)
            sample[col] = float(v) if np.isfinite(float(v)) else np.nan
        for var in met_vars:
            arr = ctrl_met.get(var)
            sample[var] = float(arr[it, ibox]) if arr is not None else np.nan
        rows.append(sample)
    return pd.DataFrame(rows)


def _make_equal_freq_bins(series, n_bins):
    valid = series.dropna()
    if len(valid) < n_bins * 2:
        return pd.array([np.nan] * len(series), dtype="Float64"), np.array([])
    edges = np.nanpercentile(valid.values, np.linspace(0, 100, n_bins + 1))
    edges = np.unique(edges)
    if len(edges) < 3:
        return pd.array([np.nan] * len(series), dtype="Float64"), np.array([])
    midpoints = 0.5 * (edges[:-1] + edges[1:])
    idx = np.digitize(series.values, edges[:-1], right=False) - 1
    idx = np.clip(idx, 0, len(midpoints) - 1).astype(float)
    idx[series.isna().values] = np.nan
    return idx, midpoints


def _fmt_tick(v):
    """Format tick label to ≤2 significant digits, compact notation."""
    if not np.isfinite(v):
        return ""
    s = f"{v:.2g}"
    # avoid scientific notation for small-range axes
    return s


def plot_dCRE_met_binned_stacked(sample_df, pair_name):
    if sample_df.empty:
        print(f"  No samples for {pair_name} — skipping Figure 6.")
        return

    n_vars  = len(_MET_VAR_DEFS)
    n_rows  = (n_vars + 1) // 2   # ceil(8/2) = 4
    fig, axes = plt.subplots(n_rows, 2, figsize=(14, n_rows * 4.8))
    axes_flat = axes.flatten()
    bar_width = 0.72
    n_total   = len(sample_df)

    for iv, (var, var_name, units) in enumerate(_MET_VAR_DEFS):
        ax = axes_flat[iv]
        idx, mids = _make_equal_freq_bins(sample_df[var], N_MET_BINS)
        n_actual  = len(mids)
        if n_actual == 0:
            ax.set_title(f"{var_name}\n(no valid data)", fontsize=10)
            continue

        x_bins  = np.arange(n_actual)
        df_work = sample_df.copy()
        df_work["_bin"] = idx
        grp = df_work.dropna(subset=["_bin"]).groupby("_bin", observed=True)

        pos_bot = np.zeros(n_actual); neg_bot = np.zeros(n_actual)
        for comp, color, clabel in _STACK_SPECS:
            vals  = np.array([grp[comp].mean().get(float(b), np.nan)
                              for b in range(n_actual)])
            pos_v = np.where(np.isfinite(vals) & (vals > 0), vals, 0.0)
            neg_v = np.where(np.isfinite(vals) & (vals < 0), vals, 0.0)
            ax.bar(x_bins, pos_v, bar_width, bottom=pos_bot,
                   color=color, alpha=0.85, label=clabel)
            ax.bar(x_bins, neg_v, bar_width, bottom=neg_bot,
                   color=color, alpha=0.85)
            pos_bot += pos_v; neg_bot += neg_v

        tot_mean = np.array([grp["DCRE_total"].mean().get(float(b), np.nan)
                             for b in range(n_actual)])
        tot_sem  = np.array([
            (grp["DCRE_total"].std().get(float(b), np.nan) /
             np.sqrt(max(int(grp.size().get(float(b), 0)), 1)))
            for b in range(n_actual)])
        ax.errorbar(x_bins, tot_mean, yerr=tot_sem,
                    fmt="D", color="black", ms=5, capsize=4,
                    linewidth=1.2, zorder=6,
                    label=r"$\Delta$CRE$_{\rm total}$ ±SE")

        n_per_bin = np.array([int(grp.size().get(float(b), 0))
                              for b in range(n_actual)])
        y_low = ax.get_ylim()[0]
        for xi, nn in enumerate(n_per_bin):
            ax.text(xi, y_low, f"n={nn}", ha="center", va="bottom",
                    fontsize=5.5, color="dimgray")

        ax.axhline(0, color="k", linewidth=0.8, zorder=2)
        ax.set_ylabel(r"W m$^{-2}$", fontsize=9)
        ax.set_xlabel(f"{var_name} {units}".strip(), fontsize=9)
        ax.set_xticks(x_bins)
        # ≤2 significant digits; vertical velocity in cm/s already converted
        ax.set_xticklabels([_fmt_tick(m) for m in mids],
                           rotation=35, ha="right", fontsize=7.5)
        ax.set_title(var_name, fontsize=10)
        ax.grid(True, alpha=0.3, axis="y", zorder=0)

    # legend panel (last slot if odd number of variables)
    if n_vars < len(axes_flat):
        ax_leg = axes_flat[n_vars]
        ax_leg.axis("off")
        handles = [Patch(color=c, alpha=0.85, label=lbl)
                   for _, c, lbl in _STACK_SPECS]
        handles.append(
            Line2D([0], [0], marker="D", color="black", ms=5, linewidth=1.2,
                   linestyle="-", label=r"$\Delta$CRE$_{\rm total}$ ±SE"))
        ax_leg.legend(handles=handles, loc="center", fontsize=10,
                      frameon=True, title="dCRE component")
    else:
        # add legend to last panel
        axes_flat[-1].legend(fontsize=7, frameon=True, loc="best")

    ctrl = DECOMP_PAIRS[pair_name]["control"]
    scheme = DECOMP_PAIRS[pair_name]["scheme"]
    fig.suptitle(
        f"dCRE components vs. met regime — {pair_name} [{scheme}] (pixel-wise)\n"
        f"All 4 Lagrangian boxes  |  n = {n_total} box-timestep samples  |  "
        f"met from {ctrl} (control)  |  10 equal-freq. bins\n"
        r"Vertical velocity in cm s$^{-1}$; x-axis ≤2 sig. digits",
        fontsize=10)
    fig.tight_layout()

    stub = pair_name.replace("-", "_")
    for ext in ("pdf", "png"):
        outname = os.path.join(OUTDIR, f"dCRE_met_binned_{stub}.{ext}")
        fig.savefig(outname, dpi=200, bbox_inches="tight")
        print(f"  Saved {outname}")
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
# Figure 7: pixel-wise vs mean-based comparison  (all 4 pairs)
# ══════════════════════════════════════════════════════════════════════════════

def plot_decomp_comparison_barplot(box_decomp_pw, domain_decomp_pw,
                                   box_decomp_mb, domain_decomp_mb,
                                   pair_name):
    df_pw_box = box_decomp_pw.get(pair_name, pd.DataFrame())
    df_mb_box = box_decomp_mb.get(pair_name, pd.DataFrame())
    df_pw_dom = domain_decomp_pw.get(pair_name, pd.DataFrame())
    df_mb_dom = domain_decomp_mb.get(pair_name, pd.DataFrame())

    fig, axes = plt.subplots(3, 2, figsize=(14, 13))
    axes_flat = axes.flatten()
    bar_width = 0.35
    x = np.arange(len(_DCRE_BAR_COLS))

    methods = [
        ("pixel-wise", df_pw_box, df_pw_dom, _METHOD_COLORS["pixel-wise"]),
        ("mean-based", df_mb_box, df_mb_dom, _METHOD_COLORS["mean-based"]),
    ]

    for ibox in range(4):
        ax = axes_flat[ibox]
        for im, (label, df_box, _, color) in enumerate(methods):
            if df_box.empty:
                continue
            sub = df_box[df_box["box"] == ibox]
            if sub.empty:
                continue
            means  = [sub[c].mean() if c in sub else np.nan for c in _DCRE_BAR_COLS]
            stds   = [sem_series(sub[c]) if c in sub else np.nan for c in _DCRE_BAR_COLS]
            offset = (im - 0.5) * bar_width
            ax.bar(x + offset, means, bar_width, yerr=stds, capsize=4,
                   error_kw={"linewidth": 1.2}, color=color, alpha=0.85, label=label)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_title(BOX_LABELS[ibox], color=BOX_COLORS[ibox],
                     fontsize=12, fontweight="bold")
        ax.set_ylabel(r"W m$^{-2}$")
        ax.set_xticks(x)
        ax.set_xticklabels(_DCRE_BAR_LABELS, rotation=30, ha="right", fontsize=9)
        ax.grid(True, alpha=0.3, axis="y")
        if ibox == 0:
            ax.legend(fontsize=9, frameon=True, loc="best")

    ax_dom = axes_flat[4]
    for im, (label, _, df_dom, color) in enumerate(methods):
        if df_dom.empty:
            continue
        means  = [df_dom[c].mean() if c in df_dom else np.nan for c in _DCRE_BAR_COLS]
        stds   = [sem_series(df_dom[c]) if c in df_dom else np.nan for c in _DCRE_BAR_COLS]
        offset = (im - 0.5) * bar_width
        ax_dom.bar(x + offset, means, bar_width, yerr=stds, capsize=4,
                   error_kw={"linewidth": 1.2}, color=color, alpha=0.85, label=label)
    ax_dom.axhline(0, color="k", lw=0.8)
    ax_dom.set_title("Whole domain", fontsize=12, fontweight="bold")
    ax_dom.set_ylabel(r"W m$^{-2}$")
    ax_dom.set_xticks(x)
    ax_dom.set_xticklabels(_DCRE_BAR_LABELS, rotation=30, ha="right", fontsize=9)
    ax_dom.grid(True, alpha=0.3, axis="y")

    ax_leg = axes_flat[5]; ax_leg.axis("off")
    mp = [Patch(color=_METHOD_COLORS["pixel-wise"], alpha=0.85,
                label="Pixel-wise (original)\nColumn-by-column"),
          Patch(color=_METHOD_COLORS["mean-based"], alpha=0.85,
                label="Mean-based (new)\nBox-mean scalars; LWP-weighted re")]
    cl = [Line2D([0], [0], color=s[1], lw=2.5, label=s[3]) for s in _DCRE_COMP_SPECS]
    ax_leg.legend(handles=mp + cl, loc="center", fontsize=9, frameon=True,
                  title="Method  /  component")
    ax_leg.text(0.5, 0.06,
                "Error bars = ±1 SE  |  Mean-based: DCRE_CF + DCRE_A ≡ total (exact)",
                ha="center", va="bottom", fontsize=8,
                transform=ax_leg.transAxes, style="italic", color="0.4")

    scheme = DECOMP_PAIRS[pair_name]["scheme"]
    fig.suptitle(
        f"dCRE decomposition method comparison — {pair_name} [{scheme}]\n"
        "Pixel-wise (blue) vs mean-based (red)",
        fontsize=13, y=1.01)
    fig.tight_layout()
    stub = pair_name.replace("-", "_")
    for ext in ("png", "pdf"):
        outname = os.path.join(OUTDIR, f"dCRE_method_comparison_{stub}.{ext}")
        fig.savefig(outname, dpi=200, bbox_inches="tight")
        print(f"  Saved {outname}")
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
# Figure 8 (NEW): domain-mean ML vs ARG side-by-side comparison
# ══════════════════════════════════════════════════════════════════════════════

def plot_ML_ARG_domain_comparison(domain_decomp_pw, domain_decomp_mb):
    """2-row × 2-col figure.

    Rows    : org background (top)  |  mid background (bottom)
    Columns : pixel-wise (left)     |  mean-based (right)
    Within each panel: grouped bars, one group per dCRE component.
      Solid fill = ML,  hatched '///' = ARG.
    Error bars = ±1 SE across time steps.
    """
    backgrounds = ["org", "mid"]
    methods     = [
        ("pixel-wise", domain_decomp_pw),
        ("mean-based", domain_decomp_mb),
    ]

    fig, axes = plt.subplots(len(backgrounds), 2,
                             figsize=(16, 5.5 * len(backgrounds)),
                             squeeze=False)
    plt.subplots_adjust(hspace=0.38, wspace=0.25)

    bar_width = 0.35
    x = np.arange(len(_DCRE_BAR_COLS))

    for irow, bg in enumerate(backgrounds):
        grp_info = BACKGROUND_GROUPS[bg]
        ml_pname  = grp_info["ML"]
        arg_pname = grp_info["ARG"]

        for icol, (method_label, decomp_dict) in enumerate(methods):
            ax = axes[irow, icol]

            df_ml  = decomp_dict.get(ml_pname,  pd.DataFrame())
            df_arg = decomp_dict.get(arg_pname, pd.DataFrame())

            for df, label, hatch, alpha in [
                (df_ml,  "ML",  "",    0.85),
                (df_arg, "ARG", "///", 0.55),
            ]:
                if df.empty:
                    continue
                means  = [df[c].mean() if c in df else np.nan for c in _DCRE_BAR_COLS]
                stds   = [sem_series(df[c]) if c in df else np.nan for c in _DCRE_BAR_COLS]
                offset = (-0.5 if label == "ML" else 0.5) * bar_width
                color  = "#4393c3" if bg == "org" else "#d6604d"
                ax.bar(x + offset, means, bar_width,
                       yerr=stds, capsize=4, error_kw={"linewidth": 1.1},
                       color=color, hatch=hatch, alpha=alpha, label=label)

            ax.axhline(0, color="k", lw=0.8)
            ax.set_xticks(x)
            ax.set_xticklabels(_DCRE_BAR_LABELS, rotation=30, ha="right", fontsize=9)
            ax.set_ylabel(r"W m$^{-2}$", fontsize=10)
            ax.set_title(
                f"{grp_info['label']}  |  {method_label}\n"
                f"ML: {ml_pname}   ARG: {arg_pname}",
                fontsize=9)
            ax.grid(True, alpha=0.3, axis="y")
            if irow == 0 and icol == 0:
                ax.legend(fontsize=9, frameon=True)

    fig.suptitle(
        "Domain-mean dCRE decomposition: ML vs ARG\n"
        "Solid fill = ML  |  Hatched = ARG  |  Error bars = ±1 SE across time",
        fontsize=14)

    for ext in ("png", "pdf"):
        outname = os.path.join(OUTDIR, f"dCRE_ML_ARG_domain_comparison.{ext}")
        fig.savefig(outname, dpi=200, bbox_inches="tight")
        print(f"  Saved {outname}")
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
# Dispersion (β–Nc) consistency check  (NEW)
#
# Recovers the gamma-PSD shape factor   β = r_e / r_v   with
#     r_v = (3·LWP / (4π ρ_w N_c))^(1/3)        (column volume-mean radius)
# from the model's own RE_QC, LWP_TOT and column-integrated N_c (QNDROP).  β is
# the same "beyond-Twomey r_e" factor that DCRE_A_re isolates; here we expose its
# N_c dependence directly so it can be compared against the observed dispersion
# effect (Liu & Daum 2002, Nature 419, 580; Wang et al. 2025, Nat. Commun.).
#
# For τ ∝ β^{-1} LWP^{a} N_c^{1/3}, the dispersion modifies the Twomey number term
# by a factor (1 − 3b), with b ≡ d ln β / d ln N_c.  The fractional offset to the
# Twomey term is therefore  (−3b)·100%.  Wang et al. report |b| ≈ 0.024 (≈7%).
# Note the cloud fraction cancels in r_v (LWP and N_c both scale with cf), so
# grid-box column LWP and N_c are used directly.  For a gamma PSD ε = 1/√(μ+1)
# and β = (1+2ε²)^{2/3}/(1+ε²)^{1/3}; Morrison sets μ via the Martin et al. (1994)
# N_c relation (module_mp_morr_two_moment.F), identical in the ML and ARG runs.
# ══════════════════════════════════════════════════════════════════════════════

_EPS_GRID  = np.linspace(0.0, 1.2, 600)
_BETA_GRID = ((1.0 + 2.0 * _EPS_GRID**2)**(2.0/3.0) /
              (1.0 + _EPS_GRID**2)**(1.0/3.0))


def beta_to_epsilon(beta):
    """Invert β = (1+2ε²)^(2/3)/(1+ε²)^(1/3) for the relative dispersion ε.

    _BETA_GRID is monotonically increasing in ε so a 1-D interpolation suffices.
    β < 1 (unphysical for a gamma PSD) returns NaN."""
    beta = np.asarray(beta, dtype=float)
    eps  = np.interp(np.clip(beta, _BETA_GRID[0], _BETA_GRID[-1]),
                     _BETA_GRID, _EPS_GRID)
    return np.where(beta >= _BETA_GRID[0] - 1e-6, eps, np.nan)


def _beta_nc_pixels(d, cf_min=CF_MIN):
    """In-cloud (β, N_c[m⁻²], ε) pixel arrays from an extract_fields_for_dcre dict."""
    cld  = d["cf"] > cf_min
    LWP  = d["LWP"]; Nc = d["Nc"]; re_m = d["re"] * 1.0e-6   # µm → m
    ok   = (cld & np.isfinite(re_m) & (re_m > 1.0e-7) &
            np.isfinite(LWP) & (LWP > EPS) &
            np.isfinite(Nc)  & (Nc  > EPS))
    if not np.any(ok):
        empty = np.array([])
        return empty, empty, empty
    LWP = LWP[ok]; Nc = Nc[ok]; re_m = re_m[ok]
    r_v  = (3.0 * LWP / (4.0 * np.pi * RHO_W * Nc))**(1.0/3.0)
    beta = re_m / r_v
    eps  = beta_to_epsilon(beta)
    return beta, Nc, eps


def compute_beta_nc_consistency(file_by_case_time, max_samples_per_step=1500):
    """Per-case β–N_c statistics for the dispersion consistency check.

    Reads each case's wrfout once; accumulates exact OLS regression sums of
    ln β on ln N_c (memory-independent of record length) plus a capped random
    subsample of (ln N_c, ln β) per timestep for the scatter plot."""
    rng     = np.random.default_rng(0)
    stats   = {}
    samples = {}
    for case_name in case_dirs:
        S = dict(n=0, Sx=0.0, Sy=0.0, Sxx=0.0, Sxy=0.0, Sbeta=0.0, Seps=0.0)
        lnNc_s = []; lnbeta_s = []
        for t, f in sorted(file_by_case_time[case_name].items()):
            try:
                with xr.open_dataset(f, decode_times=False) as ds:
                    d = extract_fields_for_dcre(squeeze_time(ds))
            except Exception as e:
                print(f"    [β–Nc] {case_name} {t}: {e}")
                continue
            beta, Nc, eps = _beta_nc_pixels(d)
            if beta.size == 0:
                continue
            good = np.isfinite(beta) & (beta > 0) & np.isfinite(eps)
            beta = beta[good]; Nc = Nc[good]; eps = eps[good]
            if beta.size == 0:
                continue
            x = np.log(Nc); y = np.log(beta)
            S["n"]     += int(x.size)
            S["Sx"]    += float(x.sum());      S["Sy"]  += float(y.sum())
            S["Sxx"]   += float((x*x).sum());  S["Sxy"] += float((x*y).sum())
            S["Sbeta"] += float(beta.sum());   S["Seps"] += float(eps.sum())
            if x.size > max_samples_per_step:
                sel = rng.choice(x.size, max_samples_per_step, replace=False)
                x = x[sel]; y = y[sel]
            lnNc_s.append(x); lnbeta_s.append(y)
        stats[case_name]   = S
        samples[case_name] = (
            np.concatenate(lnNc_s)   if lnNc_s   else np.array([]),
            np.concatenate(lnbeta_s) if lnbeta_s else np.array([]))
        if S["n"] > 0:
            print(f"  {case_name}: β–Nc n={S['n']}  "
                  f"<β>={S['Sbeta']/S['n']:.3f}  <ε>={S['Seps']/S['n']:.3f}")
    return stats, samples


def _ols_slope(S):
    """OLS slope b, intercept, and mean(y) from accumulated regression sums."""
    n = S["n"]
    if n < 10:
        return np.nan, np.nan, np.nan
    mx = S["Sx"] / n; my = S["Sy"] / n
    cov = S["Sxy"] / n - mx * my
    var = S["Sxx"] / n - mx * mx
    b   = cov / var if var > EPS else np.nan
    return b, (my - b * mx if np.isfinite(b) else np.nan), my


def report_dispersion_consistency(stats):
    """Build per-case and per-pair β–Nc tables, print a summary, save CSVs."""
    rows = []
    for case_name, S in stats.items():
        b, _, _ = _ols_slope(S)
        n = S["n"]
        rows.append(dict(
            case=case_name, n=n,
            beta_mean=(S["Sbeta"]/n if n else np.nan),
            eps_mean =(S["Seps"]/n  if n else np.nan),
            lnNc_mean=(S["Sx"]/n    if n else np.nan),
            lnbeta_mean=(S["Sy"]/n  if n else np.nan),
            b_dlnbeta_dlnNc=b,
            twomey_offset_pct=-3.0 * b * 100.0))
    df_case = pd.DataFrame(rows).set_index("case")

    pair_rows = []
    for pname, pair in DECOMP_PAIRS.items():
        c = stats.get(pair["control"]); p = stats.get(pair["perturbed"])
        if not c or not p or c["n"] < 10 or p["n"] < 10:
            continue
        dlnbeta = p["Sy"]/p["n"] - c["Sy"]/c["n"]
        dlnNc   = p["Sx"]/p["n"] - c["Sx"]/c["n"]
        b_pair  = dlnbeta / dlnNc if abs(dlnNc) > EPS else np.nan
        pair_rows.append(dict(
            pair=pname, scheme=pair["scheme"], background=pair["background"],
            dlnbeta=dlnbeta, dlnNc=dlnNc, b_pair=b_pair,
            twomey_offset_pct=-3.0 * b_pair * 100.0))
    df_pair = (pd.DataFrame(pair_rows).set_index("pair")
               if pair_rows else pd.DataFrame())

    print("\n  ── Dispersion (β–Nc) consistency ─────────────────────────────────")
    print("  b = dlnβ/dlnNc ;  Twomey offset = −3b·100%  "
          "(Wang et al. 2025: |b|≈0.024, ≈7%)")
    with pd.option_context("display.float_format", lambda v: f"{v:.4g}"):
        print(df_case[["n", "beta_mean", "eps_mean",
                       "b_dlnbeta_dlnNc", "twomey_offset_pct"]].to_string())
        if not df_pair.empty:
            print("\n  Per-pair (perturbed − control; sign matches DCRE_A_re):")
            print(df_pair[["scheme", "background", "dlnbeta", "dlnNc",
                           "b_pair", "twomey_offset_pct"]].to_string())

    p1 = os.path.join(OUTDIR, "dispersion_beta_Nc_by_case.csv")
    df_case.to_csv(p1); print(f"  Saved {p1}")
    if not df_pair.empty:
        p2 = os.path.join(OUTDIR, "dispersion_beta_Nc_by_pair.csv")
        df_pair.to_csv(p2); print(f"  Saved {p2}")
    return df_case, df_pair


def plot_dispersion_consistency(stats, samples, df_case, df_pair):
    """Figure 9: (a) per-case ln β vs ln N_c scatter + OLS fit,
                 (b) Twomey offset −3b per case and per pair vs Wang et al."""
    # ── 9a: per-case scatter + fit ────────────────────────────────────────────
    order = ML_CASES + ARG_CASES
    _nr = -(-len(order) // 4)
    fig, axes = plt.subplots(_nr, 4, figsize=(20, 4.5 * _nr), sharey=True,
                             squeeze=False)
    for ax, case_name in zip(axes.flatten(), order):
        x, y = samples.get(case_name, (np.array([]), np.array([])))
        if x.size:
            ax.scatter(x, y, s=3, alpha=0.12, color="#4393c3", edgecolors="none")
            b, b0, _ = _ols_slope(stats[case_name])
            if np.isfinite(b):
                xs = np.linspace(np.nanpercentile(x, 1),
                                 np.nanpercentile(x, 99), 50)
                ax.plot(xs, b0 + b * xs, color="#d6604d", lw=2.0)
                ax.annotate(f"b={b:.3f}\n−3b={-3*b*100:.1f}%",
                            xy=(0.04, 0.96), xycoords="axes fraction",
                            ha="left", va="top", fontsize=9,
                            bbox=dict(boxstyle="round", fc="white", alpha=0.8))
        ax.set_title(case_name, fontsize=10)
        ax.grid(True, alpha=0.3)
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$\ln N_c$ (column, m$^{-2}$)")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\ln\beta$,  $\beta=r_e/r_v$")
    fig.suptitle(
        "Dispersion consistency: PSD shape factor β vs droplet number\n"
        "Top = ML, bottom = ARG  |  red = OLS fit  |  "
        r"Twomey offset $=-3\,d\ln\beta/d\ln N_c$  (Wang et al. 2025 $\approx$ 7%)",
        fontsize=12)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        outn = os.path.join(OUTDIR, f"dispersion_beta_Nc_scatter.{ext}")
        fig.savefig(outn, dpi=200, bbox_inches="tight"); print(f"  Saved {outn}")
    plt.close(fig)

    # ── 9b: Twomey offset bars ────────────────────────────────────────────────
    fig, (axc, axp) = plt.subplots(1, 2, figsize=(15, 5))
    cases  = list(df_case.index)
    off_c  = df_case["twomey_offset_pct"].values
    colors = ["#4393c3" if c in ML_CASES else "#d6604d" for c in cases]
    axc.bar(np.arange(len(cases)), off_c, color=colors, alpha=0.85)
    axc.axhline(7.2, color="k", ls="--", lw=1.2, label="Wang et al. ≈ 7%")
    axc.axhline(0, color="k", lw=0.8)
    axc.set_xticks(np.arange(len(cases)))
    axc.set_xticklabels(cases, rotation=35, ha="right", fontsize=8)
    axc.set_ylabel(r"Twomey offset $-3b$ (%)")
    axc.set_title("Per-case dispersion offset")
    axc.legend(fontsize=8); axc.grid(True, alpha=0.3, axis="y")

    if not df_pair.empty:
        pairs = list(df_pair.index)
        offp  = df_pair["twomey_offset_pct"].values
        pcol  = ["#4393c3" if df_pair.loc[p, "scheme"] == "ML" else "#d6604d"
                 for p in pairs]
        axp.bar(np.arange(len(pairs)), offp, color=pcol, alpha=0.85)
        axp.axhline(7.2, color="k", ls="--", lw=1.2, label="Wang et al. ≈ 7%")
        axp.axhline(0, color="k", lw=0.8)
        axp.set_xticks(np.arange(len(pairs)))
        axp.set_xticklabels(pairs, rotation=20, ha="right", fontsize=8)
        axp.set_ylabel(r"Twomey offset $-3b$ (%)")
        axp.set_title("Per-pair (perturbed − control)")
        axp.legend(fontsize=8); axp.grid(True, alpha=0.3, axis="y")
    else:
        axp.axis("off")
    fig.suptitle("Modelled dispersion effect relative to the Twomey number term",
                 fontsize=12)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        outn = os.path.join(OUTDIR, f"dispersion_twomey_offset.{ext}")
        fig.savefig(outn, dpi=200, bbox_inches="tight"); print(f"  Saved {outn}")
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

# ══════════════════════════════════════════════════════════════════════════════
# Figure 10 (v22): ML vs ARG domain comparison, per hour and over all times
#
# Two panels, pixel-wise only, following replot_dcre_ml_arg_comparison.py.  The
# mean-based route leaves up to 1.45 W m-2 in an uninterpretable covariance term
# and returns NaN sub-terms for both ARG pairs, so it is not drawn.
#
# Whether the four ARG sub-terms are measurements depends on which ARG tree is
# in use, so it is DETECTED rather than assumed.  In WRF_dm, TAU_QC_TOT and
# RE_QC are identically zero, R_tau degenerates and the apportionment falls back
# to the analytic 1/3 : 2/3 Twomey split; those bars are then drawn grey and
# annotated so the figure cannot be misread.  In WRF_dm_v2 both fields are
# populated, R_tau is a real ratio, and the bars are drawn in the panel colour
# like any other measurement.
# ══════════════════════════════════════════════════════════════════════════════

# sub-terms whose ARG values are the analytic fallback when R_tau degenerates
_ARG_FALLBACK_COLS = {"DCRE_A_Nc", "DCRE_A_LWP", "DCRE_A_re", "DCRE_A_cov"}


def arg_terms_are_fallback(df, tol=1e-3):
    """True if this record's R_tau carries no information.

    The apportionment is analytic exactly when tau does not respond, i.e. when
    R_tau is absent or pinned at 1.  Checking the data beats hard-coding the
    tree, because the same script now serves both WRF_dm and WRF_dm_v2.
    """
    if "R_tau" not in df:
        return True
    r = pd.to_numeric(df["R_tau"], errors="coerce").dropna()
    return len(r) == 0 or bool((r - 1.0).abs().max() < tol)

# (panel tag, background key, panel title, colour) -- pair names come from
# BACKGROUND_GROUPS so the two stay in step.
_FIG10_PANELS = [("a", "org",  "org background",  "#4393c3"),
                 ("b", "mid",  "mid background",  "#d6604d"),
                 ("c", "mid2", "mid2 background", "#4a3aa7")]


def _slice_window(df, t0, t1):
    """Rows with t0 <= time < t1.  t0/t1 None means unbounded."""
    if df is None or df.empty:
        return pd.DataFrame()
    idx = pd.to_datetime(df.index)
    keep = np.ones(len(df), dtype=bool)
    if t0 is not None:
        keep &= (idx >= t0)
    if t1 is not None:
        keep &= (idx < t1)
    return df.loc[keep]


def plot_ML_ARG_domain_comparison_window(domain_decomp_pw, t0, t1,
                                         stem, window_label):
    """One two-panel ML-vs-ARG figure for the time window [t0, t1).

    Returns True if the figure was drawn, False if the window was too thin in
    any panel (the caller then skips that hour rather than shipping a figure
    whose error bars mean nothing).
    """
    panels = []
    for tag, bg, title, color in _FIG10_PANELS:
        grp = BACKGROUND_GROUPS[bg]
        df_ml  = _slice_window(domain_decomp_pw.get(grp["ML"]),  t0, t1)
        df_arg = _slice_window(domain_decomp_pw.get(grp["ARG"]), t0, t1)
        panels.append((tag, title, color, df_ml, df_arg))

    n_min = min(min(len(p[3]), len(p[4])) for p in panels)
    if n_min < 1:
        print(f"  {window_label}: no samples in one of the panels — skipped")
        return False

    fig, axes = plt.subplots(1, len(panels), figsize=(7.5 * len(panels), 6.2),
                             squeeze=False)
    axes = axes[0]
    plt.subplots_adjust(wspace=0.22)
    x  = np.arange(len(_DCRE_BAR_COLS))
    bw = 0.36

    any_fallback = False
    for ax, (tag, title, color, df_ml, df_arg) in zip(axes, panels):
        # detected per panel: WRF_dm degenerates, WRF_dm_v2 does not
        fb = arg_terms_are_fallback(df_arg)
        any_fallback = any_fallback or fb
        for df, lab, hatch, alpha, off in ((df_ml,  "ML-ACT",  "",    0.90, -0.5),
                                           (df_arg, "ARG-ACT", "///", 0.55,  0.5)):
            means, errs, cols = [], [], []
            for c in _DCRE_BAR_COLS:
                v = (pd.to_numeric(df[c], errors="coerce")
                     if c in df else pd.Series(dtype=float))
                means.append(v.mean() if len(v) else np.nan)
                errs.append(sem_series(v) if len(v) else np.nan)
                cols.append("#b0b0b0"
                            if (fb and lab.startswith("ARG")
                                and c in _ARG_FALLBACK_COLS)
                            else color)
            ax.bar(x + off * bw, means, bw, yerr=errs, capsize=4,
                   error_kw={"linewidth": 1.1},
                   color=cols, hatch=hatch, alpha=alpha,
                   edgecolor="k", linewidth=0.6, label=lab)

        ax.axhline(0, color="k", lw=0.9)
        ax.set_xticks(x)
        ax.set_xticklabels(_DCRE_BAR_LABELS, rotation=25, ha="right", fontsize=13)
        ax.set_ylabel(r"$\Delta$CRE  [W m$^{-2}$]", fontsize=15)
        ax.tick_params(axis="y", labelsize=13)
        ax.set_title(f"({tag}) {title}", fontsize=17, loc="left")
        ax.grid(True, alpha=0.3, axis="y")

        # shade the region occupied by the non-measured ARG sub-terms
        ax.axvspan(1.5, len(_DCRE_BAR_COLS) - 0.5, color="0.5", alpha=0.06, zorder=0)

        # Per-panel legend in THIS panel's colour: a single shared legend would
        # have to pick one hue and would misrepresent the other panel.
        _lh = [Patch(facecolor=color, edgecolor="k", alpha=0.90, label="ML-ACT"),
               Patch(facecolor=color, edgecolor="k", alpha=0.55, hatch="///",
                     label="ARG-ACT")]
        if fb:
            _lh.append(Patch(facecolor="#b0b0b0", edgecolor="k", alpha=0.55,
                             hatch="///", label="ARG-ACT (analytic fallback)"))
        ax.legend(handles=_lh, fontsize=11, frameon=True, loc="upper left")

        tot_ml  = (pd.to_numeric(df_ml["DCRE_total"], errors="coerce").mean()
                   if "DCRE_total" in df_ml else np.nan)
        tot_arg = (pd.to_numeric(df_arg["DCRE_total"], errors="coerce").mean()
                   if "DCRE_total" in df_arg else np.nan)
        ratio = (tot_arg / tot_ml) if (np.isfinite(tot_ml) and tot_ml != 0) else np.nan
        ax.text(0.015, 0.03,
                f"ML {tot_ml:+.2f}   ARG {tot_arg:+.2f} W m$^{{-2}}$"
                f"   (ARG/ML = {ratio:.2f}),  n={len(df_ml)}/{len(df_arg)} times",
                transform=ax.transAxes, fontsize=12, va="bottom",
                bbox=dict(fc="white", alpha=0.85, ec="0.7"))

    fig.text(0.5, -0.035,
             "Colour denotes the background (blue = org, red = mid); fill denotes "
             "the activation treatment."
             + ("  Grey bars are the analytic fallback, not a measurement."
                if any_fallback else
                "  All ARG sub-terms are measured: TAU_QC_TOT and RE_QC are "
                "populated in this ARG tree."),
             ha="center", fontsize=12)
    fig.suptitle("Domain-mean $\\Delta$CRE decomposition for the $3\\times$ sulfate "
                 f"perturbation: ML-ACT vs ARG-ACT (pixel-wise)\n{window_label}",
                 fontsize=18, y=1.02)

    for ext in ("png", "pdf"):
        out = os.path.join(OUTDIR, f"{stem}.{ext}")
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"  Saved {out}")
    plt.close(fig)
    return True


def plot_ML_ARG_domain_comparison_hourly(domain_decomp_pw):
    """One Figure-10 per hour, then one over every shared time."""
    times = []
    for pname in DECOMP_PAIRS:
        df = domain_decomp_pw.get(pname)
        if df is not None and not df.empty:
            times.extend(pd.to_datetime(df.index).tolist())
    if not times:
        print("  No domain decomposition rows — skipping Figure 10.")
        return
    tmin, tmax = min(times), max(times)

    hour = pd.Timestamp(tmin).floor("h")
    end  = pd.Timestamp(tmax)
    n_drawn = 0
    while hour <= end:
        nxt = hour + pd.Timedelta(minutes=HOURLY_BIN_MINUTES)
        # only draw an hour that is populated in every panel
        counts = []
        for _, bg, _, _ in _FIG10_PANELS:
            grp = BACKGROUND_GROUPS[bg]
            counts.append(len(_slice_window(domain_decomp_pw.get(grp["ML"]), hour, nxt)))
            counts.append(len(_slice_window(domain_decomp_pw.get(grp["ARG"]), hour, nxt)))
        if counts and min(counts) >= MIN_TIMES_PER_HOUR:
            label = (f"hour {hour:%H:%M}–{nxt:%H:%M} UTC on {hour:%Y-%m-%d}")
            if plot_ML_ARG_domain_comparison_window(
                    domain_decomp_pw, hour, nxt,
                    f"dCRE_ML_ARG_domain_comparison_HH{hour:%H}", label):
                n_drawn += 1
        else:
            print(f"  hour {hour:%H:%M}: only {min(counts) if counts else 0} time(s) "
                  f"in the thinnest panel (< {MIN_TIMES_PER_HOUR}) — skipped")
        hour = nxt

    # the final figure over every shared time
    label = (f"all shared times, {tmin:%Y-%m-%d %H:%M}–{tmax:%H:%M} UTC "
             f"(intersection across all 8 cases)")
    plot_ML_ARG_domain_comparison_window(
        domain_decomp_pw, None, None,
        "dCRE_ML_ARG_domain_comparison_ALL", label)
    print(f"  Figure 10: {n_drawn} hourly figure(s) + 1 all-times figure")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 11 (v22): half-hourly stacked-bar TIME SERIES of the domain-mean dCRE
#
# Same stacking as the met-binned Figure 6 -- components stack up and down about
# zero by sign, black diamonds are DCRE_total +/- SE within the bin -- but the
# x-axis is time in STACK_BIN_MINUTES bins rather than a met regime.
# ══════════════════════════════════════════════════════════════════════════════

def plot_dCRE_halfhour_stacked_timeseries(domain_decomp_pw, pair_name):
    df = domain_decomp_pw.get(pair_name)
    if df is None or df.empty:
        print(f"  No domain rows for {pair_name} — skipping Figure 11.")
        return

    df = df.copy()
    df.index = pd.to_datetime(df.index)
    bin_start = df.index.floor(f"{STACK_BIN_MINUTES}min")
    grp = df.groupby(bin_start)

    n_per = grp.size()
    keep  = n_per[n_per >= MIN_TIMES_PER_STACK].index
    if len(keep) == 0:
        print(f"  Every {STACK_BIN_MINUTES}-min bin for {pair_name} is below "
              f"{MIN_TIMES_PER_STACK} samples — skipping Figure 11.")
        return

    starts = pd.DatetimeIndex(sorted(keep))
    half   = pd.Timedelta(minutes=STACK_BIN_MINUTES) / 2
    centres = starts + half
    # bar width in matplotlib date units (days), with a small gap between bars
    width = (STACK_BIN_MINUTES / (24.0 * 60.0)) * 0.88

    fig, ax = plt.subplots(figsize=(max(11.0, 1.15 * len(starts) + 4.0), 6.4))

    pos_bot = np.zeros(len(starts))
    neg_bot = np.zeros(len(starts))
    for comp, color, clabel in _STACK_SPECS:
        if comp in df:
            vals = np.array([grp[comp].mean().get(s, np.nan) for s in starts],
                            dtype=float)
        else:
            vals = np.full(len(starts), np.nan)
        pos_v = np.where(np.isfinite(vals) & (vals > 0), vals, 0.0)
        neg_v = np.where(np.isfinite(vals) & (vals < 0), vals, 0.0)
        ax.bar(centres, pos_v, width=width, bottom=pos_bot,
               color=color, alpha=0.85, label=clabel, zorder=3)
        ax.bar(centres, neg_v, width=width, bottom=neg_bot,
               color=color, alpha=0.85, zorder=3)
        pos_bot += pos_v
        neg_bot += neg_v

    tot_mean = np.array([grp["DCRE_total"].mean().get(s, np.nan) for s in starts],
                        dtype=float)
    tot_sem = np.array(
        [(grp["DCRE_total"].std().get(s, np.nan) /
          np.sqrt(max(int(n_per.get(s, 0)), 1))) for s in starts], dtype=float)
    ax.errorbar(centres, tot_mean, yerr=tot_sem, fmt="D", color="black",
                ms=5, capsize=4, linewidth=1.2, zorder=6,
                label=r"$\Delta$CRE$_{\rm total}$ ±SE")

    y_low = ax.get_ylim()[0]
    for c, s in zip(centres, starts):
        ax.text(c, y_low, f"n={int(n_per.get(s, 0))}", ha="center", va="bottom",
                fontsize=6.5, color="dimgray", zorder=7)

    ax.axhline(0, color="k", linewidth=0.9, zorder=4)
    ax.set_ylabel(r"$\Delta$CRE  [W m$^{-2}$]", fontsize=14)
    ax.set_xlabel("time (UTC)", fontsize=14)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.set_xticks(centres)
    ax.set_xticklabels([f"{s:%H:%M}" for s in starts], rotation=35,
                       ha="right", fontsize=11)
    ax.tick_params(axis="y", labelsize=12)
    ax.grid(True, alpha=0.3, axis="y", zorder=0)
    ax.legend(fontsize=10, frameon=True, ncol=2, loc="best")

    scheme = DECOMP_PAIRS[pair_name]["scheme"]
    bg     = DECOMP_PAIRS[pair_name]["background"]
    fig.suptitle(
        f"Domain-mean $\\Delta$CRE decomposition vs. time — {pair_name} "
        f"[{scheme}, {bg} background] (pixel-wise)\n"
        f"one stacked bar per {STACK_BIN_MINUTES} min  |  "
        f"{len(df)} output times over "
        f"{df.index.min():%Y-%m-%d %H:%M}–{df.index.max():%H:%M} UTC  |  "
        f"bar label n = output times in the bin",
        fontsize=12)
    fig.tight_layout()

    stub = pair_name.replace("-", "_")
    for ext in ("png", "pdf"):
        out = os.path.join(OUTDIR, f"dCRE_timeseries_stacked_{stub}.{ext}")
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"  Saved {out}")
    plt.close(fig)


def main():
    # ── load files ────────────────────────────────────────────────────────────
    all_files = {}
    for case_name, case_dir in case_dirs.items():
        files = get_case_files(case_dir)
        if not files:
            raise RuntimeError(f"No wrfout files for {case_name}: {case_dir}")
        all_files[case_name] = files
        print(f"{case_name}: {len(files)} wrfout files")

    # ── restrict every case to the time steps SHARED by all ML + ARG cases ────
    shared_times = set.intersection(
        *[set(t for t, _ in files) for files in all_files.values()])
    if not shared_times:
        raise RuntimeError("No time steps shared across all ML and ARG cases.")
    if TIME_END is not None:
        n_late = sum(1 for t in shared_times if t > TIME_END)
        shared_times = {t for t in shared_times if t <= TIME_END}
        if not shared_times:
            raise RuntimeError(f"No shared time steps at or before {TIME_END}.")
        print(f"Analysis window ends {TIME_END:%Y-%m-%d %H:%M} "
              f"(WRF_TEND); dropped {n_late} later shared step(s)")
    for case_name in all_files:
        all_files[case_name] = [(t, f) for (t, f) in all_files[case_name]
                                if t in shared_times]
    s_sorted = sorted(shared_times)
    print(f"\nShared time steps across all {len(all_files)} cases: "
          f"{len(shared_times)} "
          f"({s_sorted[0]:%Y-%m-%d %H:%M} .. {s_sorted[-1]:%H:%M})")

    # ── trace boxes for all 8 cases ───────────────────────────────────────────
    print("\nTracing Lagrangian boxes for all cases ...")
    case_results = {}; file_by_case_time = {}
    for case_name, files in all_files.items():
        states, file_by_time, grid = trace_one_case(case_name, files)
        case_results[case_name]      = (states, file_by_time, grid)
        file_by_case_time[case_name] = file_by_time

    _, _, ref_grid = case_results[LAUNCHED_TRACE_CASE]
    ny = ref_grid["ny"]; nx = ref_grid["nx"]
    dx = ref_grid["dx"]; dy = ref_grid["dy"]

    # ── box/domain time series (all 8 cases) ──────────────────────────────────
    print("\nComputing box + domain time series (all cases) ...")
    box_ts, domain_ts = compute_box_and_domain_time_series(case_results)

    # ── dCRE decompositions for all 4 pairs ───────────────────────────────────
    print("\nComputing dCRE decompositions (pixel-wise + mean-based, all pairs) ...")
    (box_decomp_pw, domain_decomp_pw,
     box_decomp_mb, domain_decomp_mb) = compute_dcre_decomp_time_series(
        case_results, file_by_case_time, ny, nx)

    # ── save CSVs ─────────────────────────────────────────────────────────────
    for pair_name in DECOMP_PAIRS:
        stub = pair_name.replace("-", "_")
        for label, df in [
            ("pixelwise", box_decomp_pw.get(pair_name, pd.DataFrame())),
            ("meanbased", box_decomp_mb.get(pair_name, pd.DataFrame())),
        ]:
            if not df.empty:
                p = os.path.join(OUTDIR, f"dcre_{label}_{stub}_by_box.csv")
                df.to_csv(p); print(f"  Saved {p}")
        for label, df in [
            ("pixelwise", domain_decomp_pw.get(pair_name, pd.DataFrame())),
            ("meanbased", domain_decomp_mb.get(pair_name, pd.DataFrame())),
        ]:
            if not df.empty:
                p = os.path.join(OUTDIR, f"dcre_{label}_{stub}_domain.csv")
                df.to_csv(p); print(f"  Saved {p}")

    # ── LWP colour range ──────────────────────────────────────────────────────
    global_vmin, global_vmax = compute_global_lwp_range(file_by_case_time)

    # ── per-timestep launched boxes (ML winds) ────────────────────────────────
    print(f"\nTracing per-timestep launched boxes ({LAUNCHED_TRACE_CASE} wind) ...")
    launched_states = trace_all_launched_boxes(
        all_files[LAUNCHED_TRACE_CASE], ny, nx, dx, dy)

    print("\nComputing LWP + dCRE for launched trajectories ...")
    traj_lwp, traj_dcre = compute_launched_traj_data(
        launched_states, all_files, file_by_case_time, ny, nx)

    # ── per-timestep Figures 1 & 2 (ML cases / pairs only) ───────────────────
    all_plot_times = sorted(set().union(
        *[set(file_by_case_time[c].keys()) for c in ML_CASES]))
    if MAKE_MAP_FIGURES:
        print(f"\nGenerating Figures 1 & 2 for {len(all_plot_times)} times ...")
        for t in all_plot_times:
            plot_one_time(t, case_results, file_by_case_time, box_ts, domain_ts,
                          box_decomp_pw, domain_decomp_pw,
                          box_decomp_mb, domain_decomp_mb,
                          global_vmin, global_vmax)
            plot_additional_figure(t, case_results, file_by_case_time, box_ts,
                                   domain_ts, global_vmin, global_vmax)
    else:
        print(f"\nSkipping per-timestep Figures 1 & 2 ({len(all_plot_times)} times, "
              f"~{2 * len(all_plot_times) * 2.8 / 1024:.1f} GB).  "
              f"Set WRF_MAP_FIGURES=1 to generate them.")

    # ── Figure 3: time-averaged dCRE barplot (all 4 pairs) ───────────────────
    print("\nGenerating Figure 3: time-averaged dCRE bar plot (all 4 pairs) ...")
    plot_dCRE_timeavg_barplot(box_decomp_pw)

    # ── Figure 4: LWP-stratified dCRE (ML pairs only) ────────────────────────
    print("\nGenerating Figure 4: LWP-stratified dCRE (ML pairs) ...")
    plot_lwp_stratified_dcre_barplot(traj_lwp, traj_dcre)

    # ── Figure 5: domain dCRE barplot (all 4 pairs) ──────────────────────────
    print("\nGenerating Figure 5: domain dCRE bar plot (all 4 pairs) ...")
    plot_domain_dcre_barplot(domain_decomp_pw)

    # ── Figure 6: met-binned dCRE (all 4 pairs, 8 met variables) ─────────────
    print("\nGenerating Figure 6: dCRE vs met regime (all 4 pairs) ...")
    for pair_name in DECOMP_PAIRS:
        print(f"  Building sample df for {pair_name} ...")
        sample_df = build_sample_dataframe(box_decomp_pw, box_ts, pair_name)
        print(f"    {len(sample_df)} box-timestep samples")
        plot_dCRE_met_binned_stacked(sample_df, pair_name)

    # ── Figure 7: pixel-wise vs mean-based comparison (all 4 pairs) ──────────
    print("\nGenerating Figure 7: pixel-wise vs mean-based (all 4 pairs) ...")
    for pair_name in DECOMP_PAIRS:
        plot_decomp_comparison_barplot(box_decomp_pw, domain_decomp_pw,
                                       box_decomp_mb, domain_decomp_mb,
                                       pair_name)

    # ── Figure 8: ML vs ARG domain comparison ────────────────────────────────
    print("\nGenerating Figure 8: ML vs ARG domain comparison ...")
    plot_ML_ARG_domain_comparison(domain_decomp_pw, domain_decomp_mb)

    # ── Figure 10 (v22): hourly + all-times ML vs ARG domain comparison ──────
    print("\nGenerating Figure 10: ML vs ARG domain comparison, hourly + all times ...")
    plot_ML_ARG_domain_comparison_hourly(domain_decomp_pw)

    # ── Figure 11 (v22): half-hourly stacked dCRE time series (per pair) ─────
    print(f"\nGenerating Figure 11: {STACK_BIN_MINUTES}-min stacked dCRE time "
          f"series (all 4 pairs) ...")
    for pair_name in DECOMP_PAIRS:
        plot_dCRE_halfhour_stacked_timeseries(domain_decomp_pw, pair_name)

    # ── Figure 9: dispersion β–Nc consistency check ──────────────────────────
    print("\nGenerating Figure 9: dispersion (β–Nc) consistency check ...")
    beta_stats, beta_samples = compute_beta_nc_consistency(file_by_case_time)
    df_beta_case, df_beta_pair = report_dispersion_consistency(beta_stats)
    plot_dispersion_consistency(beta_stats, beta_samples,
                                df_beta_case, df_beta_pair)

    print("\nAll done.")


if __name__ == "__main__":
    main()
