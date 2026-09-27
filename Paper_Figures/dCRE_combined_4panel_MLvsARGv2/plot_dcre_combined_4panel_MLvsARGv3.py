#!/usr/bin/env python3
"""
plot_dcre_combined_4panel_MLvsARGv3.py

v3: identical data and layout to v2, but every font is 1.35x larger and the
panel titles are reduced to the bare letters (a)-(d); the descriptive text
that followed the letter now belongs in the caption.

As plot_dcre_combined_4panel.py, but the left column shows BOTH schemes: in every
half-hour bin the ML-ACT stack and the ARG-ACT stack are drawn side by side
(solid vs '///' hatched), so the two can be compared bin by bin rather than only
in the time-mean bars on the right.  The right column is unchanged.

The left column is twice as wide as the right, because a time series needs the
horizontal room and a six-bar group does not.

    (a) top left     org background, half-hourly stacked dCRE, ML | ARG
    (b) top right    org background, ML-ACT vs ARG-ACT               [bar group]
    (c) bottom left  mid background, half-hourly stacked dCRE, ML | ARG
    (d) bottom right mid background, ML-ACT vs ARG-ACT               [bar group]

Panels are lettered in reading order (left to right, then down).

Everything is read from the CSVs written by
trace_east_boundary_subdomains_lwp_v22.py -- dcre_pixelwise_<pair>_domain.csv --
so this is a pure replot and never touches a wrfout file.

The window ends at 2017-07-15 13:00 (--tend), applied when the CSV is loaded so
the stacks and the time-mean bars cover the same times.  The x-axis is ticked on
the half-hour bin edges and padded slightly past the last one, so the closing
13:00 label is inside the frame.

v2: the ARG cases are now the WRF_dm_v2 re-runs.  That tree adds
re_cloud/re_ice/re_snow to the morr_two_moment package in Registry.EM_COMMON,
which WRF_dm omitted, so TAU_QC_TOT and RE_QC are populated instead of
identically zero.  R_tau is therefore a real ratio and the ARG sub-terms
(Nc / LWP / re / cov) are MEASUREMENTS, not the analytic 1/3 : 2/3 Twomey
fallback that every earlier version had to use.  The ML and ARG component
breakdowns can now be compared directly, not just their stack totals.

Whether the fallback applies is detected from the data (R_tau pinned at 1)
rather than assumed, so pointing --outdir at an old WRF_dm run still greys the
affected bars and restores the caveat automatically.

  python plot_dcre_combined_4panel_MLvsARGv2.py \
      [--outdir dcre_figures_v23_argv2_1300]
"""
import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# every font size below is scaled by this; 1.0 reproduces the
# original layout
FS = 1.5

# ── bar panels (right column) ────────────────────────────────────────────────
BAR_COLS = ["DCRE_total", "DCRE_CF", "DCRE_A_Nc",
            "DCRE_A_LWP", "DCRE_A_re", "DCRE_A_cov"]
BAR_LABELS = [r"$\Delta$CRE$_{\rm total}$",
              r"$\Delta$CRE$_{\rm CF}$",
              r"$\Delta$CRE$_{N_d}$",
              r"$\Delta$CRE$_{\rm LWP}$",
              r"$\Delta$CRE$_{r_e}$",
              r"$\Delta$CRE$_{\rm cov}$"]
ARG_FALLBACK = {"DCRE_A_Nc", "DCRE_A_LWP", "DCRE_A_re", "DCRE_A_cov"}


def arg_terms_are_fallback(df, tol=1e-3):
    """True if this record's R_tau carries no information.

    The apportionment is analytic exactly when tau does not respond, i.e. when
    R_tau is absent or pinned at 1 -- the WRF_dm case.  Detected rather than
    hard-coded so the same script serves both ARG trees.
    """
    if "R_tau" not in df:
        return True
    r = pd.to_numeric(df["R_tau"], errors="coerce").dropna()
    return len(r) == 0 or bool((r - 1.0).abs().max() < tol)

# ── stacked time-series panels (left column) ─────────────────────────────────
STACK_SPECS = [
    ("DCRE_CF",    "#1f77b4", r"$\Delta$CRE$_{\rm CF}$"),
    ("DCRE_A_Nc",  "#ff7f0e", r"$\Delta$CRE$_{N_d}$"),
    ("DCRE_A_LWP", "#2ca02c", r"$\Delta$CRE$_{\rm LWP}$"),
    ("DCRE_A_re",  "#d62728", r"$\Delta$CRE$_{r_e}$"),
    ("DCRE_A_cov", "#9467bd", r"$\Delta$CRE$_{\rm cov}$"),
]
STACK_BIN_MINUTES   = 30
MIN_TIMES_PER_STACK = 2

# row -> (background key, panel title, colour, ML stub, ARG stub, ARG pair label)
ROWS = [
    ("org", "org background", "#4393c3", "org_3aer_org", "ARG_org_3aer_ARG_org",
     "ARG_org_3aer-ARG_org"),
    ("mid", "mid background", "#d6604d", "mid_3aer_mid", "ARG_mid_3aer_ARG_mid",
     "ARG_mid_3aer-ARG_mid"),
]


def sem(s):
    n = s.count()
    return s.std() / np.sqrt(n) if n > 1 else np.nan


def load(outdir, stub, tend=None):
    f = os.path.join(outdir, f"dcre_pixelwise_{stub}_domain.csv")
    if not os.path.exists(f):
        raise SystemExit(f"missing {f}")
    df = pd.read_csv(f)
    if "time" in df.columns:
        df["time"] = pd.to_datetime(df["time"])
        df = df.set_index("time").sort_index()
        # Clipped here rather than in each panel, so the half-hourly stacks and
        # the time-mean bars/totals are computed over exactly the same window.
        if tend is not None:
            df = df.loc[df.index <= tend]
            if df.empty:
                raise SystemExit(f"{f}: no rows at or before {tend}")
    return df


def _binned(df, starts):
    """Per-bin component means, total mean, total SE and count for one pair."""
    grp = df.groupby(df.index.floor(f"{STACK_BIN_MINUTES}min"))
    n_per = grp.size()
    comps = {}
    for comp, _, _ in STACK_SPECS:
        comps[comp] = (np.array([grp[comp].mean().get(s, np.nan) for s in starts],
                                float)
                       if comp in df else np.full(len(starts), np.nan))
    tot = np.array([grp["DCRE_total"].mean().get(s, np.nan) for s in starts], float)
    err = np.array([(grp["DCRE_total"].std().get(s, np.nan) /
                     np.sqrt(max(int(n_per.get(s, 0)), 1))) for s in starts], float)
    n = np.array([int(n_per.get(s, 0)) for s in starts])
    return comps, tot, err, n


def draw_stacked_timeseries_pair(ax, df_ml, df_arg, tag, title):
    """Half-hourly stacked dCRE for BOTH schemes, ML and ARG side by side.

    Bins come from the ML record; the two runs share the output times, so the
    ARG stack for a bin covers the same instants.  Solid = ML, '///' = ARG, to
    match the fill convention of the bar panels on the right.
    """
    ml_bins  = set(df_ml.index.floor(f"{STACK_BIN_MINUTES}min"))
    arg_bins = set(df_arg.index.floor(f"{STACK_BIN_MINUTES}min"))
    starts = pd.DatetimeIndex(sorted(ml_bins & arg_bins))
    counts = df_ml.groupby(df_ml.index.floor(f"{STACK_BIN_MINUTES}min")).size()
    starts = pd.DatetimeIndex([s for s in starts
                               if counts.get(s, 0) >= MIN_TIMES_PER_STACK])

    centre = starts + pd.Timedelta(minutes=STACK_BIN_MINUTES) / 2
    # two bars per bin, each a bit under half the bin so they do not touch
    width = (STACK_BIN_MINUTES / (24.0 * 60.0)) * 0.42
    off = pd.Timedelta(minutes=STACK_BIN_MINUTES) * 0.22

    for df, hatch, alpha, shift, is_ml in ((df_ml,  "",    0.88, -off, True),
                                           (df_arg, "///", 0.60,  off, False)):
        comps, tot, err, n = _binned(df, starts)
        xs = centre + shift
        pos_bot = np.zeros(len(starts))
        neg_bot = np.zeros(len(starts))
        for comp, color, clabel in STACK_SPECS:
            vals = comps[comp]
            pos_v = np.where(np.isfinite(vals) & (vals > 0), vals, 0.0)
            neg_v = np.where(np.isfinite(vals) & (vals < 0), vals, 0.0)
            ax.bar(xs, pos_v, width=width, bottom=pos_bot, color=color,
                   alpha=alpha, hatch=hatch, edgecolor="k", linewidth=0.35,
                   label=clabel if is_ml else None, zorder=3)
            ax.bar(xs, neg_v, width=width, bottom=neg_bot, color=color,
                   alpha=alpha, hatch=hatch, edgecolor="k", linewidth=0.35,
                   zorder=3)
            pos_bot += pos_v
            neg_bot += neg_v
        ax.errorbar(xs, tot, yerr=err,
                    fmt="D" if is_ml else "o",
                    color="black",
                    markerfacecolor="black" if is_ml else "white",
                    ms=7, capsize=3.5, linewidth=1.2, zorder=6,
                    label=(r"$\Delta$CRE$_{\rm total}$ ±SE (ML)" if is_ml
                           else r"$\Delta$CRE$_{\rm total}$ ±SE (ARG)"))
        if is_ml:
            n_ml = n

    # the larger legend needs room of its own below the deepest stack, so the
    # lower limit is pushed down by a fixed share of the data range
    lo, hi = ax.get_ylim()
    lo = lo - 0.58 * (hi - lo)
    ax.set_ylim(lo, hi)
    y_low = lo

    ax.axhline(0, color="k", linewidth=0.9, zorder=4)
    ax.set_ylabel(r"$\Delta$CRE  [W m$^{-2}$]", fontsize=28.4)
    ax.set_xlabel("time (UTC)", fontsize=28.4)
    # Tick the bin EDGES, not the bin centres.  The bars are centred in their
    # bins, so autoscaling on them left the axis ending at the last sample
    # (12:58) with the closing edge of the last bin (13:00) neither drawn nor
    # labelled; edge ticks also read correctly for binned data.
    edges = starts.append(pd.DatetimeIndex(
        [starts[-1] + pd.Timedelta(minutes=STACK_BIN_MINUTES)]))
    ax.set_xticks(edges)
    ax.set_xticklabels([f"{e:%H:%M}" for e in edges], rotation=35,
                       ha="right", fontsize=22.3)
    # a little air past the last edge, so the closing label (13:00) sits
    # inside the frame rather than flush against the spine
    bin_w = pd.Timedelta(minutes=STACK_BIN_MINUTES)
    ax.set_xlim(edges[0] - 0.15 * bin_w, edges[-1] + 0.35 * bin_w)
    ax.tick_params(axis="y", labelsize=24)
    ax.grid(True, alpha=0.3, axis="y", zorder=0)
    ax.set_title(f"({tag})", fontsize=32.4, loc="left")

    handles = [Patch(facecolor=c, alpha=0.88, label=lbl) for _, c, lbl in STACK_SPECS]
    handles += [
        Patch(facecolor="0.75", edgecolor="k", alpha=0.88, label="ML-ACT"),
        Patch(facecolor="0.75", edgecolor="k", alpha=0.60, hatch="///",
              label="ARG-ACT"),
        Line2D([0], [0], marker="D", color="black", ms=7, linewidth=1.2,
               label=r"$\Delta$CRE$_{\rm total}$ (ML)"),
        Line2D([0], [0], marker="o", color="black", markerfacecolor="white",
               ms=7, linewidth=1.2, label=r"$\Delta$CRE$_{\rm total}$ (ARG)"),
    ]
    ax.legend(handles=handles, fontsize=23, frameon=True, ncol=3, loc="lower left")


def draw_ml_arg_bars(ax, df_ml, df_arg, tag, title, color):
    """Grouped ML vs ARG bars of the dCRE components on `ax`."""
    x  = np.arange(len(BAR_COLS))
    bw = 0.36
    fb = arg_terms_are_fallback(df_arg)
    for df, lab, hatch, alpha, off in ((df_ml,  "ML-ACT",  "",    0.90, -0.5),
                                       (df_arg, "ARG-ACT", "///", 0.55,  0.5)):
        means, errs, cols = [], [], []
        for c in BAR_COLS:
            v = (pd.to_numeric(df[c], errors="coerce")
                 if c in df else pd.Series(dtype=float))
            means.append(v.mean() if len(v) else np.nan)
            errs.append(sem(v) if len(v) else np.nan)
            cols.append("#b0b0b0"
                        if (fb and lab.startswith("ARG") and c in ARG_FALLBACK)
                        else color)
        ax.bar(x + off * bw, means, bw, yerr=errs, capsize=4,
               error_kw={"linewidth": 1.1}, color=cols, hatch=hatch,
               alpha=alpha, edgecolor="k", linewidth=0.6, label=lab)

    ax.axhline(0, color="k", lw=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(BAR_LABELS, rotation=25, ha="right", fontsize=24.3)
    ax.set_ylabel(r"$\Delta$CRE  [W m$^{-2}$]", fontsize=28.4)
    ax.tick_params(axis="y", labelsize=24)
    ax.set_title(f"({tag})", fontsize=32.4, loc="left")
    ax.grid(True, alpha=0.3, axis="y")
    ax.axvspan(1.5, len(BAR_COLS) - 0.5, color="0.5", alpha=0.06, zorder=0)

    # Legend drawn in THIS panel's colour: one shared legend would have to pick
    # a single hue and would then misrepresent the other row.
    _lh = [Patch(facecolor=color, edgecolor="k", alpha=0.90, label="ML-ACT"),
           Patch(facecolor=color, edgecolor="k", alpha=0.55, hatch="///",
                 label="ARG-ACT")]
    if fb:
        _lh.append(Patch(facecolor="#b0b0b0", edgecolor="k", alpha=0.55,
                         hatch="///", label="ARG-ACT (analytic fallback)"))
    ax.legend(handles=_lh, fontsize=24, frameon=True, loc="upper left")

    # v3: no ML/ARG time-mean summary box; those numbers are quoted in the text.
    # A little headroom above the bars keeps the legend off the ARG total bar.
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + 0.36 * (hi - lo))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="dcre_figures_v23_argv2_1300")
    ap.add_argument("--stem", default="dCRE_combined_4panel_MLvsARGv3")
    ap.add_argument("--tend", default="2017-07-15 13:00",
                    help="end of the window, for BOTH the time series and the "
                         "time-mean bars; empty string keeps the whole CSV")
    a = ap.parse_args()
    tend = pd.Timestamp(a.tend) if a.tend else None

    fig = plt.figure(figsize=(24.0, 13.0))
    # left column twice as wide as the right
    gs = GridSpec(2, 2, figure=fig, width_ratios=[2, 1],
                  hspace=0.42, wspace=0.20,
                  top=0.955, bottom=0.085, left=0.065, right=0.985)

    # reading order: (a) top-left, (b) top-right, (c) bottom-left, (d) bottom-right
    tags = [["a", "b"], ["c", "d"]]

    for irow, (bg, title, color, ml_stub, arg_stub, arg_label) in enumerate(ROWS):
        df_ml  = load(a.outdir, ml_stub, tend)
        df_arg = load(a.outdir, arg_stub, tend)

        # left: both schemes' stacked time series, side by side per bin
        draw_stacked_timeseries_pair(fig.add_subplot(gs[irow, 0]),
                                     df_ml, df_arg, tags[irow][0], title)

        # right: ML vs ARG bars for the same background
        draw_ml_arg_bars(fig.add_subplot(gs[irow, 1]), df_ml, df_arg,
                         tags[irow][1], title, color)

    for ext in ("png", "pdf", "eps"):
        out = os.path.join(a.outdir, f"{a.stem}.{ext}")
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"  wrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
