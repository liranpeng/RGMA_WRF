#!/usr/bin/env python3
"""
plot_dcre_combined_4panel.py

One 2x2 figure combining the two ARG stacked time series with the two panels of
the all-times ML-vs-ARG domain comparison.  The left column is twice as wide as
the right, because a time series needs the horizontal room and a six-bar group
does not.

    (a) top left     ARG_org_3aer-ARG_org, half-hourly stacked dCRE  [time series]
    (b) top right    org background, ML-ACT vs ARG-ACT               [bar group]
    (c) bottom left  ARG_mid_3aer-ARG_mid, half-hourly stacked dCRE  [time series]
    (d) bottom right mid background, ML-ACT vs ARG-ACT               [bar group]

Panels are lettered in reading order (left to right, then down).

Everything is read from the CSVs written by
trace_east_boundary_subdomains_lwp_v22.py -- dcre_pixelwise_<pair>_domain.csv --
so this is a pure replot and never touches a wrfout file.

The window ends at 2017-07-15 13:00 (--tend).  The clip is applied when the CSV
is loaded, so the half-hourly stacks on the left and the time-mean bars and
totals on the right are all computed over the same times.  The generator caps
its shared-time set at the same instant, so on freshly written CSVs this is a
no-op; it still matters when replotting older CSVs that run past 13:00.  The
x-axis is ticked on the half-hour bin edges and padded slightly past the last
one, so the closing 13:00 label is inside the frame.

NOTE the ARG sub-terms (Nc / LWP / re / cov) in the right-hand panels are NOT
measurements.  TAU_QC_TOT and RE_QC are identically zero throughout the WRF_dm
build, so R_tau degenerates to 1 and the apportionment falls back to the
analytic 1/3 : 2/3 Twomey split.  They are drawn grey and annotated.  Only the
total, CF and A bars are informative for ARG.  The left-hand time series are ARG
pairs too, so the same caveat applies to their Nc / LWP / re stacking: the split
between those three components is analytic, while their sum and the total are
measured.

  python plot_dcre_combined_4panel.py [--outdir dcre_figures_v22]
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

# ── bar panels (right column) ────────────────────────────────────────────────
BAR_COLS = ["DCRE_total", "DCRE_CF", "DCRE_A_Nc",
            "DCRE_A_LWP", "DCRE_A_re", "DCRE_A_cov"]
BAR_LABELS = [r"$\Delta$CRE$_{\rm total}$",
              r"$\Delta$CRE$_{\rm CF}$",
              r"$\Delta$CRE$_{A,Nc}$",
              r"$\Delta$CRE$_{A,LWP}$",
              r"$\Delta$CRE$_{A,re}$",
              r"$\Delta$CRE$_{A,cov}$"]
ARG_FALLBACK = {"DCRE_A_Nc", "DCRE_A_LWP", "DCRE_A_re", "DCRE_A_cov"}

# ── stacked time-series panels (left column) ─────────────────────────────────
STACK_SPECS = [
    ("DCRE_CF",    "#1f77b4", r"$\Delta$CRE$_{\rm CF}$"),
    ("DCRE_A_Nc",  "#ff7f0e", r"$\Delta$CRE$_{A,Nc}$"),
    ("DCRE_A_LWP", "#2ca02c", r"$\Delta$CRE$_{A,LWP}$"),
    ("DCRE_A_re",  "#d62728", r"$\Delta$CRE$_{A,re}$"),
    ("DCRE_A_cov", "#9467bd", r"$\Delta$CRE$_{A,cov}$"),
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


def draw_stacked_timeseries(ax, df, tag, pair_label):
    """Half-hourly stacked bars of the domain dCRE decomposition on `ax`."""
    bin_start = df.index.floor(f"{STACK_BIN_MINUTES}min")
    grp = df.groupby(bin_start)
    n_per = grp.size()
    keep = n_per[n_per >= MIN_TIMES_PER_STACK].index
    starts = pd.DatetimeIndex(sorted(keep))
    centres = starts + pd.Timedelta(minutes=STACK_BIN_MINUTES) / 2
    width = (STACK_BIN_MINUTES / (24.0 * 60.0)) * 0.88

    pos_bot = np.zeros(len(starts))
    neg_bot = np.zeros(len(starts))
    for comp, color, clabel in STACK_SPECS:
        vals = (np.array([grp[comp].mean().get(s, np.nan) for s in starts], float)
                if comp in df else np.full(len(starts), np.nan))
        pos_v = np.where(np.isfinite(vals) & (vals > 0), vals, 0.0)
        neg_v = np.where(np.isfinite(vals) & (vals < 0), vals, 0.0)
        ax.bar(centres, pos_v, width=width, bottom=pos_bot,
               color=color, alpha=0.85, label=clabel, zorder=3)
        ax.bar(centres, neg_v, width=width, bottom=neg_bot,
               color=color, alpha=0.85, zorder=3)
        pos_bot += pos_v
        neg_bot += neg_v

    tot = np.array([grp["DCRE_total"].mean().get(s, np.nan) for s in starts], float)
    err = np.array([(grp["DCRE_total"].std().get(s, np.nan) /
                     np.sqrt(max(int(n_per.get(s, 0)), 1))) for s in starts], float)
    ax.errorbar(centres, tot, yerr=err, fmt="D", color="black", ms=5,
                capsize=4, linewidth=1.2, zorder=6,
                label=r"$\Delta$CRE$_{\rm total}$ ±SE")

    y_low = ax.get_ylim()[0]
    for c, s in zip(centres, starts):
        ax.text(c, y_low, f"n={int(n_per.get(s, 0))}", ha="center", va="bottom",
                fontsize=7, color="dimgray", zorder=7)

    ax.axhline(0, color="k", linewidth=0.9, zorder=4)
    ax.set_ylabel(r"$\Delta$CRE  [W m$^{-2}$]", fontsize=14)
    ax.set_xlabel("time (UTC)", fontsize=14)
    # Tick the bin EDGES, not the bin centres.  The bars are centred in their
    # bins, so autoscaling on them left the axis ending at the last sample
    # (12:58) with the closing edge of the last bin (13:00) neither drawn nor
    # labelled; edge ticks also read correctly for binned data.
    edges = starts.append(pd.DatetimeIndex(
        [starts[-1] + pd.Timedelta(minutes=STACK_BIN_MINUTES)]))
    ax.set_xticks(edges)
    ax.set_xticklabels([f"{e:%H:%M}" for e in edges], rotation=35,
                       ha="right", fontsize=11)
    # a little air past the last edge, so the closing label (13:00) sits
    # inside the frame rather than flush against the spine
    bin_w = pd.Timedelta(minutes=STACK_BIN_MINUTES)
    ax.set_xlim(edges[0] - 0.15 * bin_w, edges[-1] + 0.35 * bin_w)
    ax.tick_params(axis="y", labelsize=12)
    ax.grid(True, alpha=0.3, axis="y", zorder=0)
    ax.set_title(f"({tag}) {pair_label}: pixel-wise $\\Delta$CRE per "
                 f"{STACK_BIN_MINUTES} min", fontsize=16, loc="left")
    ax.legend(fontsize=9.5, frameon=True, ncol=3, loc="lower left")


def draw_ml_arg_bars(ax, df_ml, df_arg, tag, title, color):
    """Grouped ML vs ARG bars of the dCRE components on `ax`."""
    x  = np.arange(len(BAR_COLS))
    bw = 0.36
    for df, lab, hatch, alpha, off in ((df_ml,  "ML-ACT",  "",    0.90, -0.5),
                                       (df_arg, "ARG-ACT", "///", 0.55,  0.5)):
        means, errs, cols = [], [], []
        for c in BAR_COLS:
            v = (pd.to_numeric(df[c], errors="coerce")
                 if c in df else pd.Series(dtype=float))
            means.append(v.mean() if len(v) else np.nan)
            errs.append(sem(v) if len(v) else np.nan)
            cols.append("#b0b0b0"
                        if (lab.startswith("ARG") and c in ARG_FALLBACK) else color)
        ax.bar(x + off * bw, means, bw, yerr=errs, capsize=4,
               error_kw={"linewidth": 1.1}, color=cols, hatch=hatch,
               alpha=alpha, edgecolor="k", linewidth=0.6, label=lab)

    ax.axhline(0, color="k", lw=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(BAR_LABELS, rotation=25, ha="right", fontsize=12)
    ax.set_ylabel(r"$\Delta$CRE  [W m$^{-2}$]", fontsize=14)
    ax.tick_params(axis="y", labelsize=12)
    ax.set_title(f"({tag}) {title}: ML-ACT vs ARG-ACT", fontsize=16, loc="left")
    ax.grid(True, alpha=0.3, axis="y")
    ax.axvspan(1.5, len(BAR_COLS) - 0.5, color="0.5", alpha=0.06, zorder=0)

    # Legend drawn in THIS panel's colour: one shared legend would have to pick
    # a single hue and would then misrepresent the other row.
    ax.legend(handles=[
        Patch(facecolor=color, edgecolor="k", alpha=0.90, label="ML-ACT"),
        Patch(facecolor=color, edgecolor="k", alpha=0.55, hatch="///",
              label="ARG-ACT"),
        Patch(facecolor="#b0b0b0", edgecolor="k", alpha=0.55, hatch="///",
              label="ARG-ACT (analytic fallback)"),
    ], fontsize=9.5, frameon=True, loc="upper left")

    tot_ml  = pd.to_numeric(df_ml["DCRE_total"], errors="coerce").mean()
    tot_arg = pd.to_numeric(df_arg["DCRE_total"], errors="coerce").mean()
    ratio = tot_arg / tot_ml if (np.isfinite(tot_ml) and tot_ml != 0) else np.nan
    ax.text(0.015, 0.03,
            f"ML {tot_ml:+.2f}   ARG {tot_arg:+.2f} W m$^{{-2}}$\n"
            f"(ARG/ML = {ratio:.2f}),  n={len(df_ml)} times",
            transform=ax.transAxes, fontsize=11, va="bottom",
            bbox=dict(fc="white", alpha=0.85, ec="0.7"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="dcre_figures_v22")
    ap.add_argument("--stem", default="dCRE_combined_4panel")
    ap.add_argument("--tend", default="2017-07-15 13:00",
                    help="end of the window, for BOTH the time series and the "
                         "time-mean bars; empty string keeps the whole CSV")
    a = ap.parse_args()
    tend = pd.Timestamp(a.tend) if a.tend else None

    fig = plt.figure(figsize=(24.0, 13.0))
    # left column twice as wide as the right
    gs = GridSpec(2, 2, figure=fig, width_ratios=[2, 1],
                  hspace=0.42, wspace=0.20,
                  top=0.90, bottom=0.11, left=0.055, right=0.985)

    # reading order: (a) top-left, (b) top-right, (c) bottom-left, (d) bottom-right
    tags = [["a", "b"], ["c", "d"]]

    for irow, (bg, title, color, ml_stub, arg_stub, arg_label) in enumerate(ROWS):
        df_ml  = load(a.outdir, ml_stub, tend)
        df_arg = load(a.outdir, arg_stub, tend)

        # left: the ARG pair's stacked time series
        draw_stacked_timeseries(fig.add_subplot(gs[irow, 0]), df_arg,
                                tags[irow][0], arg_label)

        # right: ML vs ARG bars for the same background
        draw_ml_arg_bars(fig.add_subplot(gs[irow, 1]), df_ml, df_arg,
                         tags[irow][1], title, color)

    t0 = load(a.outdir, ROWS[0][4], tend).index
    fig.suptitle(
        "Pixel-wise $\\Delta$CRE decomposition for the $3\\times$ sulfate "
        "perturbation, aggregated over the full domain\n"
        "left: ARG-ACT half-hourly evolution   |   "
        "right: time-mean ML-ACT vs ARG-ACT components   |   "
        f"all shared times, {t0.min():%Y-%m-%d %H:%M}–{t0.max():%H:%M} UTC "
        "(intersection across all 8 cases)",
        fontsize=19, y=0.985)
    # below the axes, not on top of the bottom-left x-label; bbox_inches="tight"
    # at save time expands the canvas to include it
    fig.text(0.5, 0.03,
             "Colour in the right-hand panels denotes the background "
             "(blue = org, red = mid); fill denotes the activation treatment.  "
             "Grey bars, and the Nc / LWP / re split of the left-hand stacks, "
             "are the analytic fallback, not a measurement.",
             ha="center", fontsize=12)

    for ext in ("png", "pdf", "eps"):
        out = os.path.join(a.outdir, f"{a.stem}.{ext}")
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"  wrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
