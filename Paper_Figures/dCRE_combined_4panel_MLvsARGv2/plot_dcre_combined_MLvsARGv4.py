#!/usr/bin/env python3
"""
plot_dcre_combined_MLvsARGv4.py

v4 = v3 (large fonts, letter-only titles, no summary box, N_d notation) plus:
(--backgrounds org,mid --stem dCRE_combined_4panel_MLvsARGv5 gives v5: the same
figure without mid2 anywhere -- two background rows and composites pooled over
org and mid only; adding --composite-backgrounds org,mid,mid2 gives v6: org and
mid rows, but the composites pooled over all three backgrounds)

  * a third background row, mid2 (the moisture midpoint between org and mid):
        (e) mid2, half-hourly stacked dCRE, ML | ARG      (f) mid2 bars
  * a bottom row of STATE COMPOSITES, built exactly like (a)/(c)/(e) but with
    the x-axis binned in cloud state instead of time:
        (g) vs N_d      (h) vs LWP      (i) vs cloud fraction
    The sampling unit is one 30-min window of one pair (the same windows that
    make one bar group in (a)/(c)/(e)).  The 30-min means of all three
    backgrounds are pooled -- 8 windows x 3 backgrounds = 24 samples per
    scheme -- and sorted into fixed bins of the pair-mean state.  In each bin
    the ML and ARG stacks are drawn side by side (solid | hatched) with the
    total +-SE across the samples in the bin and the sample count n=ML|ARG.

    State per pair and time = mean of the 1aer and 3aer runs of
      N_d  : QNDROP, in-cloud, air-mass weighted below 2 km      [cm-3]
      LWP  : LWP_TOT, mean over cloudy columns                    [g m-2]
      CF   : CLDCOL, cloudy-column fraction of the domain         [-]
    read from the consequences-figure cache (state_cache/*.npz, produced by
    plot_arg_vs_ml_consequences_rwp_v11.py) and joined to the dCRE CSVs on the
    output time.  No wrfout is read.

Inputs (all in --outdir, default dcre_csv):
    dcre_pixelwise_{org,mid,mid2}_3aer_{org,mid,mid2}_domain.csv          ML
    dcre_pixelwise_ARG_{org,mid,mid2}_3aer_ARG_{org,mid,mid2}_domain.csv  ARG
  and --state-cache (default state_cache/arg_vs_ml_consequences_rwp_v11_series_CCN2.npz)

Outputs: <stem>.{png,pdf,eps} plus three tables next to the figure:
    <stem>_state_merged_2min.csv    every output time, dCRE terms + state, all pairs
    <stem>_state_samples_30min.csv  the pooled 30-min samples behind (g)-(i)
    <stem>_state_composites.csv     per-bin means / SE / n drawn in (g)-(i)
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
              r"$\Delta$CRE$_{N_d}$",
              r"$\Delta$CRE$_{\rm LWP}$",
              r"$\Delta$CRE$_{r_e}$",
              r"$\Delta$CRE$_{\rm cov}$"]
ARG_FALLBACK = {"DCRE_A_Nc", "DCRE_A_LWP", "DCRE_A_re", "DCRE_A_cov"}


def arg_terms_are_fallback(df, tol=1e-3):
    if "R_tau" not in df:
        return True
    r = pd.to_numeric(df["R_tau"], errors="coerce").dropna()
    return len(r) == 0 or bool((r - 1.0).abs().max() < tol)

# ── stacked panels ───────────────────────────────────────────────────────────
STACK_SPECS = [
    ("DCRE_CF",    "#1f77b4", r"$\Delta$CRE$_{\rm CF}$"),
    ("DCRE_A_Nc",  "#ff7f0e", r"$\Delta$CRE$_{N_d}$"),
    ("DCRE_A_LWP", "#2ca02c", r"$\Delta$CRE$_{\rm LWP}$"),
    ("DCRE_A_re",  "#d62728", r"$\Delta$CRE$_{r_e}$"),
    ("DCRE_A_cov", "#9467bd", r"$\Delta$CRE$_{\rm cov}$"),
]
STACK_BIN_MINUTES   = 30
MIN_TIMES_PER_STACK = 2

# row -> (background key, colour, ML stub, ARG stub, ML cases, ARG cases)
ROWS = [
    ("org",  "#4393c3", "org_3aer_org",   "ARG_org_3aer_ARG_org",
     ("org", "org_3aer"),   ("ARG_org", "ARG_org_3aer")),
    ("mid",  "#d6604d", "mid_3aer_mid",   "ARG_mid_3aer_ARG_mid",
     ("mid", "mid_3aer"),   ("ARG_mid", "ARG_mid_3aer")),
    ("mid2", "#e08214", "mid2_3aer_mid2", "ARG_mid2_3aer_ARG_mid2",
     ("mid2", "mid2_3aer"), ("ARG_mid2", "ARG_mid2_3aer")),
]

# state composites: (column in the state table, cache series, panel x-label,
#                    fixed bin edges, tick-label format)
STATE_SPECS = [
    ("Nd",  "QNDROP",  r"$N_d$ bin  [cm$^{-3}$]",
     [0, 15, 30, 45, 60, np.inf], "{:.0f}"),
    ("LWP", "LWP_TOT", r"LWP bin  [g m$^{-2}$]",
     [0, 75, 150, np.inf], "{:.0f}"),
    ("CF",  "CLDCOL",  r"cloud-fraction bin  [-]",
     [0, 0.5, 0.8, 0.95, 1.0000001], "{:.2f}"),
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
        if tend is not None:
            df = df.loc[df.index <= tend]
            if df.empty:
                raise SystemExit(f"{f}: no rows at or before {tend}")
    return df


def load_state(npz_path):
    """{case: DataFrame[time -> Nd, LWP, CF]} from the consequences cache."""
    z = np.load(npz_path, allow_pickle=True)
    data = z["data"].item()
    out = {}
    for case, (times, series) in data.items():
        out[case] = pd.DataFrame(
            {col: np.asarray(series[key], float) for col, key, *_ in STATE_SPECS},
            index=pd.DatetimeIndex(times)).sort_index()
    return out


def pair_state(state, cases):
    """Pair-mean state: average of the control and perturbed runs, on the
    times both have."""
    a, b = (state[c] for c in cases)
    both = a.index.intersection(b.index)
    return (a.loc[both] + b.loc[both]) / 2.0


def _binned(df, starts):
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


def _stack_pair(ax, xs, comps, width, hatch, alpha, is_ml):
    """Draw one scheme's stacked components at positions xs; returns nothing."""
    pos_bot = np.zeros(len(xs))
    neg_bot = np.zeros(len(xs))
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


def _stack_handles():
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
    return handles


def _stack_legend(ax, fontsize=23, ncol=3):
    ax.legend(handles=_stack_handles(), fontsize=fontsize, frameon=True,
              ncol=ncol, loc="lower left")


def draw_stacked_timeseries_pair(ax, df_ml, df_arg, tag, legend=True):
    ml_bins  = set(df_ml.index.floor(f"{STACK_BIN_MINUTES}min"))
    arg_bins = set(df_arg.index.floor(f"{STACK_BIN_MINUTES}min"))
    starts = pd.DatetimeIndex(sorted(ml_bins & arg_bins))
    counts = df_ml.groupby(df_ml.index.floor(f"{STACK_BIN_MINUTES}min")).size()
    starts = pd.DatetimeIndex([s for s in starts
                               if counts.get(s, 0) >= MIN_TIMES_PER_STACK])

    centre = starts + pd.Timedelta(minutes=STACK_BIN_MINUTES) / 2
    width = (STACK_BIN_MINUTES / (24.0 * 60.0)) * 0.42
    off = pd.Timedelta(minutes=STACK_BIN_MINUTES) * 0.22

    for df, hatch, alpha, shift, is_ml in ((df_ml,  "",    0.88, -off, True),
                                           (df_arg, "///", 0.60,  off, False)):
        comps, tot, err, n = _binned(df, starts)
        xs = centre + shift
        _stack_pair(ax, xs, comps, width, hatch, alpha, is_ml)
        ax.errorbar(xs, tot, yerr=err,
                    fmt="D" if is_ml else "o",
                    color="black",
                    markerfacecolor="black" if is_ml else "white",
                    ms=7, capsize=3.5, linewidth=1.2, zorder=6)
        if is_ml:
            n_ml = n

    lo, hi = ax.get_ylim()
    # room for the legend when it is drawn (no per-bin n labels since v6)
    lo = lo - (0.40 if legend else 0.03) * (hi - lo)
    ax.set_ylim(lo, hi)

    ax.axhline(0, color="k", linewidth=0.9, zorder=4)
    ax.set_ylabel(r"$\Delta$CRE  [W m$^{-2}$]", fontsize=28.4)
    ax.set_xlabel("time (UTC)", fontsize=28.4)
    edges = starts.append(pd.DatetimeIndex(
        [starts[-1] + pd.Timedelta(minutes=STACK_BIN_MINUTES)]))
    ax.set_xticks(edges)
    ax.set_xticklabels([f"{e:%H:%M}" for e in edges], rotation=35,
                       ha="right", fontsize=22.3)
    bin_w = pd.Timedelta(minutes=STACK_BIN_MINUTES)
    ax.set_xlim(edges[0] - 0.15 * bin_w, edges[-1] + 0.35 * bin_w)
    ax.tick_params(axis="y", labelsize=24)
    ax.grid(True, alpha=0.3, axis="y", zorder=0)
    ax.set_title(f"({tag})", fontsize=32.4, loc="left")
    if legend:
        _stack_legend(ax)


def draw_ml_arg_bars(ax, df_ml, df_arg, tag, color):
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

    _lh = [Patch(facecolor=color, edgecolor="k", alpha=0.90, label="ML-ACT"),
           Patch(facecolor=color, edgecolor="k", alpha=0.55, hatch="///",
                 label="ARG-ACT")]
    if fb:
        _lh.append(Patch(facecolor="#b0b0b0", edgecolor="k", alpha=0.55,
                         hatch="///", label="ARG-ACT (analytic fallback)"))
    ax.legend(handles=_lh, fontsize=24, frameon=True, loc="upper left")

    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + 0.36 * (hi - lo))


# ── state composites ─────────────────────────────────────────────────────────

def window_samples(df, tag_scheme, tag_bg):
    """30-min means of one pair: dCRE terms and pair-mean state; one row per
    window with at least MIN_TIMES_PER_STACK output times."""
    key = df.index.floor(f"{STACK_BIN_MINUTES}min")
    grp = df.groupby(key)
    n = grp.size()
    cols = ["DCRE_total"] + [c for c, _, _ in STACK_SPECS] + [c for c, *_ in STATE_SPECS]
    g = grp[cols].mean()
    g["n_times"] = n
    g = g[g["n_times"] >= MIN_TIMES_PER_STACK]
    g.insert(0, "background", tag_bg)
    g.insert(0, "scheme", tag_scheme)
    g.index.name = "window_start"
    return g


def _bin_labels(edges, fmt):
    labs = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        if not np.isfinite(hi) or hi >= 1.0000001 and edges[-1] == 1.0000001 and hi == edges[-1]:
            labs.append((r"$\geq$" + fmt.format(lo)) if not np.isfinite(hi)
                        else fmt.format(lo) + "–" + fmt.format(min(hi, 1.0)))
        elif lo == edges[0] and lo == 0:
            labs.append("<" + fmt.format(hi))
        else:
            labs.append(fmt.format(lo) + "–" + fmt.format(hi))
    return labs


def draw_state_composite(ax, samples, xcol, edges, fmt, tag, xlabel,
                         legend=False):
    """Stacked dCRE per state bin, ML | ARG, pooled over all backgrounds.
    Bins that hold no sample of either scheme are left out of the axis."""
    edges = np.asarray(edges, float)
    nb = len(edges) - 1
    width, off = 0.38, 0.21
    rows, stats = [], {}
    for scheme in ("ML", "ARG"):
        s = samples[samples["scheme"] == scheme]
        ib = np.digitize(s[xcol].values, edges) - 1
        comps = {c: np.full(nb, np.nan) for c, _, _ in STACK_SPECS}
        tot = np.full(nb, np.nan); err = np.full(nb, np.nan); n = np.zeros(nb, int)
        for b in range(nb):
            sel = s[ib == b]
            n[b] = len(sel)
            if n[b] == 0:
                continue
            for c, _, _ in STACK_SPECS:
                comps[c][b] = sel[c].mean()
            tot[b] = sel["DCRE_total"].mean()
            err[b] = sem(sel["DCRE_total"])
            rows.append(dict(x=xcol, bin=b, bin_lo=edges[b], bin_hi=edges[b + 1],
                             scheme=scheme, n=n[b],
                             backgrounds="/".join(sorted(set(sel["background"]))),
                             x_mean=sel[xcol].mean(),
                             DCRE_total=tot[b], DCRE_total_se=err[b],
                             **{c: comps[c][b] for c, _, _ in STACK_SPECS}))
        stats[scheme] = (comps, tot, err, n)

    keep = np.where((stats["ML"][3] + stats["ARG"][3]) > 0)[0]
    pos = np.arange(len(keep))
    n_by_scheme = {}
    for scheme, hatch, alpha, shift, is_ml in (("ML",  "",    0.88, -off, True),
                                               ("ARG", "///", 0.60,  off, False)):
        comps, tot, err, n = stats[scheme]
        comps = {c: v[keep] for c, v in comps.items()}
        tot, err, n = tot[keep], err[keep], n[keep]
        n_by_scheme[scheme] = n
        xs = pos + shift
        _stack_pair(ax, xs, comps, width, hatch, alpha, is_ml)
        ok = np.isfinite(tot)
        ax.errorbar(xs[ok], tot[ok], yerr=np.where(np.isfinite(err[ok]), err[ok], 0),
                    fmt="D" if is_ml else "o", color="black",
                    markerfacecolor="black" if is_ml else "white",
                    ms=7, capsize=3.5, linewidth=1.2, zorder=6)

    ax.axhline(0, color="k", linewidth=0.9, zorder=4)
    ax.set_xticks(pos)
    labels = _bin_labels(edges, fmt)
    ax.set_xticklabels([labels[k] for k in keep], fontsize=22.3, rotation=30,
                       ha="right")
    ax.set_xlim(-0.6, len(keep) - 0.4)
    ax.set_xlabel(xlabel, fontsize=26)
    ax.set_ylabel(r"$\Delta$CRE  [W m$^{-2}$]", fontsize=28.4)
    ax.tick_params(axis="y", labelsize=24)
    ax.grid(True, alpha=0.3, axis="y", zorder=0)
    ax.set_title(f"({tag})", fontsize=32.4, loc="left")
    if legend:
        _stack_legend(ax, fontsize=20, ncol=3)
    return pd.DataFrame(rows), n_by_scheme


def finish_composites(axes, n_lists, legend=True):
    """Common y-range for the composites, headroom for the legend.  The
    per-bin sample counts are in the *_state_composites.csv table, not drawn."""
    lo = min(a.get_ylim()[0] for a in axes)
    hi = max(a.get_ylim()[1] for a in axes)
    lo = lo - (0.40 if legend else 0.03) * (hi - lo)
    for a in axes:
        a.set_ylim(lo, hi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="dcre_csv")
    ap.add_argument("--state-cache",
                    default=os.path.join("state_cache",
                                         "arg_vs_ml_consequences_rwp_v11_series_CCN2.npz"))
    ap.add_argument("--stem", default="dCRE_combined_4panel_MLvsARGv4")
    ap.add_argument("--figdir", default=".",
                    help="where the figure and its tables are written")
    ap.add_argument("--tend", default="2017-07-15 13:00")
    ap.add_argument("--backgrounds", default="org,mid,mid2",
                    help="comma list of background rows to draw AND to pool "
                         "into the composites, e.g. org,mid (v5)")
    ap.add_argument("--composite-backgrounds", default=None,
                    help="comma list of backgrounds pooled into the composites; "
                         "defaults to --backgrounds.  v6 = --backgrounds org,mid "
                         "--composite-backgrounds org,mid,mid2")
    ap.add_argument("--composites", default="Nd,LWP,CF",
                    help="comma list of state composites to draw, in order "
                         "(subset of Nd,LWP,CF)")
    ap.add_argument("--stack-legend", default="all",
                    choices=["all", "first", "top", "none"],
                    help="component legend in every stacked panel, only in the "
                         "first row's, as ONE strip above the whole figure "
                         "(v6), or nowhere (composites carry one only under "
                         "'all')")
    a = ap.parse_args()
    tend = pd.Timestamp(a.tend) if a.tend else None
    comp_want = [c.strip() for c in a.composites.split(",") if c.strip()]
    comp_specs = [next(sp for sp in STATE_SPECS if sp[0] == c) for c in comp_want]
    want = [b.strip() for b in a.backgrounds.split(",") if b.strip()]
    pool = ([b.strip() for b in a.composite_backgrounds.split(",") if b.strip()]
            if a.composite_backgrounds else list(want))
    known = [r[0] for r in ROWS]
    for b in want + pool:
        if b not in known:
            raise SystemExit(f"unknown background {b!r}; known: {known}")
    rows = [next(r for r in ROWS if r[0] == b) for b in want]
    # every background that is drawn or pooled has to be loaded
    load_rows = [r for r in ROWS if r[0] in want or r[0] in pool]

    state = load_state(a.state_cache)

    # geometry in inches, so 2 and 3 background rows keep the same proportions
    # as the original 3-row layout (24 x 21.5 in)
    n_rows = len(rows)
    row_h, comp_h, gap_h, top_m, bot_m = 4.28, 4.7, 2.3, 0.55, 1.1
    if a.stack_legend == "top":
        top_m += 1.55                      # legend strip, two rows at 23 pt
    fig_h = top_m + row_h * n_rows + gap_h + comp_h + bot_m
    fig = plt.figure(figsize=(24.0, fig_h))
    gs_top = GridSpec(n_rows, 2, figure=fig, width_ratios=[2, 1],
                      hspace=0.50, wspace=0.20,
                      top=1 - top_m / fig_h,
                      bottom=(bot_m + comp_h + gap_h) / fig_h,
                      left=0.065, right=0.985)
    gs_bot = GridSpec(1, len(comp_specs), figure=fig, wspace=0.28,
                      top=(bot_m + comp_h) / fig_h, bottom=bot_m / fig_h,
                      left=0.065, right=0.985)

    letters = "abcdefghijklmnop"
    tags = [[letters[2 * i], letters[2 * i + 1]] for i in range(n_rows)]
    comp_tags = letters[2 * n_rows:2 * n_rows + len(comp_specs)]
    merged, samples = [], []
    for bg, color, ml_stub, arg_stub, ml_cases, arg_cases in load_rows:
        df_ml  = load(a.outdir, ml_stub, tend)
        df_arg = load(a.outdir, arg_stub, tend)

        if bg in want:
            irow = want.index(bg)
            draw_stacked_timeseries_pair(
                fig.add_subplot(gs_top[irow, 0]), df_ml, df_arg, tags[irow][0],
                legend=(a.stack_legend == "all" or
                        (a.stack_legend == "first" and irow == 0)))
            draw_ml_arg_bars(fig.add_subplot(gs_top[irow, 1]), df_ml, df_arg,
                             tags[irow][1], color)
        if bg not in pool:
            continue

        # join the pair-mean state to every output time of the dCRE record
        for df, scheme, cases in ((df_ml, "ML", ml_cases), (df_arg, "ARG", arg_cases)):
            st = pair_state(state, cases)
            m = df.join(st, how="inner")
            if len(m) != len(df):
                print(f"  note: {scheme} {bg}: {len(df) - len(m)} dCRE times "
                      f"have no state record and are dropped from the composites")
            m.insert(0, "background", bg)
            m.insert(0, "scheme", scheme)
            merged.append(m)
            samples.append(window_samples(m, scheme, bg))

    merged = pd.concat(merged)
    samples = pd.concat(samples)
    os.makedirs(a.figdir, exist_ok=True)
    merged.to_csv(os.path.join(a.figdir, f"{a.stem}_state_merged_2min.csv"))
    samples.to_csv(os.path.join(a.figdir, f"{a.stem}_state_samples_30min.csv"))

    axes, n_lists, tables = [], [], []
    if a.stack_legend == "top":
        fig.legend(handles=_stack_handles(), fontsize=23, frameon=True, ncol=5,
                   loc="upper center", bbox_to_anchor=(0.525, 0.995),
                   columnspacing=1.6, handlelength=2.2)
    comp_legend = a.stack_legend == "all"
    for icol, ((xcol, _, xlabel, edges, fmt), tag) in enumerate(
            zip(comp_specs, comp_tags)):
        ax = fig.add_subplot(gs_bot[0, icol])
        tab, nby = draw_state_composite(ax, samples, xcol, edges, fmt, tag,
                                        xlabel, legend=(comp_legend and icol == 0))
        axes.append(ax); n_lists.append(nby); tables.append(tab)
    finish_composites(axes, n_lists, legend=comp_legend)
    comp = pd.concat(tables, ignore_index=True)
    comp.to_csv(os.path.join(a.figdir, f"{a.stem}_state_composites.csv"), index=False)

    with pd.option_context("display.width", 220, "display.max_columns", 30,
                           "display.float_format", lambda v: f"{v:.1f}"):
        print(comp.to_string(index=False))

    for ext in ("png", "pdf", "eps"):
        out = os.path.join(a.figdir, f"{a.stem}.{ext}")
        fig.savefig(out, dpi=200, bbox_inches="tight")
        print(f"  wrote {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
