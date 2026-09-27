#!/usr/bin/env python3
"""
plot_dcre_argminusml_30min.py

Companion to the dCRE 4-panel figure: WHEN does each channel favour ARG-ACT?

Top row   : ARG-ACT minus ML-ACT, per 30-min bin, for every dCRE term and the
            total (thick black, +-1 SE combined across the two schemes).
            Positive = ARG-ACT cools less; a line below zero = ARG-ACT's term is
            the larger cooling.  Dotted verticals mark the first bin in which the
            LWP term and the total go negative and stay negative.
Bottom row: ARG-ACT / ML-ACT ratio of the cloud state (control runs solid,
            3aer runs dashed): in-cloud N_d, LWP and rain water path, from the
            consequences cache.  Dotted verticals mark the 2-min output at which
            the control LWP ratio drops through 1 and the RWP ratio rises
            through 1.

One column per background (default org,mid).  Same 30-min bins, same 09:04-
13:00 window, same CSVs and cache as plot_dcre_combined_MLvsARGv4.py.
"""
import argparse, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

TERMS = [("DCRE_CF",    "#1f77b4", r"$\Delta$CRE$_{\rm CF}$",   1.6),
         ("DCRE_A_Nc",  "#ff7f0e", r"$\Delta$CRE$_{N_d}$",      1.6),
         ("DCRE_A_LWP", "#2ca02c", r"$\Delta$CRE$_{\rm LWP}$",  3.2),
         ("DCRE_A_re",  "#d62728", r"$\Delta$CRE$_{r_e}$",      1.6),
         ("DCRE_A_cov", "#9467bd", r"$\Delta$CRE$_{\rm cov}$",  1.6)]
STATE = [("QNDROP",  "#ff7f0e", r"$N_d$"),
         ("LWP_TOT", "#2ca02c", "LWP"),
         ("RWP",     "#17becf", "RWP")]
ROWS = {"org":  ("org_3aer_org",   "ARG_org_3aer_ARG_org",   ("org", "org_3aer"),   ("ARG_org", "ARG_org_3aer")),
        "mid":  ("mid_3aer_mid",   "ARG_mid_3aer_ARG_mid",   ("mid", "mid_3aer"),   ("ARG_mid", "ARG_mid_3aer")),
        "mid2": ("mid2_3aer_mid2", "ARG_mid2_3aer_ARG_mid2", ("mid2", "mid2_3aer"), ("ARG_mid2", "ARG_mid2_3aer"))}
BIN = "30min"; MIN_TIMES = 2
FS = dict(title=30, label=26, tick=22, legend=20, note=18)


def load(outdir, stub, tend):
    df = pd.read_csv(os.path.join(outdir, f"dcre_pixelwise_{stub}_domain.csv"), parse_dates=["time"]).set_index("time").sort_index()
    return df[df.index <= tend]


def state_series(z, case, key, tend):
    t, s = z[case]
    x = pd.Series(np.asarray(s[key], float), index=pd.DatetimeIndex(t)).sort_index()
    return x[x.index <= tend]


def first_persistent(sign_series, want):
    """First index from which sign_series == want holds through the end."""
    ok = (sign_series == want).values
    if ok.all():
        return sign_series.index[0]
    bad = np.where(~ok)[0]
    return sign_series.index[bad[-1] + 1] if bad[-1] + 1 < len(ok) else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default="dcre_csv")
    ap.add_argument("--state-cache", default="state_cache/arg_vs_ml_consequences_rwp_v11_series_CCN2.npz")
    ap.add_argument("--backgrounds", default="org,mid,mid2",
                    help="columns; 'mid+mid2' pools the two moist pairs (2-min records concatenated)")
    ap.add_argument("--stem", default="dCRE_ARGminusML_30min")
    ap.add_argument("--tend", default="2017-07-15 13:00")
    a = ap.parse_args()
    tend = pd.Timestamp(a.tend)
    bgs = [b.strip() for b in a.backgrounds.split(",")]
    z = np.load(a.state_cache, allow_pickle=True)["data"].item()

    fig, axes = plt.subplots(2, len(bgs), figsize=(9.5 * len(bgs), 15.0),
                             gridspec_kw=dict(height_ratios=[1.35, 1.0], hspace=0.30, wspace=0.22),
                             squeeze=False)
    letters = "abcdefgh"
    rows_out = []
    for j, bg in enumerate(bgs):
        parts = bg.split("+")
        ms, rs = [], []
        for pb in parts:
            ml_stub, arg_stub, _, _ = ROWS[pb]
            m1, r1 = load(a.outdir, ml_stub, tend), load(a.outdir, arg_stub, tend)
            both1 = m1.index.intersection(r1.index)
            ms.append(m1.loc[both1]); rs.append(r1.loc[both1])
        # pooled backgrounds: concatenate their 2-min records (a bin then holds
        # both pairs' samples, n doubled); the state ratios are averaged over
        # the pooled backgrounds
        m, r = pd.concat(ms), pd.concat(rs)
        both = m.index
        key = both.floor(BIN)
        n = pd.Series(1, index=both).groupby(key).size()
        keep = n[n >= MIN_TIMES * len(parts)].index
        centre = keep + pd.Timedelta(BIN) / 2
        ml_cases = tuple(ROWS[pb][2] for pb in parts); arg_cases = tuple(ROWS[pb][3] for pb in parts)

        # ---- top: ARG - ML per term, per bin --------------------------------
        ax = axes[0, j]
        first_sig = {}
        for col, c, lab, lw in TERMS + [("DCRE_total", "black", r"$\Delta$CRE$_{\rm total}$", 3.2)]:
            d = (r[col].values - m[col].values)
            d = pd.Series(d, index=both)
            g = d.groupby(key)
            mean = g.mean().loc[keep]
            se = np.sqrt(m[col].groupby(key).var() / n + r[col].groupby(key).var() / n).loc[keep]
            if col == "DCRE_total":
                ax.errorbar(centre, mean, yerr=se, color=c, lw=lw, marker="o", ms=8, capsize=4, zorder=6, label=lab)
            else:
                ax.plot(centre, mean, color=c, lw=lw, marker="o", ms=7 if lw > 2 else 6, zorder=5, label=lab)
            # first bin from which the difference is negative by more than 1 SE and stays so
            sig = pd.Series(np.where(mean + se < 0, -1, np.where(mean - se > 0, 1, 0)), index=keep)
            first_sig[col] = first_persistent(sig, -1)
            for t, v, e in zip(mean.index, mean.values, se.values):
                rows_out.append(dict(background=bg, bin=t, term=col, ARG_minus_ML=v, SE=e, n_ML=int(n[t]), n_ARG=int(n[t])))
        ax.axhline(0, color="k", lw=1.0, zorder=3)
        notes = []
        for col, c in (("DCRE_A_LWP", "#2ca02c"), ("DCRE_total", "black")):
            t0 = first_sig[col]
            if t0 is not None:
                ax.axvline(t0 + pd.Timedelta(BIN) / 2, color=c, ls=":", lw=2.2, zorder=2)
                notes.append((f"{lab_short(col)} < 0 beyond 1 SE from the {t0:%H:%M} bin", c))
        if notes:                                   # room under the data for the notes
            lo, hi = ax.get_ylim(); ax.set_ylim(lo - 0.20 * (hi - lo), hi)
        _note_box(ax, notes)
        ax.set_title(f"({letters[j]}) {bg.replace('+', ' + ')}", fontsize=FS["title"], loc="left")
        ax.set_ylabel(r"ARG-ACT $-$ ML-ACT  [W m$^{-2}$]", fontsize=FS["label"])
        ax.tick_params(labelsize=FS["tick"]); ax.grid(True, alpha=0.3, zorder=0)
        if j == 0:
            h_top, l_top = ax.get_legend_handles_labels()

        # ---- bottom: state ratios ARG/ML --------------------------------------
        ax2 = axes[1, j]
        notes2 = []
        for key_s, c, lab in STATE:
            for member, ls, tag in ((0, "-", "control"), (1, "--", "3aer")):
                xs = [state_series(z, ac[member], key_s, tend) / state_series(z, mc[member], key_s, tend)
                      for mc, ac in zip(ml_cases, arg_cases)]
                x = pd.concat(xs, axis=1).mean(axis=1) if len(xs) > 1 else xs[0]
                ax2.plot(x.index, x.values, color=c, ls=ls, lw=2.6 if ls == "-" else 1.8,
                         label=f"{lab} ({tag})" if j == 0 else None, zorder=4)
                if tag == "control" and key_s in ("LWP_TOT", "RWP"):
                    want = -1 if key_s == "LWP_TOT" else 1
                    t0 = first_persistent(np.sign(x - 1), want)
                    if t0 is not None and t0 != x.index[0]:
                        ax2.axvline(t0, color=c, ls=":", lw=2.2, zorder=2)
                        notes2.append((f"{lab} (control) {'<' if want < 0 else '>'} 1 from {t0:%H:%M}", c))
        _note_box(ax2, notes2)
        ax2.axhline(1, color="k", lw=1.0, zorder=3)
        ax2.set_yscale("log")
        ax2.set_yticks([0.05, 0.1, 0.2, 0.5, 1, 2]); ax2.set_yticklabels(["0.05", "0.1", "0.2", "0.5", "1", "2"])
        ax2.set_ylim(0.03, 2.2)
        ax2.set_title(f"({letters[len(bgs) + j]}) {bg.replace('+', ' + ')}", fontsize=FS["title"], loc="left")
        ax2.set_ylabel("ARG-ACT / ML-ACT", fontsize=FS["label"])
        ax2.set_xlabel("time (UTC)", fontsize=FS["label"])
        ax2.tick_params(labelsize=FS["tick"]); ax2.grid(True, alpha=0.3, which="both", zorder=0)
        if j == 0:
            h_bot, l_bot = ax2.get_legend_handles_labels()

        for axx in (ax, ax2):
            edges = keep.append(pd.DatetimeIndex([keep[-1] + pd.Timedelta(BIN)]))
            axx.set_xticks(edges); axx.set_xticklabels([f"{e:%H:%M}" for e in edges], rotation=35, ha="right")
            axx.set_xlim(edges[0] - 0.15 * pd.Timedelta(BIN), edges[-1] + 0.35 * pd.Timedelta(BIN))

    fig.legend(h_top, l_top, loc="upper center", bbox_to_anchor=(0.5, 1.045), ncol=6,
               fontsize=FS["legend"], frameon=True,
               title="top row: ARG-ACT minus ML-ACT per 30-min bin;  (+) ARG-ACT cools less, (–) ARG-ACT cools more",
               title_fontsize=FS["legend"])
    fig.legend(h_bot, l_bot, loc="lower center", bbox_to_anchor=(0.5, -0.085), ncol=6,
               fontsize=FS["legend"], frameon=True,
               title="bottom row: ratio of the two treatments' domain cloud state (in-cloud $N_d$, LWP over cloudy columns, rain water path)",
               title_fontsize=FS["legend"])
    pd.DataFrame(rows_out).to_csv(f"{a.stem}_table.csv", index=False)
    for ext in ("png", "pdf"):
        fig.savefig(f"{a.stem}.{ext}", dpi=200, bbox_inches="tight"); print("  wrote", f"{a.stem}.{ext}")


def _note_box(ax, notes):
    """Stack short coloured notes in the panel's lower-left corner."""
    for i, (txt, c) in enumerate(notes):
        ax.text(0.02, 0.03 + 0.075 * i, txt, transform=ax.transAxes, color=c,
                fontsize=FS["note"], va="bottom", ha="left", zorder=8,
                bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))


def lab_short(col):
    return {"DCRE_A_LWP": "LWP term", "DCRE_total": "total"}[col]


if __name__ == "__main__":
    main()
