#!/usr/bin/env python
"""
plot_arg_relerr_map_v3.py  (v3: v2 split into TWO 2x2 figures, one per mode;
                            no per-panel stats box; --box_fs N restores it)
=============================================================================
Companion to plot_arg_failure_map_bwr.py, in the metric the activation
literature uses.

    colour         = MEDIAN RELATIVE error in the activated fraction
                     (and hence in activated number, since N_a cancels)

                         delta_rel = (fn_ARG - fn_parcel) / fn_parcel

    contour lines  = the same FAILURE RATE as the companion figure,
                     |fn_ARG - fn_parcel| > tau,  tau = 0.1

Why a separate figure.  The companion figure colours the ABSOLUTE difference
fn_ARG - fn_parcel, whose magnitude is capped by fn_parcel itself: where the
parcel model activates almost nothing, even a total misprediction is a small
absolute number.  Phinney et al. (2003), Ghan et al. (2011) and Rothenberg &
Wang (2016) all report RELATIVE error in Smax / Nd, so the absolute map cannot
be compared with them directly.  This figure closes that gap.

Conventions
-----------
* delta_rel is undefined where the parcel model does not activate, so bins are
  built only from samples with fn_parcel > FN_MIN (default 0.01).  Panels
  report how many samples survive that cut -- for the Aitken mode it is a small
  minority, which is itself the explanation for that row's low absolute
  failure rate.
* The per-bin statistic is the MEDIAN, not the mean: delta_rel is bounded below
  by -1 but unbounded above, so its mean is dragged upward by a thin tail of
  large over-activations and is not a good measure of the typical cell.  The
  stats file reports both.
* Colour scale is fixed at +/-1 (blue = ARG activates less; -1 = ARG activates
  nothing the parcel model activates; red = ARG activates more) so the two
  modes and all four panel pairs are directly comparable.  Bins above +1 are
  shown at the top colour and the bar is drawn with an arrow.

Writes
    <out_prefix>_aitken_map_filtered.{png,eps}   panels (a)-(d) of v2
    <out_prefix>_accum_map_filtered.{png,eps}    panels (e)-(h) of v2
    <out_prefix>_map_filtered_stats.txt
Each figure is lettered (a)-(d); pass --letters continue to keep v2's (e)-(h)
on the Accum figure.

Usage
-----
  $HOME/.conda/envs/liranenv_gpu/bin/python plot_arg_relerr_map.py
  ... --tau 0.1 --nbins 24 --min_cnt 20 --min_cnt_ait 10 --fn_min 0.01

The bin-occupancy threshold is per mode: --min_cnt for the Accum row, and the
looser --min_cnt_ait for the Aitken row, where only a few per cent of the pool
survives the fn_parcel > fn_min cut and the Accum threshold would blank almost
every bin.
=============================================================================
"""
import os
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from plot_allcases_arg import CASES, MODES, load_case

HERE = os.path.dirname(os.path.abspath(__file__))

AX = {
    "wbar":   (3,  r"$\bar{w}$ (cm s$^{-1}$)",                    True),
    "na_acc": (5,  r"$N_{\mathrm{a,accum}}$ (cm$^{-3}$)",          True),
    "na_ait": (4,  r"$N_{\mathrm{a,Aitken}}$ (cm$^{-3}$)",         True),
    "pres":   (1,  "pressure (hPa)",                               False),
    "tair":   (0,  "T (K)",                                        False),
}
PAIRS = [("wbar", "na_acc"), ("wbar", "na_ait"), ("na_acc", "na_ait"), ("wbar", "pres")]
LETTERS = "abcdefgh"
LEVELS = [0.10, 0.50, 0.90]

CMAP = plt.get_cmap("RdBu_r").copy()
CMAP.set_bad(alpha=0.0)


def binned_median(x, y, v, xe, ye, min_cnt):
    """Per-bin median of v, masked where the bin holds < min_cnt samples.

    Done by hand rather than with scipy so the script has no extra dependency:
    digitize into a flat bin id, sort by it, and take the median of each run."""
    nx, ny = len(xe) - 1, len(ye) - 1
    ix = np.clip(np.digitize(x, xe) - 1, 0, nx - 1)
    iy = np.clip(np.digitize(y, ye) - 1, 0, ny - 1)
    inside = (x >= xe[0]) & (x <= xe[-1]) & (y >= ye[0]) & (y <= ye[-1])
    ix, iy, v = ix[inside], iy[inside], v[inside]

    flat = ix * ny + iy
    order = np.argsort(flat, kind="stable")
    flat, vs = flat[order], v[order]
    starts = np.searchsorted(flat, np.arange(nx * ny), side="left")
    ends = np.searchsorted(flat, np.arange(nx * ny), side="right")

    M = np.full(nx * ny, np.nan)
    cnt = (ends - starts).astype(float)
    for b in np.where(cnt >= min_cnt)[0]:
        M[b] = np.median(vs[starts[b]:ends[b]])
    return (np.ma.masked_invalid(M.reshape(nx, ny)), cnt.reshape(nx, ny))


def binned_mean(x, y, v, xe, ye, min_cnt):
    cnt, _, _ = np.histogram2d(x, y, bins=[xe, ye])
    tot, _, _ = np.histogram2d(x, y, bins=[xe, ye], weights=v)
    with np.errstate(invalid="ignore", divide="ignore"):
        m = tot / cnt
    return np.ma.masked_where(cnt < min_cnt, m)


def smooth_nan(M):
    """3x3 box mean ignoring masked cells, for legible contour lines."""
    A = np.where(M.mask, np.nan, M.filled(np.nan)) if np.ma.isMaskedArray(M) else np.asarray(M)
    acc = np.zeros_like(A, float)
    num = np.zeros_like(A, float)
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            S = np.roll(np.roll(A, dx, 0), dy, 1)
            if dx == 1:
                S[0, :] = np.nan
            elif dx == -1:
                S[-1, :] = np.nan
            if dy == 1:
                S[:, 0] = np.nan
            elif dy == -1:
                S[:, -1] = np.nan
            ok = np.isfinite(S)
            acc[ok] += S[ok]
            num[ok] += 1.0
    with np.errstate(invalid="ignore", divide="ignore"):
        out = acc / num
    return np.ma.masked_invalid(np.where(np.isfinite(A), out, np.nan))


def draw_mode(X, R, F, USE, r, mode, nbins, mc, letters, out_prefix, box_fs=0):
    """One 2x2 figure for a single mode (row r of the v2 layout)."""
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 13.0), constrained_layout=True)
    for li, (kx, ky) in enumerate(PAIRS):
        ax = axes.ravel()[li]
        jx, lx, logx = AX[kx]
        jy, ly, logy = AX[ky]
        x, y = X[:, jx], X[:, jy]

        fin = np.isfinite(x) & np.isfinite(y)
        if logx:
            fin &= x > 0
        if logy:
            fin &= y > 0
        xx = np.log10(x[fin]) if logx else x[fin]
        yy = np.log10(y[fin]) if logy else y[fin]

        # bin edges from ALL samples so the axes match the companion figure
        xe = np.linspace(np.percentile(xx, 0.5), np.percentile(xx, 99.5), nbins + 1)
        ye = np.linspace(np.percentile(yy, 0.5), np.percentile(yy, 99.5), nbins + 1)

        # colour: median relative error, only over samples where it is defined
        u = USE[fin, r]
        M, _ = binned_median(xx[u], yy[u], R[fin, r][u], xe, ye, mc)
        # contours: failure rate over ALL samples, as in the companion figure
        Rate = binned_mean(xx, yy, F[fin, r], xe, ye, mc)

        ax.set_facecolor("0.86")
        pc = ax.pcolormesh(xe, ye, M.T, cmap=CMAP, vmin=-1.0, vmax=1.0,
                           shading="flat", rasterized=True)

        xc = 0.5 * (xe[:-1] + xe[1:])
        yc = 0.5 * (ye[:-1] + ye[1:])
        S = smooth_nan(Rate)
        lv = [l for l in LEVELS if np.nanmin(S) < l < np.nanmax(S)] if S.count() else []
        if lv:
            cs = ax.contour(xc, yc, S.T, levels=lv, colors="k",
                            linewidths=[1.0 + 0.5 * (l == 0.5) for l in lv],
                            linestyles="-", zorder=4)
            ax.clabel(cs, inline=True, inline_spacing=4, fontsize=18, fmt="%.2f")

        ax.text(0.5, 1.02, f"({letters[li]})", transform=ax.transAxes,
                fontsize=30, fontweight="bold", va="bottom", ha="center")
        med = float(np.median(R[fin, r][u])) if u.any() else np.nan
        # Accum panels fill the upper left, so drop the box to the bottom left there
        ty, tva = (0.03, "bottom") if mode == "Accum" else (0.97, "top")
        if box_fs > 0:
          ax.text(0.03, ty,
                f"N usable={int(u.sum()):,} / {int(fin.sum()):,}\n"
                  f"bins shown={int(M.count())}/{nbins*nbins} (>= {mc}/bin)\n"
                  f"median rel err={med:+.3f}\n"
                  f"frac ARG<parcel={np.mean(R[fin, r][u] < 0) if u.any() else np.nan:.3f}",
                  transform=ax.transAxes, va=tva, fontsize=box_fs, zorder=6,
                  bbox=dict(fc="white", alpha=0.92, ec="0.7"))
        ax.set_xlabel((r"log$_{10}$ " if logx else "") + lx, fontsize=28)
        ax.set_ylabel((r"log$_{10}$ " if logy else "") + ly, fontsize=28)
        ax.tick_params(axis="both", labelsize=24)
        # pin tick spacing: the auto locator thins ticks at this font size
        ax.xaxis.set_major_locator(MultipleLocator(1 if logx else 25))
        ax.yaxis.set_major_locator(MultipleLocator(1 if logy else 25))

    cb = fig.colorbar(pc, ax=axes.ravel().tolist(), fraction=0.035, pad=0.010,
                      aspect=40, extend="max")
    cb.ax.tick_params(labelsize=24)
    cb.set_label(r"$\delta_{\mathrm{rel}}$", fontsize=28, labelpad=8)

    for ext in ("png", "eps"):
        p = f"{out_prefix}.{ext}"
        fig.savefig(p, dpi=300, bbox_inches="tight")
        print(f"Wrote {p}", flush=True)
    plt.close(fig)


def draw(X, R, F, USE, tau, nbins, min_cnt, fn_min, npool, out_prefix, title_extra,
         letters_mode="restart", box_fs=0):
    """Same signature as v2 (tau/fn_min/npool/title_extra unused, kept for the
    call site); emits one figure per mode instead of one 2x4 figure."""
    for r, mode in enumerate(MODES):
        letters = LETTERS[4 * r:4 * r + 4] if letters_mode == "continue" else LETTERS[:4]
        draw_mode(X, R, F, USE, r, mode, nbins, min_cnt[mode], letters,
                  os.path.join(os.path.dirname(out_prefix),
                               os.path.basename(out_prefix).replace(
                                   "_map_filtered", f"_{mode.lower()}_map_filtered")),
                  box_fs=box_fs)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cases", nargs="*", default=None)
    ap.add_argument("--tau", type=float, default=0.1, help="failure threshold (contours)")
    ap.add_argument("--nbins", type=int, default=24)
    ap.add_argument("--min_cnt", type=int, default=20,
                    help="minimum usable samples for a bin to be coloured (Accum row)")
    ap.add_argument("--min_cnt_ait", type=int, default=10,
                    help="same, for the Aitken row.  Only ~6%% of the pool has "
                         "fn_parcel > fn_min in the Aitken mode, so the Accum "
                         "threshold blanks nearly the whole row; default 10.")
    ap.add_argument("--fn_min", type=float, default=0.01,
                    help="relative error is only defined where fn_parcel exceeds this")
    ap.add_argument("--out_prefix", default=os.path.join(HERE, "arg_relerr_v3"))
    ap.add_argument("--letters", choices=["restart", "continue"], default="restart",
                    help="restart: (a)-(d) on both figures; continue: (e)-(h) on Accum")
    ap.add_argument("--box_fs", type=float, default=0,
                    help="font size of the per-panel stats box; 0 (default) = no box. v2 used 11")
    ap.add_argument("--callsite", choices=["all", "old_cloud", "grow_shrink"],
                    default="all",
                    help="restrict to rows produced by one activate_ml call site. "
                         "old_cloud = the site that actually wrote FN13/FN23, so the "
                         "stored ARG is synchronous there; grow_shrink = the site that "
                         "did not, so the stored ARG is stale. Requires 'callsite' in "
                         "the npz (see arg_offline/add_callsite.py).")
    args = ap.parse_args()

    cases = [c for c in CASES if args.cases is None or c[0] in args.cases]
    want = {"old_cloud": 2, "grow_shrink": 1}.get(args.callsite)

    Xs, Ps, As, per_case = [], [], [], []
    for cdir, label in cases:
        print(f"loading {cdir} ...", flush=True)
        c = load_case(os.path.join(HERE, cdir))
        k = c["keep"]
        z = np.load(os.path.join(HERE, cdir, "online", "parcel_inputs.npz"))
        X = z["X"][k]
        P = c["fn_parcel"][k]; A = c["fn_arg"][k]
        if want is not None:
            if "callsite" not in z.files:
                raise KeyError(f"{cdir}: no 'callsite' in parcel_inputs.npz -- "
                               "run arg_offline/add_callsite.py first")
            sel = z["callsite"][k] == want
            X, P, A = X[sel], P[sel], A[sel]
            print(f"    callsite={args.callsite}: kept {int(sel.sum()):,}"
                  f" / {int(sel.size):,}", flush=True)
        Xs.append(X); Ps.append(P); As.append(A)
        per_case.append((label, int(X.shape[0])))
        print(f"    kept={X.shape[0]:,}", flush=True)
    X = np.concatenate(Xs, 0); P = np.concatenate(Ps, 0); A = np.concatenate(As, 0)
    print(f"pooled: {X.shape[0]:,}", flush=True)

    USE = P > args.fn_min                       # relative error defined here
    with np.errstate(invalid="ignore", divide="ignore"):
        R = np.where(USE, (A - P) / np.where(P > 0, P, np.nan), np.nan)
    F = (np.abs(A - P) > args.tau).astype(float)

    extra = ("\nCAVEAT: fn_ARG (FN13/FN23) is a persistent WRF array; on cells last served "
             "by the GROW_SHRINK call site it is stale, so part of this error is "
             "time-mismatch, not scheme error.")
    min_cnt = {m: (args.min_cnt_ait if m == "Aitken" else args.min_cnt) for m in MODES}
    draw(X, R, F, USE, args.tau, args.nbins, min_cnt, args.fn_min,
         X.shape[0], args.out_prefix + "_map_filtered", extra,
         letters_mode=args.letters, box_fs=args.box_fs)

    # ── stats, including the regimes the activation literature quotes ──
    W_cms, NAIT, NACC = X[:, 3], X[:, 4], X[:, 5]
    W_ms = W_cms / 100.0
    rep = args.out_prefix + "_map_filtered_stats.txt"
    with open(rep, "w") as fh:
        fh.write("=" * 104 + "\nARG RELATIVE ERROR vs THE FORTRAN PARCEL MODEL\n" + "=" * 104 + "\n\n")
        fh.write("relative error := (fn_ARG - fn_parcel) / fn_parcel, "
                 f"taken only where fn_parcel > {args.fn_min:g}\n")
        fh.write(f"failure (contours) := |fn_ARG - fn_parcel| > {args.tau:g}\n")
        fh.write(f"bins: {args.nbins} x {args.nbins}, blank below "
                 + ", ".join(f"{min_cnt[m]} usable samples ({m})" for m in MODES) + "\n\n")
        fh.write(f"{'case':<14} {'kept':>9}\n")
        for label, n in per_case:
            fh.write(f"{label:<14} {n:>9,}\n")
        fh.write(f"{'POOLED':<14} {X.shape[0]:>9,}\n\n")

        for m, mode in enumerate(MODES):
            u = USE[:, m]
            r = R[u, m]
            fh.write(f"{mode}: usable={int(u.sum()):,} ({u.mean():.3f} of pool)  "
                     f"median rel={np.median(r):+.4f}  mean rel={r.mean():+.4f}  "
                     f"frac(ARG<parcel)={np.mean(A[u, m] < P[u, m]):.4f}\n")
            fh.write(f"    fn_parcel mean={P[:, m].mean():.4f}  frac<0.01={np.mean(P[:, m] < 0.01):.4f}"
                     f"   fn_ARG mean={A[:, m].mean():.4f}  frac<0.01={np.mean(A[:, m] < 0.01):.4f}"
                     f"  frac exactly 0={np.mean(A[:, m] == 0.0):.4f}\n")

        m = 1
        u = USE[:, m]
        fh.write("\n" + "=" * 104 + "\n")
        fh.write("ACCUM, in the coordinate the literature uses: N_accum / w  (cm-3 per m s-1)\n")
        fh.write("Phinney et al. (2003): ARG fails for V < 50 cm/s and Na > 500 cm-3.\n")
        fh.write("Ghosh et al. (2025, GMD 18, 4899): the kinetically limited regime is N/w > 1e4.\n")
        fh.write("=" * 104 + "\n")
        ratio = NACC / np.maximum(W_ms, 1e-12)
        fh.write(f"{'N/w bin':>26} {'n':>8} {'fn_parcel':>10} {'fn_ARG':>9} {'abs bias':>10} "
                 f"{'median rel':>11} {'mean rel':>10} {'ARG<parcel':>11}\n")
        edges = [0, 1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6, np.inf]
        for b in range(len(edges) - 1):
            s = (ratio >= edges[b]) & (ratio < edges[b + 1])
            su = s & u
            if s.sum() < 50 or su.sum() < 10:
                continue
            fh.write(f"[{edges[b]:>10.3g},{edges[b+1]:>10.3g}) {int(s.sum()):>8,} "
                     f"{P[s, m].mean():>10.4f} {A[s, m].mean():>9.4f} "
                     f"{(A[s, m]-P[s, m]).mean():>+10.4f} {np.median(R[su, m]):>+11.4f} "
                     f"{R[su, m].mean():>+10.4f} {np.mean(A[su, m] < P[su, m]):>11.4f}\n")

        fh.write("\n" + "=" * 104 + "\n")
        fh.write("ACCUM, the named regimes\n" + "=" * 104 + "\n")
        for name, s in (("all cells", np.ones(X.shape[0], bool)),
                        ("Phinney regime: V<50 cm/s & Na_acc>500 cm-3",
                         (W_cms < 50) & (NACC > 500)),
                        ("Phinney 'good' regime: V>50 & Na_acc<500",
                         (W_cms >= 50) & (NACC <= 500))):
            su = s & u
            if su.sum() < 10:
                continue
            nd = (A[su, m] * NACC[su]).mean() / (P[su, m] * NACC[su]).mean()
            fh.write(f"  {name:<44} n={int(su.sum()):>7,}  median rel={np.median(R[su, m]):+.4f}  "
                     f"mean rel={R[su, m].mean():+.4f}  frac(ARG<parcel)={np.mean(A[su, m] < P[su, m]):.4f}  "
                     f"mean Nact ARG/parcel={nd:.4f}\n")

        fh.write("\n" + "=" * 104 + "\n")
        fh.write("ACCUM: does ARG show the unphysical 'Nact falls as Na rises'? (w < 50 cm/s)\n")
        fh.write("=" * 104 + "\n")
        lo = W_cms < 50.0
        ed = np.percentile(NACC[lo], np.linspace(0, 100, 11))
        fh.write(f"{'Na_accum bin (cm-3)':>26} {'n':>8} {'Nact_parcel':>13} {'Nact_ARG':>11} {'ARG/parcel':>11}\n")
        for b in range(10):
            s = lo & (NACC >= ed[b]) & (NACC <= ed[b+1] if b == 9 else NACC < ed[b+1])
            if s.sum() < 50:
                continue
            ndp = (P[s, m] * NACC[s]).mean(); nda = (A[s, m] * NACC[s]).mean()
            fh.write(f"[{ed[b]:>10.4g},{ed[b+1]:>10.4g}) {int(s.sum()):>8,} {ndp:>13.2f} "
                     f"{nda:>11.2f} {nda/max(ndp,1e-9):>11.4f}\n")
    print(f"Wrote {rep}", flush=True)


if __name__ == "__main__":
    main()
