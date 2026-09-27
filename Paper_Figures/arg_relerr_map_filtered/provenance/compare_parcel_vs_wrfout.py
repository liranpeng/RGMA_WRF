#!/usr/bin/env python
"""
compare_parcel_vs_wrfout.py
=========================================================================
For every updraft cell in the wrfout_d01_* files, run the OFFLINE Fortran
parcel box model (libactivate.so, fortranactivate_) on the reconstructed
emulator inputs and compare its activation fraction against the emulator
fraction WRF actually applied and stored on disk (FN11 Aitken, FN21 Accum).

This is the "ground-truth physics vs deployed emulator" test, sourced entirely
from wrfout (no rsl.error logs). It answers: on the states WRF really visited,
how far does the applied ML emulator sit from the parcel model it was meant to
emulate?

Pipeline (mirrors retrain_wrfdist/extract_and_label.py exactly)
---------------------------------------------------------------
1. Read each wrfout, reconstruct the 16 emulator inputs as module_mixactivate.F
   builds input_array (L4311-4332): tair, pressure(hPa), rh(%), wbar(cm/s),
   na*1e-6 (4 modes), rmean(um) (4 modes), kappa (4 modes).
2. Mask to updraft cells WRF evaluated (EMTAIR>0 & EMWBAR>0) that also carry
   valid physical aerosol in the two active modes (NA11,NA21 > 1 cm-3, HG>0,
   VO>0) -- the parcel ODE needs physical modes.
3. Run the parcel model on the 2 physical modes (Aitken=NA11, Accum=NA21),
   parallel over cores, dynamic dispatch. -> fn_parcel[N,2].
4. Compare fn_parcel vs the stored emulator FN11/FN21 (fn_wrf). Print stats,
   write scatter + histogram PNGs.

Incomplete wrfout files (Time dim length 0, still being written) are skipped.

--------------------------------------------------------------------------
MULTI-NODE (--phase). multiprocessing is single-node, so to use N nodes we
split into 3 phases driven by srun (see submit_parcel_compare.slurm):

  reconstruct : 1 task  -> read wrfout ONCE, save X/sg2/fn_wrf to --inputs_npz
  parcel      : K tasks -> each node loads --inputs_npz, solves cells[shard::K],
                            saves fn_parcel for its shard to --shard_npz
  merge       : 1 task  -> gather all shard npz, plot + stats

Single-node default is --phase all (reconstruct+parcel+merge in one process).

Usage:
  python compare_parcel_vs_wrfout.py --files ../wrfout_d01_2017-07-15*        # single node
  python compare_parcel_vs_wrfout.py --all --workers 128                      # single node, all files
  # multi-node phases are wired up by submit_parcel_compare.slurm
=========================================================================
"""
import os
import glob
import time
import ctypes
import argparse
import multiprocessing as mp

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

HERE     = os.path.dirname(os.path.abspath(__file__))
TEST_DIR = os.path.dirname(HERE)
DEFAULT_LIB = "/scratch/07088/tg863871/Perlm_Backup/retrain_wrfdist_stam3/parcel_lib/libactivate.so"

THIRD = 1.0 / 3.0
PI    = np.pi
MODES = ["Aitken", "Accum"]
N_MIN_CM3 = 1.0        # physical-aerosol floor for the 2 active modes
TLIMIT    = 60         # parcel per-sample wallclock limit (s)

# fields needed to reconstruct inputs + the applied emulator outputs
NEED = ["EMTAIR", "EMPRES", "EMRH", "EMWBAR",
        "NA11", "NA21", "NA12", "NA22",
        "VO11", "VO21", "VO12", "VO22",
        "SG11", "SG21", "SG12", "SG22",
        "HG11", "HG21", "HG12", "HG22",
        "FN11", "FN21"]


# ── parcel-model ctypes interface (identical to extract_and_label.py) ──────────
_lib = None
def _setup_argtypes(lib):
    fn = lib.fortranactivate_
    fn.restype = None
    Pi = ctypes.POINTER(ctypes.c_int); Pf = ctypes.POINTER(ctypes.c_float)
    fn.argtypes = [Pi, Pf, Pf, Pf, Pf, Pf, Pf, Pf, Pf, Pi, Pi, Pi, Pf, Pf]

# globals shared with workers via fork (copy-on-write; set before Pool)
G_T = G_P = G_RH = G_W = None        # [N] scalars
G_NUM = G_SIG = G_KAP = G_RA = None  # [N,2] the two physical modes

def _pool_init(lib_path):
    global _lib
    devnull = os.open(os.devnull, os.O_WRONLY)
    os.dup2(devnull, 1); os.dup2(devnull, 2); os.close(devnull)   # silence Fortran prints
    _lib = ctypes.CDLL(lib_path)
    _setup_argtypes(_lib)

def _run_idx(i):
    """Run the 2-mode parcel model for sample i. Returns (i,fn0,fn1,smax,bad,tout)."""
    num_aer = np.ascontiguousarray(G_NUM[i], np.float32)
    sigmag  = np.ascontiguousarray(G_SIG[i], np.float32)
    kappa   = np.ascontiguousarray(G_KAP[i], np.float32)
    r_aer   = np.ascontiguousarray(G_RA[i],  np.float32)
    fn_out  = np.zeros(2, np.float32)
    Pf = ctypes.POINTER(ctypes.c_float)
    c_nm = ctypes.c_int(2)
    c_t  = ctypes.c_float(float(G_T[i]));  c_p  = ctypes.c_float(float(G_P[i]))
    c_rh = ctypes.c_float(float(G_RH[i])); c_w  = ctypes.c_float(float(G_W[i]))
    c_tl = ctypes.c_int(TLIMIT); c_bad = ctypes.c_int(0); c_tout = ctypes.c_int(0)
    c_smax = ctypes.c_float(0.0)
    _lib.fortranactivate_(
        ctypes.byref(c_nm), ctypes.byref(c_t), ctypes.byref(c_p), ctypes.byref(c_rh),
        num_aer.ctypes.data_as(Pf), ctypes.byref(c_w), sigmag.ctypes.data_as(Pf),
        kappa.ctypes.data_as(Pf), r_aer.ctypes.data_as(Pf), ctypes.byref(c_tl),
        ctypes.byref(c_bad), ctypes.byref(c_tout), fn_out.ctypes.data_as(Pf),
        ctypes.byref(c_smax))
    return (i, float(fn_out[0]), float(fn_out[1]), float(c_smax.value),
            int(c_bad.value), int(c_tout.value))


def rmean_um(VO, NA, SG):
    """WRF emulator-input radius: exp(-1.5 ln^2 SG)*(3 VO/(4 pi NA))^(1/3)*1e6, 0 if NA<=0."""
    out = np.zeros_like(NA, dtype=np.float64)
    m = NA > 0.0
    aln = np.log(SG[m])
    out[m] = np.exp(-1.5 * aln * aln) * (3.0 * VO[m] / (4.0 * PI * NA[m])) ** THIRD * 1e6
    return out


def read_wrfout(path):
    """Return (X16[n,16], sg2[n,2], fn_wrf[n,2]) for updraft cells with valid
       physical aerosol, or (None,None,None) if the file has no time record yet."""
    import netCDF4
    nc = netCDF4.Dataset(path, "r"); nc.set_auto_mask(False)
    if len(nc.dimensions["Time"]) == 0:               # still being written
        nc.close(); return None, None, None
    g = lambda v: nc.variables[v][0].astype(np.float64).ravel()   # [0] drops Time
    d = {v: g(v) for v in NEED}
    nc.close()

    mask = (d["EMTAIR"] > 0.0) & (d["EMWBAR"] > 0.0)  # updraft cells WRF evaluated
    for tag in ("11", "21"):                          # valid physical aerosol
        mask &= (d[f"NA{tag}"] * 1e-6 > N_MIN_CM3) & (d[f"HG{tag}"] > 0.0) & (d[f"VO{tag}"] > 0.0)
    idx = np.where(mask)[0]
    if idx.size == 0:
        return None, None, None
    for k in d:
        d[k] = d[k][idx]

    n = idx.size
    X = np.zeros((n, 16), dtype=np.float64)
    X[:, 0] = d["EMTAIR"]; X[:, 1] = d["EMPRES"]; X[:, 2] = d["EMRH"]; X[:, 3] = d["EMWBAR"]
    X[:, 4] = d["NA11"] * 1e-6; X[:, 5] = d["NA21"] * 1e-6
    X[:, 6] = d["NA12"] * 1e-6; X[:, 7] = d["NA22"] * 1e-6
    X[:, 8]  = rmean_um(d["VO11"], d["NA11"], d["SG11"])
    X[:, 9]  = rmean_um(d["VO21"], d["NA21"], d["SG21"])
    X[:, 10] = rmean_um(d["VO12"], d["NA12"], d["SG12"])
    X[:, 11] = rmean_um(d["VO22"], d["NA22"], d["SG22"])
    X[:, 12] = d["HG11"]; X[:, 13] = d["HG21"]; X[:, 14] = d["HG12"]; X[:, 15] = d["HG22"]
    sg2   = np.stack([d["SG11"], d["SG21"]], axis=1)          # SG for the 2 active modes
    fnwrf = np.stack([d["FN11"], d["FN21"]], axis=1)          # emulator applied
    return X, sg2, fnwrf


def reconstruct(files, max_cells=0):
    """Read all wrfout files -> concatenated (X, sg2, fn_wrf) for updraft cells."""
    print("=" * 70, flush=True)
    print("Reconstructing emulator inputs + applied FN from wrfout ...", flush=True)
    Xs, sgs, fns = [], [], []
    for f in files:
        try:
            X, sg2, fnwrf = read_wrfout(f)
        except Exception as e:      # file still being written by a live WRF job
            print(f"  {os.path.basename(f)}: SKIPPED (unreadable: {e})", flush=True)
            continue
        if X is None:
            print(f"  {os.path.basename(f)}: SKIPPED (no time record / no valid updraft cells)", flush=True)
            continue
        print(f"  {os.path.basename(f)}: {len(X):,} updraft cells", flush=True)
        Xs.append(X); sgs.append(sg2); fns.append(fnwrf)
    if not Xs:
        raise RuntimeError("no readable wrfout files with valid updraft cells")
    X   = np.concatenate(Xs, 0)
    sg2 = np.concatenate(sgs, 0)
    fn_wrf = np.concatenate(fns, 0)
    if max_cells and X.shape[0] > max_cells:
        sel = np.random.default_rng(42).choice(X.shape[0], size=max_cells, replace=False)
        sel.sort()
        X, sg2, fn_wrf = X[sel], sg2[sel], fn_wrf[sel]
        print(f"Subsampled to {max_cells:,} cells (--max_cells)", flush=True)
    print(f"Total updraft cells: {X.shape[0]:,}", flush=True)
    return X, sg2, fn_wrf


def run_parcel(X, sg2, lib, workers, indices, deadline=0.0):
    """Solve the parcel model for the given global `indices`. Returns res[k,6]
       columns = (global_idx, fn0, fn1, smax, bad, tout).

       If `deadline` > 0, stop after that many seconds and return only the k cells
       that finished. `indices` must therefore be in random order, so a truncated
       run is still an unbiased sample of the pool (see main())."""
    global G_T, G_P, G_RH, G_W, G_NUM, G_SIG, G_KAP, G_RA
    G_T, G_P, G_RH, G_W = X[:, 0], X[:, 1], X[:, 2], X[:, 3]
    G_NUM = np.ascontiguousarray(X[:, 4:6])     # NA11, NA21  (#/cm3)
    G_SIG = np.ascontiguousarray(sg2)           # SG11, SG21
    G_KAP = np.ascontiguousarray(X[:, 12:14])   # HG11, HG21
    G_RA  = np.ascontiguousarray(X[:, 8:10])    # rmean Aitken, Accum (um)

    n = len(indices)
    out = np.empty((n, 6), dtype=np.float64)
    step = max(1, min(2500, n // 100))
    t0 = time.time()
    done = 0
    stopped_early = False
    with mp.Pool(workers, initializer=_pool_init, initargs=(lib,),
                 maxtasksperchild=5000) as pool:
        for r in pool.imap_unordered(_run_idx, [int(i) for i in indices], chunksize=1):
            out[done] = r; done += 1
            if done % step == 0:
                el = time.time() - t0
                rate = done / max(el, 1e-9)
                print(f"  parcel {done:,}/{n:,} ({100*done//n:3d}%)  {rate:6.1f} samp/s  "
                      f"eta {(n-done)/max(rate,1e-9)/60:5.1f} min", flush=True)
            if deadline and (time.time() - t0) > deadline:
                print(f"  deadline {deadline:.0f}s reached after {done:,}/{n:,} cells "
                      f"-- stopping, keeping completed solves", flush=True)
                stopped_early = True
                pool.terminate()          # kill in-flight workers; drop partial solves
                break
    print(f"Parcel model {'stopped at deadline' if stopped_early else 'done'} "
          f"({time.time()-t0:.0f}s, {done:,} cells solved)", flush=True)
    return out[:done]


def summarize(parcel, wrf):
    print(f"\n===== parcel vs wrfout emulator  (N={len(parcel):,}) =====", flush=True)
    for m, name in enumerate(MODES):
        p, y = parcel[:, m], wrf[:, m]
        ok = np.isfinite(p) & np.isfinite(y)
        p, y = p[ok], y[ok]
        d = y - p
        corr = np.corrcoef(p, y)[0, 1] if len(p) > 2 else np.nan
        print(f"  [{name}]  N={ok.sum():,}")
        print(f"     fn_parcel : mean={p.mean():.4f}  median={np.median(p):.4f}  "
              f">0.99={np.mean(p > 0.99):.3f}  <0.01={np.mean(p < 0.01):.3f}")
        print(f"     fn_wrf(ML): mean={y.mean():.4f}  median={np.median(y):.4f}  "
              f">0.99={np.mean(y > 0.99):.3f}  <0.01={np.mean(y < 0.01):.3f}")
        print(f"     ML-parcel : mean diff={d.mean():+.4f}  MSE={(d**2).mean():.4f}  "
              f"RMSE={np.sqrt((d**2).mean()):.4f}  corr={corr:+.3f}")
        print(f"     fraction parcel<0.5 but ML>0.99 (emulator over-activates): "
              f"{np.mean((p < 0.5) & (y > 0.99)):.3f}", flush=True)


def scatter_panel(ax, parcel, wrf, title):
    ok = np.isfinite(parcel) & np.isfinite(wrf)
    p, y = parcel[ok], wrf[ok]
    if len(p) == 0:
        ax.set_title(title + " (no data)"); return
    ax.hexbin(p, y, gridsize=80, mincnt=1, cmap="viridis",
              extent=[0, 1, 0, 1], norm=mcolors.LogNorm(vmin=1), rasterized=True)
    ax.plot([0, 1], [0, 1], "r--", lw=1)
    mse = ((y - p) ** 2).mean()
    ax.text(0.03, 0.97, f"N={len(p):,}\nMSE={mse:.2e}\nmeanP={p.mean():.2f}\nmeanML={y.mean():.2f}",
            transform=ax.transAxes, va="top", fontsize=8, bbox=dict(fc="white", alpha=0.85))
    ax.set_xlabel("fn_parcel (Fortran box model)")
    ax.set_ylabel("fn_wrf (ML emulator applied)")
    ax.set_title(title, fontsize=10)


def save_fig(fig, path_noext, label):
    """Write both PNG and EPS. Data layers are rasterized (see scatter_panel /
       merge_and_plot) because the PS backend drops alpha; axes and text stay vector."""
    for ext in ("png", "eps"):
        p = f"{path_noext}.{ext}"
        fig.savefig(p, dpi=300, bbox_inches="tight")
        print(f"Wrote {label} -> {p}", flush=True)
    plt.close(fig)


def merge_and_plot(fn_parcel, fn_wrf, good, out_prefix, attempted=None):
    N = len(good)
    na = int(attempted.sum()) if attempted is not None else N
    print(f"Pool cells: {N:,}   attempted: {na:,} ({100.0*na/max(N,1):.1f}% of pool)   "
          f"clean solves: {int(good.sum()):,} ({100.0*int(good.sum())/max(na,1):.1f}% of attempted)",
          flush=True)
    fp, fw = fn_parcel[good], fn_wrf[good]
    summarize(fp, fw)

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for m, name in enumerate(MODES):
        scatter_panel(axes[m], fp[:, m], fw[:, m], f"{name}")
    fig.suptitle("Fortran parcel model vs WRF-applied ML emulator FN "
                 "(y=x => emulator matches physics)", fontsize=12)
    fig.tight_layout()
    print("", flush=True)
    save_fig(fig, out_prefix + "_scatter", "scatter")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    bins = np.linspace(0, 1, 60)
    for m, name in enumerate(MODES):
        ax = axes[m]
        ax.hist(fp[:, m], bins=bins, density=True, alpha=0.6, color="steelblue",
                label="fn_parcel (physics)", rasterized=True)
        ax.hist(fw[:, m], bins=bins, density=True, alpha=0.6, color="tomato",
                label="fn_wrf (ML applied)", rasterized=True)
        ax.set_yscale("log")
        ax.set_xlabel("activation fraction"); ax.set_ylabel("density (log)")
        ax.set_title(name, fontsize=10); ax.legend(fontsize=8)
    fig.suptitle("Distribution: Fortran parcel model vs WRF-applied ML emulator", fontsize=12)
    fig.tight_layout()
    save_fig(fig, out_prefix + "_hist", "histogram")

    print("\n" + "=" * 70)
    print("fn_parcel = ground-truth Fortran box model on the WRF states;")
    print("fn_wrf    = ML emulator fraction WRF actually applied (FN11/FN21).")
    print("Large ML>0.99 while parcel<0.5 => the deployed emulator over-activates")
    print("on the real WRF distribution (the online saturation).")
    print("=" * 70, flush=True)


def res_to_arrays(res, N):
    """Scatter a res[?,6] (global_idx,fn0,fn1,smax,bad,tout) into aligned [N] arrays.
       `res` may cover only a subset of N (subsample / --deadline cut)."""
    fn_parcel = np.full((N, 2), np.nan)
    smax = np.full(N, np.nan); bad = np.ones(N, np.int32); tout = np.zeros(N, np.int32)
    attempted = np.zeros(N, bool)
    idx = res[:, 0].astype(np.int64)
    attempted[idx] = True
    fn_parcel[idx, 0] = np.clip(res[:, 1], 0, 1)
    fn_parcel[idx, 1] = np.clip(res[:, 2], 0, 1)
    smax[idx] = res[:, 3]; bad[idx] = res[:, 4].astype(np.int32); tout[idx] = res[:, 5].astype(np.int32)
    good = (bad == 0) & (tout == 0) & (smax > 0)
    return fn_parcel, good, attempted


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--phase", choices=["all", "reconstruct", "parcel", "merge"], default="all")
    ap.add_argument("--files", nargs="*", default=None, help="wrfout files")
    ap.add_argument("--all", action="store_true", help="use all wrfout_d01_*")
    ap.add_argument("--lib", default=DEFAULT_LIB, help="path to libactivate.so")
    ap.add_argument("--workers", type=int, default=0,
                    help="parcel workers (0 = all allocatable cores)")
    ap.add_argument("--max_cells", type=int, default=0,
                    help="0 = ALL updraft cells; >0 randomly subsamples")
    ap.add_argument("--inputs_npz", default=os.path.join(HERE, "parcel_inputs.npz"),
                    help="reconstruct writes / parcel+merge read X,sg2,fn_wrf here")
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--subsample", type=int, default=0,
                    help="parcel phase: solve a fixed seeded random subset of this "
                         "many cells (0 = all); same subset across all shards")
    ap.add_argument("--deadline", type=float, default=0.0,
                    help="parcel phase: seconds to solve before stopping and keeping "
                         "whatever finished (0 = no limit). Set below the job walltime "
                         "so every shard always writes its npz and merge can run.")
    ap.add_argument("--shard_npz", default=None, help="parcel phase output for this shard")
    ap.add_argument("--shard_glob", default=os.path.join(HERE, "parcel_shard_*.npz"),
                    help="merge phase: glob of shard npz files")
    ap.add_argument("--out_prefix", default=os.path.join(HERE, "parcel_vs_wrfout"))
    args = ap.parse_args()
    if not os.path.exists(args.lib) and args.phase in ("all", "parcel"):
        raise FileNotFoundError(args.lib)

    ncores  = len(os.sched_getaffinity(0))
    workers = args.workers if args.workers > 1 else ncores

    # ---- resolve file list (reconstruct / all) ----
    if args.phase in ("all", "reconstruct"):
        if args.files:
            files = args.files
        else:
            allf = sorted(glob.glob(os.path.join(TEST_DIR, "wrfout_d01_*")))
            if not allf:
                raise FileNotFoundError(f"no wrfout_d01_* under {TEST_DIR}")
            files = allf if args.all else allf[:1]

    # ==================== PHASE: reconstruct ====================
    if args.phase == "reconstruct":
        X, sg2, fn_wrf = reconstruct(files, args.max_cells)
        np.savez(args.inputs_npz, X=X, sg2=sg2, fn_wrf=fn_wrf)
        print(f"Wrote inputs -> {args.inputs_npz}  (N={X.shape[0]:,})", flush=True)
        return

    # ==================== PHASE: parcel (one shard) ====================
    if args.phase == "parcel":
        z = np.load(args.inputs_npz)
        X, sg2 = z["X"], z["sg2"]
        N = X.shape[0]
        # Seeded permutation: identical on every shard, so the strided slices stay
        # disjoint. Left UNSORTED so that a --deadline cut keeps a random prefix
        # (sorting would bias a truncated run toward the earliest wrfout times).
        pool_idx = np.random.default_rng(42).permutation(N)
        if args.subsample and args.subsample < N:
            pool_idx = pool_idx[:args.subsample]
        idx = pool_idx[args.shard::args.nshards]         # strided, non-overlapping
        out_npz = args.shard_npz or os.path.join(HERE, f"parcel_shard_{args.shard}.npz")
        print(f"[shard {args.shard}/{args.nshards}] {len(idx):,} of "
              f"{len(pool_idx):,} selected ({N:,} total) cells "
              f"on {workers} workers ({ncores} cores), deadline={args.deadline or 0:.0f}s",
              flush=True)
        res = run_parcel(X, sg2, args.lib, workers, idx, args.deadline)
        np.savez(out_npz, res=res)
        print(f"[shard {args.shard}] wrote {out_npz}", flush=True)
        return

    # ==================== PHASE: merge ====================
    if args.phase == "merge":
        z = np.load(args.inputs_npz)
        fn_wrf = z["fn_wrf"]; N = fn_wrf.shape[0]
        parts = sorted(glob.glob(args.shard_glob))
        if not parts:
            raise FileNotFoundError(args.shard_glob)
        res = np.concatenate([np.load(p)["res"] for p in parts], axis=0)
        print(f"Merged {len(parts)} shards -> {res.shape[0]:,} solved cells (N={N:,})", flush=True)
        fn_parcel, good, attempted = res_to_arrays(res, N)
        merge_and_plot(fn_parcel, fn_wrf, good, args.out_prefix, attempted)
        return

    # ==================== PHASE: all (single node) ====================
    X, sg2, fn_wrf = reconstruct(files, args.max_cells)
    N = X.shape[0]
    print(f"\nRunning Fortran parcel model on {N:,} cells "
          f"({workers} workers, {ncores} cores allocatable) ...", flush=True)
    idx = np.random.default_rng(42).permutation(N)   # random order: --deadline-safe
    res = run_parcel(X, sg2, args.lib, workers, idx, args.deadline)
    fn_parcel, good, attempted = res_to_arrays(res, N)
    merge_and_plot(fn_parcel, fn_wrf, good, args.out_prefix, attempted)


if __name__ == "__main__":
    main()
