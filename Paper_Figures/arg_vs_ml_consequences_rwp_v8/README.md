# arg_vs_ml_consequences_rwp_v8

Self-contained bundle to recreate `arg_vs_ml_consequences_rwp_v8.{png,eps}`:
the 5x2 time-series figure of the consequences of the ARG activation bias,
ML-ACT (solid) vs ARG-ACT (dashed, WRF_dm_v2) for org / org_3aer / mid /
mid_3aer, 09:00-13:00 UTC 15 July 2017. Panels: (a) N_acc, (b) N_ait,
(c) CCN at S=0.05%, (d) N_d, (e) r_v, (f) LWP, (g,i) rain water path and
(h,j) cloudy-column fraction, each split mid / org.

Copied 2026-09-08 from `/scratch/07088/tg863871/figures/` (script + figure,
produced 2026-08-21 23:22) and `/scratch/07088/tg863871/npzfile/` (series cache).

## Recreate

    ./run.sh                    # = python plot_arg_vs_ml_consequences_rwp_v8.py \
                                #     --out arg_vs_ml_consequences_rwp_v8 --font-scale 1.6 --tend 13:00
    ./run.sh --ccn CCN3         # other supersaturation (needs a CCN3 cache; none shipped)

Outputs and `run.log` (which includes the ML-vs-ARG summary table printed by
the script) land in this directory. Uses `~/.conda/envs/liranenv/bin/python`
(falls back to `liranenv_gpu`); needs numpy, netCDF4, matplotlib. Seconds on a
login node: nothing is read from wrfout because every case comes from the cache.

## How the data flow works

The script keeps per-case, per-output-time domain means in a cache
`<out>_series_<CCN>.npz`. The original v8 run had no v8 cache, so it SEEDED
from `arg_vs_ml_consequences_rwp_v4_series_CCN2.npz` (default `--seed-cache`):
the four ML cases under their own keys, and the four WRF_dm_v2 ARG cases
remapped from that cache's `V2_*` keys onto `ARG_*`. The v4 cache's own `ARG_*`
keys (old WRF_dm tree) are deliberately ignored. That cache is shipped here
unchanged (mtime 2026-08-18), and `run.sh` reproduces exactly that seeding.

Re-reading from wrfout (`--refresh all` or a missing case) needs the model
trees below, ~14 GB per wrfout, and hours; do it on a compute node.

    ML-ACT : /scratch/07088/tg863871/Perlm_Backup/WRF_stam3_ml_mixout_bce_lr1e-4_cosine_more_fix/test/20170714_{1,3}aer_{org,mid}
    ARG-ACT: /scratch/07088/tg863871/WRF_dm_v2/test/20170714_{org,mid}_{1,3}aer_ARG

## Layout

    plot_arg_vs_ml_consequences_rwp_v8.py       the figure script (unchanged)
    arg_vs_ml_consequences_rwp_v4_series_CCN2.npz   the data: 12 series x 8 cases (+4 old-WRF_dm keys, unused)
    run.sh                                      wrapper with the published options
    original_figure/                            png / eps produced 2026-08-21, plus the LaTeX caption
    provenance/
        plot_arg_vs_ml_consequences_rwp_v4.py   the version that wrote the cache from wrfout
        plot_arg_vs_ml_consequences_rwp_v7.py   immediate predecessor (v8 = v7 with RWP/CLDCOL panels)
        submit_plot_consequences_v6.sh          SLURM script used for the v6 wrfout read
        consequences_rwp_v4.log, conseq_v6_login.log   logs of those runs

## Reproduction check (2026-09-08)

`./run.sh` with liranenv (matplotlib 3.10.3) regenerated a PNG that is
pixel-identical to `original_figure/arg_vs_ml_consequences_rwp_v8.png`
(2893 x 4249, zero differing pixels). `--tend 13:00` is essential: without it
the curves run to each case's end (13:04-14:22) and the layout changes.

## v12 (2026-09-15): 8-panel variant

    plot_arg_vs_ml_consequences_rwp_v12.py   v8 without the Aitken-number panel (old b) and the
                                             r_v panel (old e); the remaining eight relettered
                                             (a) N_acc, (b) CCN, (c) N_d, (d) LWP, (e,f) RWP and
                                             cloud fraction mid, (g,h) the same for org. 4x2 grid.
    -> arg_vs_ml_consequences_rwp_v12.{png,eps}, run_v12.log
    Same data, seeding and flags as run.sh:
        python plot_arg_vs_ml_consequences_rwp_v12.py --out arg_vs_ml_consequences_rwp_v12 \
               --font-scale 1.6 --tend 13:00
    The ML-vs-ARG summary table in run_v12.log is identical to run.log (r_v is still
    computed and tabulated, only no longer plotted).
