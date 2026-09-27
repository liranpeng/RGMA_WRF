# parcel_vs_ARG_vs_ML_ALLCASES

Self-contained bundle to recreate `parcel_vs_ARG_vs_ML_ALLCASES.{png,eps}` and
`parcel_vs_ARG_vs_ML_ALLCASES_stats.txt`: the 2x3 pooled hexbin figure comparing
Fortran parcel-model activation fraction (ground truth), the ML emulator fraction
WRF applied (FN11/FN21), and the ARG fraction WRF computed for the same cell
(FN13/FN23), for all four 2017-07-14 test cases (1aer/3aer x org/mid).

Copied 2026-09-08 from
`/scratch/07088/tg863871/Perlm_Backup/WRF_stam3_ml_mixout_bce_lr1e-4_cosine_more_fix/test/`
(the only tree that holds the input data; the copies of the script under
`key/Online_test/code` and `Model_Backup/WRF_ML/test` are byte-identical).

## Recreate

    ./run.sh                        # all four cases
    ./run.sh --cases 20170714_1aer_org 20170714_3aer_mid   # subset

Outputs land in this directory. Uses `~/.conda/envs/liranenv_gpu/bin/python`
(falls back to `liranenv`); needs numpy 2 + matplotlib. Runs on a login node in a
few minutes; no SLURM needed.

## Layout

    plot_allcases_arg.py        the figure script (unchanged; finds cases relative to itself)
    run.sh                      wrapper that picks the conda python
    <case>/online/
        parcel_inputs.npz       X / sg2 / fn_wrf / fn_arg, 1.5M cells per case (fn_arg added by add_fn_arg.py)
        parcel_shard_0.npz      Fortran parcel-model results (idx, fn0, fn1, smax, bad, timeout)
        diag_fn_ml_offline.npz  offline emulator answer, used to flag "lost output" cells
    original_figure/            png / eps / stats as produced 2026-08-04 (reference)
    provenance/                 scripts that produced the inputs (not needed to replot):
        add_fn_arg.py               appended fn_arg (FN13/FN23) to parcel_inputs.npz
        compare_parcel_vs_wrfout.py + submit_parcel_compare_stam3.slurm   parcel-model run -> parcel_shard_*.npz
        parcel_cmp_*.out            logs of those runs
        setup_online_dirs.sh, run_all_analysis.sh, plot_allcases_lost_removed.py   related pipeline

## Sample cleaning (see docstring in plot_allcases_arg.py)

Rows kept = good parcel solve (bad==0, no timeout, smax>0)
          & not "lost" (fn_wrf==0 in both modes while offline emulator > 0.01)
          & fn_arg within [0,1] in both modes.

## v2 (2026-09-12): paper-style variant

    plot_allcases_arg_v2.py     same data/cleaning; drops the suptitle and per-panel titles,
                                panel letters centred above each panel, axis labels typeset as
                                FN_parcel / FN_ML / FN_ARG (subscripts), labels 28 pt, ticks 24 pt
    -> parcel_vs_ARG_vs_ML_ALLCASES_v2.{png,eps,_stats.txt}   (run with the same python as run.sh)
