# arg_relerr_map_filtered

Self-contained bundle to recreate `arg_relerr_map_filtered.{png,eps}` and
`arg_relerr_map_filtered_stats.txt`: the 2x4 map of the ARG scheme's MEDIAN
RELATIVE error in activated fraction, (fn_ARG - fn_parcel)/fn_parcel, against
the Fortran parcel model, binned in (wbar, na_accum), (wbar, na_Aitken),
(na_accum, na_Aitken) and (wbar, pressure). Rows = Aitken (a-d), Accum (e-h).
Black contours = failure rate |fn_ARG - fn_parcel| > 0.1.

Copied 2026-09-08 from
`/scratch/07088/tg863871/Perlm_Backup/WRF_stam3_ml_mixout_bce_lr1e-4_cosine_more_fix/test/`
(the 2026-08-06 22:47 version, the newest). `key/ARG_relerr_map/` holds a
byte-identical copy of the scripts and data but with a `data/` layout that the
script cannot run from; this bundle puts the case dirs beside the script so it
runs in place.

## Recreate

    ./run.sh                     # defaults = published settings
    ./run.sh --tau 0.1 --nbins 24 --min_cnt 20 --min_cnt_ait 10 --fn_min 0.01
    ./run.sh --cases 20170714_1aer_org 20170714_3aer_mid     # subset
    ./run.sh --callsite old_cloud   # needs 'callsite' in the npz (not present here)

Outputs land in this directory. Uses `~/.conda/envs/liranenv_gpu/bin/python`
(falls back to `liranenv`); needs numpy 2 + matplotlib. Runs on a login node in
well under a minute; no SLURM needed.

## Layout

    plot_arg_relerr_map.py      the figure script (unchanged)
    plot_allcases_arg.py        HARD DEPENDENCY: provides CASES, MODES, load_case (sample cleaning)
    run.sh                      wrapper that picks the conda python
    <case>/online/
        parcel_inputs.npz       X (state), sg2, fn_wrf, fn_arg; 1.5M cells per case
        parcel_shard_0.npz      Fortran parcel-model results (idx, fn0, fn1, smax, bad, timeout)
        diag_fn_ml_offline.npz  offline emulator answer, used to flag "lost output" cells
    original_figure/            png / eps / stats as produced 2026-08-06 (reference)
    provenance/                 not needed to replot:
        README_key_ARG_relerr_map.txt   the original write-up with the headline numbers
        add_fn_arg.py, add_fn_arg.log   appended fn_arg (FN13/FN23) to parcel_inputs.npz
        build_fn_ml_cache.py, model/mlp_wrapper.pt   built diag_fn_ml_offline.npz
        compare_parcel_vs_wrfout.py, submit_parcel_compare_stam3.slurm, parcel_cmp_*.out
                                        parcel-model run -> parcel_shard_*.npz
        plot_arg_failure_map_bwr.py     companion ABSOLUTE-difference map
        plot_arg_failure_map.py, plot_arg_density_map.py, setup_online_dirs.sh,
        run_all_analysis.sh, relerr_map.log

The four `<case>/online/` inputs are byte-identical to those in the sibling
bundle `../parcel_vs_ARG_vs_ML_ALLCASES/`; both figures use the same 56,689
kept samples. Duplicated so each folder stands alone.

## X column index (parcel_inputs.npz "X")

    0 T (K)   1 pressure (hPa)   3 wbar (cm/s)   4 na_Aitken (cm-3)   5 na_accum (cm-3)

## v2 (2026-09-12): paper-style variant

    plot_arg_relerr_map_v2.py   same data/settings; drops the suptitle and per-panel
                                titles, axis labels 14->28 pt, tick labels 12->24 pt
    -> arg_relerr_v2_map_filtered.{png,eps,_stats.txt}   (run with the same python as run.sh)

## v3 (2026-09-15): v2 split into two 2x2 figures

    plot_arg_relerr_map_v3.py   same data/settings/cleaning as v2; one figure per mode,
                                2x2 panels each lettered (a)-(d), NO per-panel stats box
    -> arg_relerr_v3_aitken_map_filtered.{png,eps}   = v2 panels (a)-(d)
       arg_relerr_v3_accum_map_filtered.{png,eps}    = v2 panels (e)-(h)
       arg_relerr_v3_map_filtered_stats.txt          (identical to the v2 stats)
    --letters continue   keeps v2's (e)-(h) on the Accum figure
    --box_fs 20          restores the per-panel stats box at that font size (stats are in the _stats.txt)
