================================================================================
ARG relative-error map against the Fortran parcel model
Collected 2026-08-04 from
  /scratch/07088/tg863871/Perlm_Backup/WRF_stam3_ml_mixout_bce_lr1e-4_cosine_more_fix/test
================================================================================

WHAT THIS IS
    figures/arg_relerr_map_filtered.png  (+ .eps, + _stats.txt)

    colour        = MEDIAN RELATIVE error of the ARG activated fraction
                        delta_rel = (fn_ARG - fn_parcel) / fn_parcel
    contours      = failure rate |fn_ARG - fn_parcel| > 0.1, at 0.1/0.25/0.5/0.75/0.9
    rows          = Aitken (a-d), Accumulation (e-h)
    columns       = na_acc vs wbar | na_ait vs wbar | na_ait vs na_acc | pres vs wbar

    Relative error is undefined where the parcel model does not activate, so
    bins use only samples with fn_parcel > 0.01; bins with too few such samples
    are left grey.  The occupancy threshold is per row: 20 for Accumulation,
    10 for Aitken (--min_cnt / --min_cnt_ait).  Only 5.7% of the pool is usable
    in the Aitken mode, so the Accum threshold blanked nearly the whole top row;
    at 10 the Aitken panels show 67/80/58/89 bins instead of 47/43/11/54.  Bins
    resting on 10-19 samples are noisier -- read them as pattern, not value.
    The per-bin statistic is the median (delta_rel is
    bounded below by -1 but unbounded above, so its mean is dragged up by a
    thin tail of over-activations); the stats file reports both.

RESULT, all four cases pooled (N = 56,689 kept samples)
    Accum:  usable 54,775 (96.6% of pool)
            median rel = -0.3611   mean rel = -0.3975   frac(ARG<parcel) = 0.9994
            mean fn_parcel 0.2112  mean fn_ARG 0.1449
    Aitken: usable  3,241 ( 5.7% of pool)
            median rel = -0.4384   mean rel = -0.5085   frac(ARG<parcel) = 1.0000
            mean fn_parcel 0.0100  mean fn_ARG 0.0071

    !! The manuscript currently quotes a pooled accumulation-mode median
    relative error of -0.503.  With all four cases included (3aer_mid was
    added 2026-08-04, see below) that number is now -0.361.  The -0.503 value
    is close to the CURRENT Aitken-mode MEAN (-0.5085), so check which
    statistic the text intends before updating it.

    The deficit deepens sharply with N/w (accumulation mode):
      N/w < 1e3        median rel -0.378   (n = 19,225)
      1e3 - 1e4        -0.349 / -0.317
      1e4 - 3e4        -0.364
      3e4 - 1e5        -0.615
      1e5 - 3e5        -0.878
      3e5 - 1e6        -0.969        (n =  1,272)
    i.e. in the kinetically limited regime (N/w > 1e4, Ghosh et al. 2025) ARG
    activates almost nothing the parcel model activates.

    Named regimes (accumulation mode):
      all cells                                n=54,775  median -0.361  Nact ARG/parcel 0.630
      Phinney regime  V<50 cm/s & Na_acc>500   n=26,329  median -0.349  Nact ARG/parcel 0.612
      Phinney 'good'  V>50      & Na_acc<500   n= 2,778  median -0.321  Nact ARG/parcel 0.698

    CAVEAT carried on the figure: fn_ARG (FN13/FN23) is a persistent WRF array
    written only by the OLD_CLOUD call site; on cells last served by the
    GROW_SHRINK call site it is stale, so part of this error is time-mismatch
    rather than scheme error.

LAYOUT
    code/       plot_arg_relerr_map.py      makes this figure
                plot_allcases_arg.py        HARD DEPENDENCY -- plot_arg_relerr_map.py
                                            does "from plot_allcases_arg import
                                            CASES, MODES, load_case"
                add_fn_arg.py               adds fn_arg to parcel_inputs.npz
                build_fn_ml_cache.py        builds diag_fn_ml_offline.npz
                plot_arg_failure_map_bwr.py the companion ABSOLUTE-difference map
                plot_arg_failure_map.py, plot_arg_density_map.py,
                plot_allcases_lost_removed.py
    model/      mlp_wrapper.pt -- deployed TorchScript emulator (mtime 2026-07-31)
    data/<case>/  parcel_inputs.npz (X, sg2, fn_wrf, fn_arg; 1,500,000 rows),
                parcel_shard_0.npz, diag_fn_ml_offline.npz, the per-case
                driver script + slurm file, and the run logs
    figures/    arg_relerr_map_filtered.{png,eps} + _stats.txt
    logs/       add_fn_arg.log

    data/ here is a byte-identical copy of key/Online_test/data -- both figures
    are built from the same inputs.  Kept duplicated so each subfolder stands
    alone; if one is ever regenerated, refresh the other.

PROVENANCE
    Same input chain as key/Online_test (see that README for full detail):
      * parcel jobs 3368245/3368247/3368246 finished 2026-08-03 20:11-20:18,
        3368822 (3aer_mid) finished 2026-08-04 03:34
      * fn_arg attached 2026-08-04 by add_fn_arg.py, row alignment verified
        exactly against the stored fn_wrf in all four cases (1,500,000 rows)
      * diag_fn_ml_offline.npz built 2026-08-04 from ml_models/mlp_wrapper.pt
    This figure is therefore the first version that includes 3aer_mid.
    The four cases are NOT over identical model-time windows: 3aer_mid's
    inputs were reconstructed from a longer wrfout record than the other three.

HOW TO REGENERATE
    Python env: $HOME/.conda/envs/liranenv_gpu

      python plot_arg_relerr_map.py                 # defaults used here
      python plot_arg_relerr_map.py --tau 0.1 --nbins 24 --min_cnt 20 \
                                    --min_cnt_ait 10 --fn_min 0.01

    Options: --cases (subset), --callsite {all,old_cloud,grow_shrink},
    --out_prefix.  Writes <out_prefix>_map_filtered.{png,eps} and
    <out_prefix>_map_filtered_stats.txt next to the script; the default prefix
    is "arg_relerr" in the script's own directory.

    Paths resolve from the script's location, so running from this directory
    will NOT find the case data -- run in the original test tree, or edit
    HERE at the top of plot_allcases_arg.py.
================================================================================
