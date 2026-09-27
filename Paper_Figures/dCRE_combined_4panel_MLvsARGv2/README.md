# dCRE_combined_4panel_MLvsARGv2

Self-contained bundle to recreate `dCRE_combined_4panel_MLvsARGv2.{png,pdf,eps}`:
the 2x2 figure of the pixel-wise shortwave cloud radiative effect change (dCRE)
decomposition, ML-ACT vs ARG-ACT (WRF_dm_v2), for the org (top) and mid (bottom)
backgrounds. Left column: half-hourly stacked components, ML solid | ARG hatched.
Right column: time-mean bars of total, CF, A_Nc, A_LWP, A_re, A_cov.

Copied 2026-09-08 from `/scratch/07088/tg863871/WRF_liran/JP/` (script) and its
output dir `dcre_figures_v23_argv2_1300/` (CSVs + figure, produced 2026-08-16).

## Recreate

    ./run.sh                              # window end 2017-07-15 13:00 (published)
    ./run.sh --tend ""                    # whole CSV record
    ./run.sh --stem my_name               # different output stem

`run.sh` runs the script with `--outdir dcre_csv` and moves the three outputs
to this directory. Pure replot: reads only the CSVs, never a wrfout file.
Uses `~/.conda/envs/liranenv/bin/python` (falls back to `liranenv_gpu`);
needs numpy, pandas, matplotlib. Seconds on a login node.

## Layout

    plot_dcre_combined_4panel_MLvsARGv2.py   the figure script (unchanged)
    run.sh                                   wrapper
    dcre_csv/                                inputs (all CSVs from the producer run)
        dcre_pixelwise_{org_3aer_org, mid_3aer_mid,
                        ARG_org_3aer_ARG_org, ARG_mid_3aer_ARG_mid}_domain.csv
                                             <- the four files this figure reads
        dcre_pixelwise_*_by_box.csv, dcre_meanbased_*.csv, dispersion_*.csv
                                             other products of the same run (not used here)
    original_figure/                         png / pdf / eps as produced 2026-08-16
    provenance/                              not needed to replot:
        trace_east_boundary_subdomains_lwp_v23_argv2.py   wrote the CSVs from wrfout
        run_dcre_v23_argv2_1300.sh           driver: WRF_TEND=13:00, WRF_ARG_DIR=WRF_dm_v2
        dcre_v23_argv2_nohup_20260815_230208.log          log of that run
        plot_dcre_combined_4panel.py, plot_dcre_combined_4panel_MLvsARG.py   earlier versions

## Source model runs (for the CSVs; NOT included, multi-TB)

    ML-ACT : /scratch/07088/tg863871/Perlm_Backup/WRF_stam3_ml_mixout_bce_lr1e-4_cosine_more_fix/test/20170714_{1,3}aer_{org,mid}
    ARG-ACT: /scratch/07088/tg863871/WRF_dm_v2/test/20170714_{org,mid}_{1,3}aer_ARG

Pairs are 3aer minus 1aer within each scheme and background. WRF_dm_v2 (unlike
WRF_dm) populates TAU_QC_TOT/RE_QC, so the ARG Nc/LWP/re/cov terms are real
measurements; the script detects an analytic-fallback record from R_tau and
greys those bars automatically.

## v3 (2026-09-16): larger fonts, letter-only panel titles

`dCRE_combined_4panel_MLvsARGv3.{png,pdf,eps}` from
`plot_dcre_combined_4panel_MLvsARGv3.py`: same data, layout and window as v2,
with every font 1.35x larger, panel titles reduced to (a)-(d) (the descriptive
text now belongs in the caption), and extra y-range headroom so the larger
legend no longer overlaps the bars. The ML/ARG time-mean summary box in (b),(d) is dropped, and the
number term is labelled N_d (not Nc) to match the consequences figure; the CSV
column is still DCRE_A_Nc.

    ./run.sh   # still builds v2; for v3:
    ~/.conda/envs/liranenv/bin/python plot_dcre_combined_4panel_MLvsARGv3.py --outdir dcre_csv \
        && mv -f dcre_csv/dCRE_combined_4panel_MLvsARGv3.{png,pdf,eps} .

## v4 (2026-09-16): + mid2 row and N_d / LWP / CF composites

`dCRE_combined_4panel_MLvsARGv4.{png,pdf,eps}` from
`plot_dcre_combined_MLvsARGv4.py` (`./run_v4.sh`).  v3 plus:

* row 3 = mid2 background, (e) stacks and (f) bars, same style;
* bottom row (g)-(i): the (a)/(c)/(e) stacks re-binned in cloud STATE instead
  of time.  Sampling unit = one 30-min window of one pair (the same windows
  as the bar groups above); the 30-min means of org, mid and mid2 are pooled
  (24 samples per scheme) and sorted into fixed bins of the pair-mean
  (1aer+3aer)/2 state: N_d (in-cloud QNDROP below 2 km, cm-3), LWP (LWP_TOT
  over cloudy columns), CF (cloudy-column fraction).  Bin edges are the
  STATE_SPECS constants at the top of the script.  ML | ARG side by side,
  total +-SE over the samples in the bin, n=ML|ARG below.

Extra inputs added to the bundle:
    dcre_csv/dcre_*mid2*.csv            from WRF_liran/JP/dcre_figures_v24_mid2_1300
                                        (its org/mid CSVs are byte-identical to v23's)
    state_cache/arg_vs_ml_consequences_rwp_v11_series_CCN2.npz
                                        per-case 2-min N_d / LWP / CF series for all
                                        12 runs (from npzfile/, made by the
                                        consequences-figure script v11)
    provenance/trace_east_boundary_subdomains_lwp_v24_mid2.py, run_dcre_v24_mid2_1300.sh,
               dcre_v24_mid2_nohup_20260824_213549.log
Tables written next to the figure:
    dCRE_combined_4panel_MLvsARGv4_state_merged_2min.csv    all times, dCRE + state, 6 pairs
    dCRE_combined_4panel_MLvsARGv4_state_samples_30min.csv  the 48 pooled samples
    dCRE_combined_4panel_MLvsARGv4_state_composites.csv     what (g)-(i) draw

## v5 (2026-09-16): v4 without mid2

`dCRE_combined_4panel_MLvsARGv5.{png,pdf,eps}` (`./run_v5.sh`): the same
script with `--backgrounds org,mid`.  Two background rows (a)-(d) and the
composites (e)-(g) pooled over org and mid only (16 samples per scheme).
State bins that hold no sample are left off the axis (so LWP shows <75 and
>=150 only, CF shows <0.50, 0.50-0.80, 0.95-1.00).  Its tables are
`dCRE_combined_4panel_MLvsARGv5_state_*.csv`.

## v6 (2026-09-16): org/mid rows, N_d and LWP composites over org+mid+mid2

`dCRE_combined_4panel_MLvsARGv6.{png,pdf,eps}` (`./run_v6.sh`): the same
script with `--backgrounds org,mid --composite-backgrounds org,mid,mid2
--composites Nd,LWP --stack-legend first`.  Rows (a)-(d) show org and mid
only; (e) N_d and (f) LWP composites are pooled over all three backgrounds
(24 samples per scheme), so the mid2 pair enters the figure only there.  The
cloud-fraction composite is not drawn, and the stacked-component legend
appears only in (a) (it is the same in every stacked panel).  Tables:
`dCRE_combined_4panel_MLvsARGv6_state_*.csv` (the composites table covers
only the two drawn composites).

Options added to the v4 script for this: `--composites` (subset/order of
Nd,LWP,CF) and `--stack-legend all|first|top|none`.

## Notation change (2026-09-17), all of v3-v6

Sub-term labels now follow the manuscript: ΔCRE_CF, ΔCRE_{N_d}, ΔCRE_LWP,
ΔCRE_{r_e}, ΔCRE_cov (the "A," prefix is dropped from the figure labels).
CSV column names are unchanged (DCRE_A_Nc, DCRE_A_LWP, DCRE_A_re, DCRE_A_cov).
v3-v6 png/pdf/eps were all re-rendered.

## 2026-09-17: no n labels, larger legends (all of v3-v6)

The per-bin "n=" counts are no longer drawn in any panel (the counts are in
the *_state_composites.csv tables and, for the time bins, are 13 for 09:00 and
15 for every later half hour).  Legend fonts: component legend 23 pt, bar-panel
legend 24 pt.  In v6 the component legend is one strip above the whole figure
(`--stack-legend top`); v3/v4/v5 keep in-panel legends with more y headroom.

## dCRE_ARGminusML_30min (2026-09-20): when does each channel favour ARG?

`plot_dcre_argminusml_30min.py` -> `dCRE_ARGminusML_30min.{png,pdf}` (+ `_table.csv`).
Top row: ARG-ACT minus ML-ACT per 30-min bin for every dCRE term and the
total (+-1 SE, variances of the two schemes added), one column per
background (org, mid, mid2).  Dotted verticals + corner notes: first bin from
which the LWP term / the total is negative by more than 1 SE and stays so.
Bottom row: ARG/ML ratio of the control (solid) and 3aer (dashed) cloud state
from state_cache/ (in-cloud N_d, LWP, RWP), log axis; verticals where the
control LWP ratio drops through 1 and the RWP ratio rises through 1.
`--backgrounds org,mid+mid2` pools the two moist pairs (2-min records
concatenated per bin, state ratios averaged) -> `dCRE_ARGminusML_30min_pooled.*`.
