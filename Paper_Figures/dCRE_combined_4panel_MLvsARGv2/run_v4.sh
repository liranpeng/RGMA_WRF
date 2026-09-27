#!/bin/bash
# Recreate dCRE_combined_4panel_MLvsARGv4.{png,pdf,eps} and its three *_state_*.csv
# tables in this directory: 3 background rows (org, mid, mid2) + N_d / LWP / CF
# composites.  Reads dcre_csv/ and state_cache/ only; no wrfout access.
cd "$(dirname "$0")"
PY=$HOME/.conda/envs/liranenv/bin/python
[ -x "$PY" ] || PY=$HOME/.conda/envs/liranenv_gpu/bin/python
"$PY" plot_dcre_combined_MLvsARGv4.py --outdir dcre_csv --figdir . "$@"
