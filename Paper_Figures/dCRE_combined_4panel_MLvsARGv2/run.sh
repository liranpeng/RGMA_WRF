#!/bin/bash
# Recreate dCRE_combined_4panel_MLvsARGv2.{png,pdf,eps} in this directory from the
# CSVs in dcre_csv/ (pure replot; no wrfout access). Window end defaults to 2017-07-15 13:00.
cd "$(dirname "$0")"
PY=$HOME/.conda/envs/liranenv/bin/python
[ -x "$PY" ] || PY=$HOME/.conda/envs/liranenv_gpu/bin/python
"$PY" plot_dcre_combined_4panel_MLvsARGv2.py --outdir dcre_csv "$@" && mv -f dcre_csv/dCRE_combined_4panel_MLvsARGv2.{png,pdf,eps} .
