#!/bin/bash
# Recreate arg_relerr_map_filtered.{png,eps} + _stats.txt in this directory.
# Defaults reproduce the published figure: --tau 0.1 --nbins 24 --min_cnt 20 --min_cnt_ait 10 --fn_min 0.01
cd "$(dirname "$0")"
PY=$HOME/.conda/envs/liranenv_gpu/bin/python
[ -x "$PY" ] || PY=$HOME/.conda/envs/liranenv/bin/python
"$PY" plot_arg_relerr_map.py "$@"
