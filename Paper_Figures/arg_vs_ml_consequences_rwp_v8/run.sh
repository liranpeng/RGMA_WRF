#!/bin/bash
# Recreate arg_vs_ml_consequences_rwp_v8.{png,eps} in this directory.
# All eight cases are read from the v4 series cache shipped here (seeded exactly as the
# original run was), so no wrfout files are touched. Full run log -> run.log
cd "$(dirname "$0")"
PY=$HOME/.conda/envs/liranenv/bin/python
[ -x "$PY" ] || PY=$HOME/.conda/envs/liranenv_gpu/bin/python
"$PY" plot_arg_vs_ml_consequences_rwp_v8.py --out arg_vs_ml_consequences_rwp_v8 --font-scale 1.6 --tend 13:00 "$@" 2>&1 | tee run.log
