#!/bin/bash
# Recreate parcel_vs_ARG_vs_ML_ALLCASES.{png,eps} + _stats.txt in this directory.
# Needs numpy + matplotlib (numpy 2 is fine; the npz caches were written with it).
cd "$(dirname "$0")"
PY=$HOME/.conda/envs/liranenv_gpu/bin/python
[ -x "$PY" ] || PY=$HOME/.conda/envs/liranenv/bin/python
"$PY" plot_allcases_arg.py "$@"
