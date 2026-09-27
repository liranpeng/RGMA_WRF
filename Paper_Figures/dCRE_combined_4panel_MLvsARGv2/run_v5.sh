#!/bin/bash
# v5 = v4 without mid2: org and mid rows only, composites pooled over org+mid.
cd "$(dirname "$0")"
PY=$HOME/.conda/envs/liranenv/bin/python
[ -x "$PY" ] || PY=$HOME/.conda/envs/liranenv_gpu/bin/python
"$PY" plot_dcre_combined_MLvsARGv4.py --outdir dcre_csv --figdir . \
      --backgrounds org,mid --stem dCRE_combined_4panel_MLvsARGv5 "$@"
