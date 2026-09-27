#!/bin/bash
# v6: org and mid rows only; composites (e) N_d and (f) LWP pooled over org, mid AND
# mid2; stacked-component legend only in (a).
cd "$(dirname "$0")"
PY=$HOME/.conda/envs/liranenv/bin/python
[ -x "$PY" ] || PY=$HOME/.conda/envs/liranenv_gpu/bin/python
"$PY" plot_dcre_combined_MLvsARGv4.py --outdir dcre_csv --figdir . \
      --backgrounds org,mid --composite-backgrounds org,mid,mid2 \
      --composites Nd,LWP --stack-legend top \
      --stem dCRE_combined_4panel_MLvsARGv6 "$@"
