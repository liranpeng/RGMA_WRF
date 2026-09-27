#!/bin/bash
#SBATCH -J conseq_v6
#SBATCH -o conseq_v6.%j.out
#SBATCH -e conseq_v6.%j.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p skx-dev
#SBATCH -t 02:00:00
#SBATCH -A TG-EES250031

set -euo pipefail

cd "$SLURM_SUBMIT_DIR"

# conda is not initialised in batch shells here -- call the interpreter directly
PY=$HOME/.conda/envs/liranenv/bin/python

# v6 = v5 panels + the four new mid2 runs (ML 20170714_{1,3}aer_mid2,
# ARG 20170714_mid2_{1,3}aer_ARG).  The eight original cases seed from the
# v4 cache; only the mid2 files are read off Lustre.
$PY plot_arg_vs_ml_consequences_rwp_v6.py --workers 48 --font-scale 1.6
