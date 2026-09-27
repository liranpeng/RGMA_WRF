#!/bin/bash
#SBATCH -J parcel_cmp
#SBATCH -A TG-EES250031
#SBATCH -p skx
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH -o %x_%j.out
#SBATCH -e %x_%j.err
# ---------------------------------------------------------------------------
# Stampede3 port of the NERSC parcel-vs-emulator comparison.
#   Run the OFFLINE Fortran parcel box model (libactivate.so) on the WRF states
#   this case actually visited, and compare against the ML fraction WRF applied
#   (FN11/FN21), producing parcel_vs_wrfout_<case>_scatter.png/eps + _hist.
#
# sbatch this FROM inside a case's online/ dir, e.g.:
#   cd .../test/20170714_1aer_org/online && sbatch submit_parcel_compare_stam3.slurm
# compare_parcel_vs_wrfout.py lives here; TEST_DIR = the parent case dir, so
# --all globs THIS case's wrfout_d01_*.
#
# Single node (NSHARDS=1). SKX = 48 physical cores / 96 hardware threads; the
# Pool uses all allocatable cores. Measured ~2 samp/s (48 workers) .. ~3.5 samp/s
# (96 HT) per node on these stiff near-saturation cells. DEADLINE solves a seeded
# RANDOM subset of the updraft pool for that many seconds and stops, so the shard
# always writes and merge can plot. 12000 s -> ~30k-42k solved cells (plenty for
# the hexbin scatter). Raise -N / DEADLINE for more points.
# ---------------------------------------------------------------------------
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"

module load gcc/13.2.0                # libgfortran.so.5 for libactivate.so
export OMP_NUM_THREADS=1              # 1 thread/worker; parallelism is across samples

PY=$HOME/.conda/envs/liranenv/bin/python
LIB=/scratch/07088/tg863871/Perlm_Backup/retrain_wrfdist_stam3/parcel_lib/libactivate.so
SCRIPT="$SLURM_SUBMIT_DIR/compare_parcel_vs_wrfout.py"
INPUTS="$SLURM_SUBMIT_DIR/parcel_inputs.npz"

CASE=$(basename "$(dirname "$SLURM_SUBMIT_DIR")")
OUT_PREFIX="$SLURM_SUBMIT_DIR/parcel_vs_wrfout_${CASE}"

export NSHARDS=1
export DEADLINE=12000                 # seconds of parcel solving (< walltime)
export MAXCELLS=1500000               # cap the reconstructed pool npz size

rm -f parcel_shard_*.npz

if [ -f "$INPUTS" ]; then
    echo "=== reuse existing $INPUTS (skip reconstruct)  $(date) ==="
else
    echo "=== phase reconstruct ($CASE)  $(date) ==="
    srun -N 1 -n 1 --cpu-bind=none "$PY" -u "$SCRIPT" --phase reconstruct \
        --all --max_cells "$MAXCELLS" --inputs_npz "$INPUTS"
fi

echo "=== phase parcel (1 node, deadline=${DEADLINE}s)  $(date) ==="
srun -N 1 -n 1 --cpu-bind=none "$PY" -u "$SCRIPT" --phase parcel \
    --inputs_npz "$INPUTS" --nshards "$NSHARDS" --shard 0 \
    --deadline "$DEADLINE" --shard_npz "parcel_shard_0.npz" --lib "$LIB"

echo "=== phase merge  $(date) ==="
srun -N 1 -n 1 --cpu-bind=none "$PY" -u "$SCRIPT" --phase merge \
    --inputs_npz "$INPUTS" --shard_glob "$SLURM_SUBMIT_DIR/parcel_shard_*.npz" \
    --out_prefix "$OUT_PREFIX"

echo "=== done  $(date) ==="
