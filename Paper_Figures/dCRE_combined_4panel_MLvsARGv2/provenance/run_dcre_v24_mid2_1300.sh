#!/bin/bash
# ─────────────────────────────────────────────────────────────────────────────
# run_dcre_v24_mid2_1300.sh
#
# Regenerates the dCRE record with the ARG cases taken from WRF_dm_v2 instead
# of WRF_dm, window pinned to 2017-07-15 13:00, then redraws the combined
# 4-panel ML-vs-ARG figure from the CSVs it just wrote.
#
#   sbatch run_dcre_v24_mid2_1300.sh      # batch (skx)
#   ./run_dcre_v24_mid2_1300.sh           # nohup on the current node
#
# WHY WRF_dm_v2
#   WRF_dm's Registry omitted re_cloud/re_ice/re_snow from the morr_two_moment
#   package, so TAU_QC_TOT and RE_QC were identically zero and the ARG
#   Nc/LWP/re/cov apportionment had to fall back to the analytic 1/3:2/3
#   Twomey split.  WRF_dm_v2 populates both fields, so those sub-terms become
#   real measurements.  The plotting code detects this from R_tau rather than
#   assuming it, so pointing WRF_ARG_DIR back at WRF_dm restores the old
#   greyed-out treatment automatically.
#
# WINDOW
#   13:00 is inside every record: the eight ML/old-ARG runs reach 14:02+, and
#   the four WRF_dm_v2 runs reach 13:04-13:20.  The shared-time intersection is
#   therefore capped by WRF_TEND, not by a run that stopped short.
#
# Output directory: $WRF_OUTDIR (default ./dcre_figures_v24_mid2_1300).
# ─────────────────────────────────────────────────────────────────────────────
#SBATCH -J dcre_v24_mid2
#SBATCH -A TG-EES250031
#SBATCH -p skx
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 08:00:00
#SBATCH -o dcre_v24_mid2.%j.out
#SBATCH -e dcre_v24_mid2.%j.err

set -euo pipefail

export WORKDIR="/scratch/07088/tg863871/WRF_liran/JP"
export GEN="trace_east_boundary_subdomains_lwp_v24_mid2.py"
export PLOT2="plot_dcre_combined_6panel_mid2.py"
export PY="/scratch/07088/tg863871/conda_copy/envs/liranenv/bin/python"
export WRF_TEND="${WRF_TEND:-2017-07-15_13:00:00}"
export WRF_ARG_DIR="${WRF_ARG_DIR:-/scratch/07088/tg863871/WRF_dm_v2/test}"
export WRF_OUTDIR="${WRF_OUTDIR:-${WORKDIR}/dcre_figures_v24_mid2_1300}"

cd "$WORKDIR"
mkdir -p "$WRF_OUTDIR"

run_all() {
    set -euo pipefail
    cd "$WORKDIR"
    echo "================================================================"
    echo "  Generator  : ${GEN}"
    echo "  ARG tree   : ${WRF_ARG_DIR}"
    echo "  Window end : ${WRF_TEND}"
    echo "  Output dir : ${WRF_OUTDIR}"
    echo "  Started    : $(date)"
    echo "================================================================"
    "$PY" -u "$GEN"
    echo
    echo "===== replotting the combined 4-panel ML vs ARG-v2 figure ====="
    "$PY" -u "$PLOT2" --outdir "$WRF_OUTDIR"
    echo "Finished at $(date)."
}

STAMP=$(date +%Y%m%d_%H%M%S)

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    LOG="${WORKDIR}/dcre_v24_mid2_slurm_${SLURM_JOB_ID}.log"
    run_all 2>&1 | tee "$LOG"
else
    LOG="${WORKDIR}/dcre_v24_mid2_nohup_${STAMP}.log"
    PIDFILE="${WORKDIR}/dcre_v24_mid2_nohup.pid"
    nohup bash -c "$(declare -f run_all); run_all" > "$LOG" 2>&1 &
    echo "$!" > "$PIDFILE"
    echo "  Launched with nohup, PID $(cat "$PIDFILE")"
    echo "  Log     : ${LOG}"
    echo "  Results : ${WRF_OUTDIR}"
    echo "  Stop    : kill \$(cat ${PIDFILE})"
fi
