#!/bin/bash
# Run a CalvingMIP steady-state experiment (exp1 or exp3) followed by the
# perturbation experiment (exp2 or exp4) that restarts from it.
#
# Usage:
#     scripts/run_calvingmip_chain.sh <output-folder> <ss-exp> <pert-exp> [<dx-km>]
#
# Examples:
#     scripts/run_calvingmip_chain.sh tmp/yelmo-calvingmip-2026-05-18 exp1 exp2
#     scripts/run_calvingmip_chain.sh tmp/yelmo-calvingmip-2026-05-18 exp3 exp4 10
#
# Notes:
#   - Run from the yelmo repository root (where ./runme lives).
#   - Requires libyelmo/bin/yelmo_calving.x (`make calving`).
#   - Each experiment's stdout is written to <output-folder>/<exp>/run.log.
#   - <output-folder> may be relative; an absolute form is used for
#     yelmo.restart so the perturbation run resolves it correctly.
#   - <dx-km>, if given, overrides ctl.dx for both runs.

set -euo pipefail

if [[ $# -lt 3 || $# -gt 4 ]]; then
    echo "Usage: $0 <output-folder> <ss-exp> <pert-exp> [<dx-km>]" >&2
    exit 2
fi

OUTFLDR="$1"
SS_EXP="$2"
PERT_EXP="$3"
DX="${4:-}"

DX_ARG=()
if [[ -n "$DX" ]]; then
    DX_ARG=("ctl.dx=$DX")
fi

case "$OUTFLDR" in
    /*) ABS_OUTFLDR="$OUTFLDR" ;;
    *)  ABS_OUTFLDR="$(pwd)/$OUTFLDR" ;;
esac

# --- Steady-state experiment (exp1 or exp3) ---------------------------------
./runme -e calving -o "$OUTFLDR/$SS_EXP" -n par/yelmo_calvingmip.nml \
        -p ctl.exp="$SS_EXP" ycalv.calv_flt_method="exp1" "${DX_ARG[@]}"

(
    cd "$ABS_OUTFLDR/$SS_EXP"
    ./yelmo_calving.x yelmo_calvingmip.nml > run.log 2>&1
)

# --- Perturbation experiment (exp2 or exp4) restarting from steady state ---
RESTART="$ABS_OUTFLDR/$SS_EXP/yelmo_restart.nc"

./runme -e calving -o "$OUTFLDR/$PERT_EXP" -n par/yelmo_calvingmip.nml \
        -p ctl.exp="$PERT_EXP" ctl.time_end=1000 ctl.dtt=1 ctl.dt2D_out=100 \
           yelmo.restart="$RESTART" yelmo.restart_z_bed=True yelmo.restart_H_ice=True \
           ycalv.calv_flt_method="exp2" "${DX_ARG[@]}"

(
    cd "$ABS_OUTFLDR/$PERT_EXP"
    ./yelmo_calving.x yelmo_calvingmip.nml > run.log 2>&1
)
