#!/bin/bash
#
# run_trc_eismint.sh — EISMINT-moving with the Lagrangian particle tracer (%trc)
# and its gridded stats, activated entirely from the command line via `runme -p`.
#
# No dedicated namelist file is needed: the tracer's full &trc schema lives in
# input/tracer_defaults.nml (a copy synced from tracer/input/), so only the
# handful of switches below are overridden. The base par carries a dormant &trc
# group (just `stats`) because `runme -p` overrides existing keys but does not
# create them.
#
# Usage:  scripts/run_trc_eismint.sh [OUTDIR]
#
# The tracer's gridded stats (trc_count on depth_norm, trc_depth_iso on time_iso)
# are written into <OUTDIR>/yelmo_restart.nc via the restart path. The benchmark
# driver's 2D writer emits a fixed variable list, so they do not appear in
# yelmo2D.nc unless a driver requests them by name.
#
# Note on time_iso: the isochrone targets come from ytrc%time_iso (default in
# input/yelmo_defaults.nml). Arrays cannot be set with `runme -p` (a comma list
# is read as an ensemble), so adjust them in the par/defaults if needed. The
# defaults are paleo deposition times; on a short forward EISMINT run they match
# no particles, so trc_depth_iso is mostly missing while trc_count is populated.

set -e

OUTDIR=${1:-output/trc-eismint}

make benchmarks

runme -r -e benchmarks -o "${OUTDIR}" -n par/yelmo_EISMINT_moving.nml \
    -p ctrl.time_end=2000 ctrl.dt2D_out=1000 \
       ytrc.use_euler=True ytrc.use_tracer=True ytrc.calc_age=True \
       ytrc.tracer_nml=yelmo_EISMINT_moving.nml \
       trc.stats=True

echo "Done. Tracer gridded stats are in ${OUTDIR}/yelmo_restart.nc"
echo "  ncdump -h ${OUTDIR}/yelmo_restart.nc | grep -E 'depth_norm|trc_count|trc_depth_iso'"
