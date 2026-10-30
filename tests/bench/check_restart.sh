#!/bin/bash
# check_restart.sh
#
# Restart round trip: ISLAND4 (32 km) runs 0-200 yr continuously, with a restart at 100 yr; a second run
# restarts from it and runs 100-200 yr. The restarts at 200 yr must agree bit for bit in every variable,
# including pc_tau_max and the window pc_taus it is taken from. Fails on any difference.
#
# Usage (from the yelmo root, after `make bench`; needs runme and h5diff):
#   tests/bench/check_restart.sh [outdir] [fixture]
# Defaults: outdir = output/check_restart, fixture = input/bench/island4-32km.nc (make_fixture.jl island4 dx_km=32)

set -u
OUT=$(mkdir -p "${1:-output/check_restart}" && cd "${1:-output/check_restart}" && pwd)
FIX=$(cd "$(dirname "${2:-input/bench/island4-32km.nc}")" && pwd)/$(basename "${2:-input/bench/island4-32km.nc}")
EXE=${EXE:-libyelmo/bin/yelmo_bench.x}
PAR=par/yelmo_bench_ISLAND4.nml
COMMON="bench.fixture=$FIX yelmo.grid_name=ISLAND4-32KM"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-2}

run() {  # run <dir> <params...>
    d=$1; shift
    rm -rf "$d"
    runme -r -e $EXE --omp "$OMP_NUM_THREADS" -n $PAR -o "$d" -p $COMMON "$@" > /dev/null 2>&1
    while pgrep -f "$(basename $EXE)" > /dev/null; do sleep 1; done
}

run "$OUT/cont"  ctrl.time_end=200 ctrl.dt_restart=100
run "$OUT/chunk" ctrl.time_init=100 ctrl.time_end=200 yelmo.restart_H_ice=True yelmo.restart_z_bed=True \
                 yelmo.restart="$OUT/cont/yelmo_restart_100.nc"

A=$OUT/cont/yelmo_restart.nc; B=$OUT/chunk/yelmo_restart.nc
[ -f "$A" ] && [ -f "$B" ] || { echo "FAIL: restart files missing in $OUT"; exit 1; }
ncdump -h "$A" | grep -q "pc_taus" || { echo "FAIL: pc_taus not in the restart file"; exit 1; }

diffs=$(h5diff -r "$A" "$B" 2>&1 | awk '/^dataset:/{d=$2} /^[0-9]+ differences found/{if ($1>0) print "  " d ": " $1 " differences"}')
if [ -n "$diffs" ]; then
    echo "FAIL: restart cont vs chunk differ:"; echo "$diffs"; exit 1
fi
echo "PASS: restart cont vs chunk identical (pc_tau_max, pc_taus and all other variables)"
