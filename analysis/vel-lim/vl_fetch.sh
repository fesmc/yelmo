#!/bin/bash
# Fetch vel-lim test runs from Levante (yelmo-vellim/output/vl/<run>) into ./data/<run>.
P=$(cd $(dirname $0) && pwd)/data; R=/work/ba1442/robinson/models/yelmo-vellim/output/vl
for r in "${@:-dev_clip br_clip br_drag br_drag6k br_drag8k br_drag10k}"; do for k in $r; do mkdir -p $P/$k
  for f in timesteps.nc yelmo_ts.nc yelmo.nc; do scp -q levante:$R/$k/$f $P/$k/ 2>&1 | grep -v module; done
  ssh -o BatchMode=yes levante "grep -E '^ssa: *[CX] ' $R/$k/out.out" > $P/$k/ssa_end.txt 2>/dev/null; done; done
