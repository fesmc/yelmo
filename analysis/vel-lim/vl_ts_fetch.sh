#!/bin/bash
# Extract the centreline (yc=0) of the 10-yr 2D output of the ts_* reruns on Levante and fetch it.
P=$(cd $(dirname $0) && pwd)/data; R=/work/ba1442/robinson/models/yelmo-vellim/output/vl
ssh -o BatchMode=yes levante "bash -lc 'module load nco; cd $R; for r in ts_*; do ncks -O -v uxy_bar,ux_bar,uxy_b,H_ice,f_grnd,z_srf -d yc,0.0 \$r/yelmo.nc \$r/centreline.nc; done'" 2>&1 | grep -v "module\|X11\|bashrc"
for r in ts_clip ts_drag5k ts_drag6k ts_drag8k ts_drag10k ts_drag15k ts_drag20k ts_drag30k ts_drag50k; do mkdir -p $P/$r; scp -q levante:$R/$r/centreline.nc $P/$r/; done
