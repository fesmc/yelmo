#!/bin/bash
# Benchmark and test suite (Robinson et al., GMD 2020, plus later benchmarks).
# Run from the yelmo root after `configme` and `runme config init`.
# Jobs are submitted with runme (-s) to the queue alias `q` of the cluster set
# in .runme/config.toml (see `runme queues`); use `-r` alone to run locally.

fldr='output/yelmo-bench'

q='short'                       # queue alias of your cluster, e.g. short (pik_hpc2024), 12h (awi_albedo), shared (dkrz_levante)
runopt="-rs -q ${q}"

### BENCHMARK TESTS ###

make benchmarks

# EISMINT1 moving margin, EISMINT2 EXPA and EXPF
runme ${runopt} -w 01:00:00 -e benchmarks -o ${fldr}/moving -n par/yelmo_EISMINT_moving.nml
runme ${runopt} -w 05:00:00 -e benchmarks -o ${fldr}/expa   -n par/yelmo_EISMINT_expa.nml
runme ${runopt} -w 05:00:00 -e benchmarks -o ${fldr}/expf   -n par/yelmo_EISMINT_expf.nml -p ctrl.dt2D_out=1000

# EISMINT1-moving margin with DIVA solvers
# Note ydyn.solver="diva-noslip" is broken, as well as too high friction values, eg ydyn.beta_const=1e6
runme ${runopt} -w 01:00:00 -e benchmarks -o ${fldr}/moving-diva-noslip -n par/yelmo_EISMINT_moving.nml -p ydyn.solver="diva-noslip" ctrl.time_end=30e3 ctrl.dt2D_out=200
runme ${runopt} -w 01:00:00 -e benchmarks -o ${fldr}/moving-diva -n par/yelmo_EISMINT_moving.nml -p ydyn.solver="diva" ydyn.beta_method=0 ydyn.beta_const=1e4 ctrl.time_end=30e3 ctrl.dt2D_out=200

# EISMINT1 EXPA-like domain with sloping bed and shelves possible for testing symmetry (not part of GMD suite of tests)
# In tests/yelmo_benchmarks.f90, set the hard-coded parameter `sym_test_with_shelves=.TRUE.` and recompile.
# Fixed topo and fixed beta; fixed topo; dynamic simulation:
runme -r -e benchmarks -o ${fldr}/moving-float-diva-fixed-beta -n par/yelmo_EISMINT_moving.nml -p ctrl.time_end=1e3 ctrl.dt2D_out=200 ydyn.solver=diva ytopo.topo_fixed=True ydyn.beta_method=0 ydyn.beta_const=1e4
runme -r -e benchmarks -o ${fldr}/moving-float-diva-fixed-topo -n par/yelmo_EISMINT_moving.nml -p ctrl.time_end=1e3 ctrl.dt2D_out=200 ydyn.solver=diva ytopo.topo_fixed=True
runme -r -e benchmarks -o ${fldr}/moving-float-diva -n par/yelmo_EISMINT_moving.nml -p ctrl.time_end=1e3 ctrl.dt2D_out=200 ydyn.solver=diva

# Ensemble of HALFAR simulations with various values of
# dx to test numerical convergence with analytical solution
runme ${runopt} -w 05:00:00 -e benchmarks -n par/yelmo_HALFAR.nml -o ${fldr}/halfar -p ctrl.dx=0.5,1.0,2.0,3.0,4.0,5.0,8.0

# Ensemble of EISMINT1-moving simulations with various values of
# dx and pc_eps to test adaptive timestepping
runme ${runopt} -w 05:00:00 -e benchmarks -n par/yelmo_EISMINT_moving.nml -o ${fldr}/moving_dts \
    -p ctrl.time_end=25e3 yelmo.log_timestep=True ytherm.method="fixed" ctrl.dx=5.0,10.0,25.0,50.0,60.0 yelmo.pc_eps=1e-2,1e-1,1e0


### INITMIP TESTS ###

make initmip

# Antarctica present-day initialization (single run, 32 km).
# For LGM forcing, set ctrl.set_nm="set_ant_lgm".
# For higher resolutions, change yelmo.grid_name to "ANT-16KM" or "ANT-8KM" (4 km not yet supported).
runme ${runopt} -w 05:00:00 -e initmip -n par/yelmo_initmip.nml -o ${fldr}/initmip-ant-32km \
    -p ctrl.dtt=5 ctrl.time_end=1e3 \
       ctrl.set_nm="set_ant_pd" yelmo.log_timestep=True \
       ydyn.solver="diva" yelmo.domain="Antarctica" yelmo.grid_name="ANT-32KM"

# Antarctica resolution ensemble (32 / 16 / 8 km).
runme ${runopt} -w 05:00:00 -e initmip -n par/yelmo_initmip.nml -o ${fldr}/initmip-ant-ens \
    -p ctrl.dtt=5 ctrl.time_end=1e3 \
       ctrl.set_nm="set_ant_pd" yelmo.log_timestep=True \
       ydyn.solver="diva" yelmo.domain="Antarctica" \
       yelmo.grid_name="ANT-32KM","ANT-16KM","ANT-8KM"

# Greenland present-day initialization (single run, 32 km).
# For higher resolutions, change yelmo.grid_name to "GRL-16KM", "GRL-8KM", or "GRL-4KM".
runme ${runopt} -w 05:00:00 -e initmip -n par/yelmo_initmip.nml -o ${fldr}/initmip-grl-32km \
    -p ctrl.dtt=5 ctrl.time_end=1e3 \
       ctrl.set_nm="set_grl_pd" yelmo.log_timestep=True \
       ydyn.solver="diva" yelmo.domain="Greenland" yelmo.grid_name="GRL-32KM"

# Greenland resolution ensemble (32 / 16 / 8 / 4 km).
runme ${runopt} -w 05:00:00 -e initmip -n par/yelmo_initmip.nml -o ${fldr}/initmip-grl-ens \
    -p ctrl.dtt=5 ctrl.time_end=1e3 \
       ctrl.set_nm="set_grl_pd" yelmo.log_timestep=True \
       ydyn.solver="diva" yelmo.domain="Greenland" \
       yelmo.grid_name="GRL-32KM","GRL-16KM","GRL-8KM","GRL-4KM"

# Solver stability (Robinson et al., 2022): Greenland with DIVA, one run and all resolutions
runme ${runopt} -w 05:00:00 -e initmip -n par/yelmo_initmip.nml -o ${fldr}/grl-diva-test \
    -p ctrl.dtt=5 ctrl.time_end=1e3 ctrl.set_nm="set_grl_pd" yelmo.domain="Greenland" \
       yelmo.log_timestep=True ydyn.solver="diva" yelmo.grid_name="GRL-16KM"
runme ${runopt} -w 05:00:00 -e initmip -n par/yelmo_initmip.nml -o ${fldr}/grl-diva \
    -p ctrl.dtt=5 ctrl.time_end=1e3 ctrl.set_nm="set_grl_pd" yelmo.domain="Greenland" \
       yelmo.log_timestep=True ydyn.solver="diva" yelmo.grid_name="GRL-32KM","GRL-16KM","GRL-8KM","GRL-4KM"

# OpenMP scaling (GRL-8KM, 1-32 threads)
make clean
make initmip openmp=1
gridname='GRL-8KM'
for nt in 1 2 4 8 16 32; do
    runme ${runopt} -w 01:00:00 --omp ${nt} -e initmip -n par/yelmo_initmip.nml -o ${fldr}/openmp/${gridname}-omp$(printf "%02d" ${nt}) \
        -p yelmo.grid_name=${gridname} ctrl.dtt=5 ctrl.time_end=1e3 ctrl.set_nm="set_grl_pd" \
           yelmo.domain="Greenland" yelmo.log_timestep=True
done
make clean

### MISMIP TESTS ###

make mismip

# MISMIP3D rate factor (RF) experiment, three grounding-line treatments at four resolutions:
runme ${runopt} -w 24:00:00 -e mismip -n par/yelmo_MISMIP3D.nml -o ${fldr}/mismip/default -p ctrl.experiment="RF" ydyn.beta_gl_scale=0 ydyn.beta_gl_stag=0 ctrl.dx=2.5,5.0,10.0,20.0
runme ${runopt} -w 24:00:00 -e mismip -n par/yelmo_MISMIP3D.nml -o ${fldr}/mismip/subgrid -p ctrl.experiment="RF" ydyn.beta_gl_scale=0 ydyn.beta_gl_stag=3 ctrl.dx=2.5,5.0,10.0,20.0
runme ${runopt} -w 24:00:00 -e mismip -n par/yelmo_MISMIP3D.nml -o ${fldr}/mismip/scaling -p ctrl.experiment="RF" ydyn.beta_gl_scale=2 ydyn.beta_gl_stag=3 ctrl.dx=2.5,5.0,10.0,20.0

### TROUGH, MISMIP+ and SLAB-S06 (yelmo_trough.x) ###

make trough

# MISMIP+
runme ${runopt} -w 05:00:00 -e trough -n par/yelmo_MISMIP+.nml -o ${fldr}/mismip+

# MISMIP+ ensemble (hybrid, diva):
runme ${runopt} -w 05:00:00 -e trough -n par/yelmo_MISMIP+.nml -o ${fldr}/mismip+-solver -p ydyn.solver="hybrid","diva"

# TROUGH-F17 (Feldmann and Levermann, 2017)
runme ${runopt} -w 05:00:00 -e trough -n par/yelmo_TROUGH-F17.nml -o ${fldr}/trough

# ssa solver; resolutions 1 and 2 km; cf_ref ensemble with beta_u0=100:
runme ${runopt} -w 05:00:00 -e trough -n par/yelmo_TROUGH-F17.nml -o ${fldr}/trough-ssa -p ydyn.solver="ssa"
runme ${runopt} -w 48:00:00 -e trough -n par/yelmo_TROUGH-F17.nml -o ${fldr}/trough-dx1 -p ctrl.dx=1.0
runme ${runopt} -w 24:00:00 -e trough -n par/yelmo_TROUGH-F17.nml -o ${fldr}/trough-dx2 -p ctrl.dx=2.0
runme ${runopt} -w 05:00:00 -e trough -n par/yelmo_TROUGH-F17.nml -o ${fldr}/trough-u0.100-cf -p ydyn.beta_u0=100 ytill.cf_ref=5.0,10.0,20.0

# SLAB-S06 (Schoof, 2006): one test run, then the dx ensemble
runme -r -e trough -n par/yelmo_SLAB-S06.nml -o ${fldr}/slab06-test -p ctrl.dx=4 ydyn.ssa_iter_max=10
runme ${runopt} -w 05:00:00 -e trough -n par/yelmo_SLAB-S06.nml -o ${fldr}/slab06 -p ctrl.dx=0.5,1,2,4,8

# Constant viscosity
runme ${runopt} -w 05:00:00 -e trough -n par/yelmo_SLAB-S06.nml -o ${fldr}/slab06-visc -p ctrl.dx=1,2,4,8 ydyn.visc_method=0 ydyn.visc_const=2e6

# Different friction-law exponent
runme ${runopt} -w 05:00:00 -e trough -n par/yelmo_SLAB-S06.nml -o ${fldr}/slab06-q0.5 -p ctrl.dx=0.5,1,2,4,8 ydyn.beta_q=0.5

### ISMIP-HOM ###

make ismiphom

runme ${runopt} -w 01:00:00 -e ismiphom -n par/yelmo_ISMIPHOM.nml -o ${fldr}/ismiphom/expa/diva   -p ctrl.experiment="EXPA" ctrl.L=5,10,20,40,80,160 ydyn.solver="diva"
runme ${runopt} -w 01:00:00 -e ismiphom -n par/yelmo_ISMIPHOM.nml -o ${fldr}/ismiphom/expa/hybrid -p ctrl.experiment="EXPA" ctrl.L=5,10,20,40,80,160 ydyn.solver="hybrid"
runme ${runopt} -w 01:00:00 -e ismiphom -n par/yelmo_ISMIPHOM.nml -o ${fldr}/ismiphom/expc/diva   -p ctrl.experiment="EXPC" ctrl.L=5,10,20,40,80,160 ydyn.solver="diva"
runme ${runopt} -w 01:00:00 -e ismiphom -n par/yelmo_ISMIPHOM.nml -o ${fldr}/ismiphom/expc/hybrid -p ctrl.experiment="EXPC" ctrl.L=5,10,20,40,80,160 ydyn.solver="hybrid"

# Experiment F (EXPF1 no slip, EXPF2 slip ratio 1; L is fixed to 100 km by the driver)
runme ${runopt} -w 05:00:00 -e ismiphom -n par/yelmo_ISMIPHOM.nml -o ${fldr}/ismiphom/expf -p ctrl.experiment="EXPF1","EXPF2" ctrl.time_end=2000 ctrl.dtt=10 ctrl.dt2D_out=100

### CalvingMIP ###

make calving

runme ${runopt} -w 05:00:00 -e calving -n par/yelmo_calvingmip.nml -o ${fldr}/calvmip-exp01 -p ctl.exp="exp1" ctl.dt2D_out=200 ctl.time_end=10e3

### AGE TESTS ###

# tests/test_icetemp.f90 (not run through runme). Arguments: nz enth_cr solver experiment.
# experiment="eismint" uses t_end = 300 kyr, dt = 0.5 yr, dt_out = 1 kyr,
# T_pmp_beta = 9.7e-8 K Pa^-1 and init_eismint_summit(smb=0.5).
# Precision: set `wp = dp` in src/yelmo_defs.f90 (default `wp = sp`) and recompile.
# Run 3 experiments (nz=12,32,52), corresponding to 10, 30 and 50 internal ice points.
# All output goes to nc files in the folder "output/".
make clean
make icetemp
./libyelmo/bin/test_icetemp.x 12 1e-3 temp eismint
./libyelmo/bin/test_icetemp.x 32 1e-3 temp eismint
./libyelmo/bin/test_icetemp.x 52 1e-3 temp eismint
# Then set `wp = sp` in src/yelmo_defs.f90 and run the nz=32 case again.
make clean
make icetemp
./libyelmo/bin/test_icetemp.x 32 1e-3 temp eismint
