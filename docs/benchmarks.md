# Benchmarks

Yelmo ships with a suite of standard ice-sheet benchmark experiments used for
verification and regression testing. They span analytic and idealized cases
(EISMINT, slab, trough) as well as community intercomparison protocols
(MISMIP3D, initMIP). Each benchmark has a dedicated test program under
[`tests/`](https://github.com/fesmc/yelmo/tree/main/tests) and a parameter
file under `par-gmd/` or `par/`.

## General workflow

1. **Build** the test program: `make <target>`.
2. **Run** it via `runme`, choosing the executable alias, namelist, and output
   directory:
   ```bash
   ./runme -r -e <alias> -n <namelist> -o <output-dir>
   ```
3. **Override** namelist values inline as needed with `-p group.key=value`.

For ensemble runs (one execution per parameter combination), use
[`jobrun`](https://github.com/alex-robinson/runner) from the `runner` Python
module — see `README_BATCH.md` for installation. Comma-separated values in
`-p` are interpreted by `jobrun` as ensemble dimensions:

```bash
jobrun ./runme -r -e <alias> -n <namelist> -o <ens-dir> -p key=val1,val2,val3
```

The sections below cover each benchmark in the current Yelmo test suite.
For additional variants (alternative solvers, basal-friction sweeps, OpenMP
scaling tests, and benchmarks not yet documented here such as HALFAR,
ISMIPHOM, MISMIP+, and CalvingMIP), see
[`run_batch_gmd.sh`](https://github.com/fesmc/yelmo/blob/main/run_batch_gmd.sh)
in the repository root.

## EISMINT1-moving

EISMINT Phase 1 moving-margin experiment ([Huybrechts and Payne, 1996](https://doi.org/10.3189/S0260305500013197)).
Idealized circular dome on a flat bed with a prescribed surface mass-balance
field that drives transient margin advance and retreat toward a quasi
steady state. Exercises the SIA flow law, mass conservation, and basic
margin handling.

```bash
make clean
make benchmarks
./runme -r -e benchmarks -n par-gmd/yelmo_EISMINT_moving.nml -o output/benchmarks/eismint-moving
```

## EISMINT2-expa

EISMINT Phase 2 Experiment A ([Payne et al., 2000](https://doi.org/10.3189/172756500781832891)).
Steady-state thermomechanical experiment with a fixed margin on a flat bed,
designed to test the coupling between ice dynamics and internal thermodynamics
— in particular the formation of cold-ice cores and basal-temperate regions.

```bash
make clean
make benchmarks
./runme -r -e benchmarks -n par-gmd/yelmo_EISMINT_expa.nml -o output/benchmarks/eismint-expa
```

## MISMIP3D

Marine Ice Sheet Model Intercomparison Project, 3D version
([Pattyn et al., 2013](https://doi.org/10.3189/2013JoG12J129)). Tests
grounding-line dynamics on a marine bed with retrograde slope, using a SSA
or DIVA momentum balance. Two experiments are defined in the namelist via
the `ctrl.experiment` field: `Stnd` (standard, advance to steady state) and
`RF` (reverse forcing, to test reversibility).

Standard experiment:

```bash
make clean
make mismip
./runme -r -e mismip -n par-gmd/yelmo_MISMIP3D.nml -o output/benchmarks/mismip3d-stnd -p ctrl.experiment=Stnd
```

Reverse-forcing experiment:

```bash
./runme -r -e mismip -n par-gmd/yelmo_MISMIP3D.nml -o output/benchmarks/mismip3d-rf -p ctrl.experiment=RF
```

### Resolution ensembles

To probe grid resolution and the treatment of the grounding line, run
ensembles over `ctrl.dx` with different basal-stress scaling /
staggering options at the grounding line:

```bash
# Default grounding-line treatment
jobrun ./runme -r -e mismip -n par-gmd/yelmo_MISMIP3D.nml -o output/benchmarks/mismip3d-default \
    -p ydyn.beta_gl_scale=0 ydyn.beta_gl_stag=0 ctrl.dx=2.5,5.0,10.0,20.0

# Subgrid grounding-line interpolation (beta_gl_stag=3)
jobrun ./runme -r -e mismip -n par-gmd/yelmo_MISMIP3D.nml -o output/benchmarks/mismip3d-subgrid \
    -p ydyn.beta_gl_scale=0 ydyn.beta_gl_stag=3 ctrl.dx=2.5,5.0,10.0,20.0

# Subgrid + basal-stress scaling (beta_gl_scale=2)
jobrun ./runme -r -e mismip -n par-gmd/yelmo_MISMIP3D.nml -o output/benchmarks/mismip3d-scaling \
    -p ydyn.beta_gl_scale=2 ydyn.beta_gl_stag=3 ctrl.dx=2.5,5.0,10.0,20.0
```

## slab

1D ice-slab benchmark on a tilted bed ([Schoof, 2006](https://doi.org/10.1017/S0022112006002977)),
used to test SSA / DIVA solver behaviour, basal-friction parameterizations,
and numerical convergence on a configuration with an analytic reference. The
domain is dispatched inside the `yelmo_trough` program based on the
`ctrl.domain` field of the namelist (`SLAB-S06`), so the `trough` build
target and executable alias are used here.

```bash
make clean
make trough
./runme -r -e trough -n par/yelmo_SLAB-S06.nml -o output/benchmarks/slab
```

Note: the `make slab` target builds a separate slab program
(`yelmo_slab.x`, used in Robinson, 2022) and is **not** the SLAB-S06
benchmark.

## trough-f17

Idealized trough geometry ([Feldmann and Levermann, 2017](https://doi.org/10.5194/tc-11-1745-2017))
designed to test buttressing and grounding-line behaviour in a confined,
narrowing channel. Used to probe stability and lateral-stress
representation in marine outlet glaciers.

```bash
make clean
make trough
./runme -r -e trough -n par/yelmo_TROUGH-F17.nml -o output/benchmarks/trough-f17
```

### Variants and ensembles

Alternative solver, resolution, and basal-friction settings used to
probe sensitivity:

```bash
# SSA solver instead of the default
./runme -r -e trough -n par/yelmo_TROUGH-F17.nml -o output/benchmarks/trough-f17-ssa -p ydyn.solver=ssa

# Higher resolution (1 km and 2 km)
./runme -r -e trough -n par/yelmo_TROUGH-F17.nml -o output/benchmarks/trough-f17-dx1 -p ctrl.dx=1.0
./runme -r -e trough -n par/yelmo_TROUGH-F17.nml -o output/benchmarks/trough-f17-dx2 -p ctrl.dx=2.0

# cf_ref ensemble (3 values) with adjusted beta_u0
jobrun ./runme -r -e trough -n par/yelmo_TROUGH-F17.nml -o output/benchmarks/trough-f17-cf \
    -p ydyn.beta_u0=100 ytill.cf_ref=5.0,10.0,20.0
```

## initmip-grl

Greenland initialization benchmark following the
initMIP-Greenland protocol ([Goelzer et al., 2018](https://doi.org/10.5194/tc-12-1433-2018)).
Spins the ice sheet up toward a present-day steady state using the
optimization-based initialization scheme (`equil_method = "opt"`) with
present-day boundary conditions. Default settings: 5-yr outer timestep,
1000-yr simulation, steady-state forcing.

Greenland is supported at 32, 16, 8, and 4 km resolution. Select the grid
by uncommenting the desired line:

```bash
grid=GRL-32KM
#grid=GRL-16KM
#grid=GRL-8KM
#grid=GRL-4KM

make clean
make initmip
./runme -r -e initmip -n par/yelmo_initmip.nml -o output/initmip-grl-$grid \
    -p ctrl.dtt=5 ctrl.time_end=1e3 ctrl.time_equil=100 \
       ctrl.set_nm=set_grl_pd yelmo.log_timestep=True \
       ydyn.solver=diva yelmo.domain=Greenland yelmo.grid_name=$grid
```

To run all four resolutions as an ensemble:

```bash
jobrun ./runme -r -e initmip -n par/yelmo_initmip.nml -o output/initmip-grl-ens \
    -p ctrl.dtt=5 ctrl.time_end=1e3 ctrl.time_equil=100 \
       ctrl.set_nm=set_grl_pd yelmo.log_timestep=True \
       ydyn.solver=diva yelmo.domain=Greenland \
       yelmo.grid_name=GRL-32KM,GRL-16KM,GRL-8KM,GRL-4KM
```

## initmip-ant

Antarctica initialization benchmark following the
initMIP-Antarctica protocol ([Seroussi et al., 2019](https://doi.org/10.5194/tc-13-1441-2019)).
Same optimization-based spin-up as initmip-grl, but on the Antarctic domain
with `set_nm = "set_ant_pd"`.

Antarctica is currently supported at 32, 16, and 8 km resolution. The 4 km
configuration is not yet available.

```bash
grid=ANT-32KM
#grid=ANT-16KM
#grid=ANT-8KM

make clean
make initmip
./runme -r -e initmip -n par/yelmo_initmip.nml -o output/initmip-ant-$grid \
    -p ctrl.dtt=5 ctrl.time_end=1e3 ctrl.time_equil=100 \
       ctrl.set_nm=set_ant_pd yelmo.log_timestep=True \
       ydyn.solver=diva yelmo.domain=Antarctica yelmo.grid_name=$grid
```

To run all three resolutions as an ensemble:

```bash
jobrun ./runme -r -e initmip -n par/yelmo_initmip.nml -o output/initmip-ant-ens \
    -p ctrl.dtt=5 ctrl.time_end=1e3 ctrl.time_equil=100 \
       ctrl.set_nm=set_ant_pd yelmo.log_timestep=True \
       ydyn.solver=diva yelmo.domain=Antarctica \
       yelmo.grid_name=ANT-32KM,ANT-16KM,ANT-8KM
```
