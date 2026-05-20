# Basal friction optimization

A simple optimization program has developed that attempts
to optimize the basal friction field applied in Yelmo so
that the errors between simulated and observed ice thickness
are minimized.

Program: `tests/yelmo_opt.f90`
To compile: `make opt`
To run:

```bash
./runme -rs -e opt -o output/test -n par/yelmo_Antarctica_opt.nml
```

The program consists of the following steps:

## 1. Spin-up a steady-state ice sheet with constant forcing and fixed topography

For this step, the restart parameter should be set to `yelmo.restart='none'`, to
ensure that the spin-up is performed with the current parameters. Currently,
the program is hard-coded to spin-up the ice sheet for 20 kyr using SIA only,
followed by another 10 kyr using the solver of choice, as seen in the following lines of code:

```bash
call yelmo_update_equil_external(yelmo1,hyd1,cf_ref,time_init,time_tot=20e3,topo_fixed=.TRUE.,dt=5.0,ssa_vel_max=0.0)
call yelmo_update_equil_external(yelmo1,hyd1,cf_ref,time_init,time_tot=10e3, topo_fixed=.TRUE.,dt=1.0,ssa_vel_max=5000.0)
```

Note that this spin-up is obtained with a fixed topography set to the present-day
observed fields (`H_ice`,`z_bed`). After the spin-up finishes, a restart file is
written in the output directory with the name `yelmo_restart.nc`. The simulation
will terminate at this point.

## 2. Optimization

The restart file from Step 1 should be saved somewhere convenient for the model
(like in the `input` folder). Then the restart parameter should be set to that
location `yelmo.restart='PATH_TO_RESTART.nc'`. This will ensure that the spin-up
step is skipped, and instead the program will start directly with the optimization
iterations.

The optimization method follows Pollard and DeConto (2012), in that the basal
friction coefficient is scaled as a function of the error in elevation. Here we
do not modify `beta` directly, however, we assume that `beta = cf_ref * lambda_bed * N_eff * f(u)`.
`lambda_bed`, `N_eff` and `f(u)` are all controlled by parameter choices in the .nml file
like normal. Thus we are left with a unitless field `cf_ref`, which for any given
friction law varies within the range of about [0:1]. When `cf_ref=1.0`, sliding will
diminish to near zero, and `cf_ref~0.0` (near, but not zero) will give fast sliding.
This gives a convenient range for optimization.  

Parameters that control the total run time are hard coded:

- `qmax`: number of total iterations to run, where `qmax-1` is the number of optimization steps,
during which `cf_ref` is updated, and the last step is a steady-state run with `cf_ref` held constant.
- `time_iter`: time to run the model for each iteration before updating `cf_ref`.
- `time_steady`: Time to run the model to steady state with `cf_ref` held constant (last iteration step).

So, the program runs for, e.g., `time_iter=500` years with a given initial field of `cf_ref`
(with `C_bed` and `beta` updating every time step to follow changes in `u/v` and `N_eff`). At the end
of `time_iter`, the error in ice thickness is determined and used to update `cf_ref`
via the function `update_cf_ref_thickness_simple`. The model is again run for `time_iter` years
and the process is repeated.

Two important parameters control the optimization process: `tau` and `H_scale`.
The optimization works best when the ice shelves are relaxed to the reference
(observed) ice thickness in the beginning of the simulation, and then gradually
allowed to freely evolve. `tau` is the time scale of relaxation, which is applied
in Yelmo as `yelmo1%tpo%par%topo_rel_tau`. A lower value of `tau` means that the
ice shelves are more tightly held to the observed thickness. Likewise, `H_scale` controls the scaling
of the ice thickness error, which determines how to modify `cf_ref` at each iteration.
A higher value of `H_scale` means that changes to `cf_ref` will be applied more slowly.

These parameters are designed to change over time with the simulation. `tau` is
set to `rel_tau1` from the start of the simulation until `rel_time1`. Between `rel_time1`
and `rel_time2`, `tau` is linearly scaled from the value of `rel_tau1` to `rel_tau2`.
Or, if `rel_q > 1`, then the scaling is non-linear with an exponent of `rel_q` (this helps
maintain small values of `tau` longer which seems to help keep errors low).
Once `rel_time2` is reached, relaxation in the model is disabled, and the ice shelves
are allowed to freely evolve. Analogously, `H_scale` is modified the same way:
it is constant at the value of `scale_H1` until `scale_time1`,
linearly scaled between `scale_time1` and `scale_time2`,
and then constant thereafter at the value of `scale_H2`. Increasing the value of
`H_scale` over time helps to avoid oscillations in the optimization procedure as
`cf_ref` approaches the best fit.

Finally, after `qmax-1` iterations or `time=(qmax-1)*time_iter`, `cf_ref` is held
constant, and the simulation runs for `time_steady` years to equilibrate the model
with the current conditions. This step minimizes drift in the final result and
confirms that the optimized `cf_ref` field works well.
