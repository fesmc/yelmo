# Passive tracers

The passive-tracer subsystem (`calc_ytrc` in
[`src/yelmo_tracers.f90`](https://github.com/fesmc/yelmo/blob/main/src/yelmo_tracers.f90),
group `&ytrc`) traces the deposition time of the ice, from which the depths of
isochronal layers follow. It does not feed back on the flow. Three backends
can run in any combination:

| Backend | Switch | Method | Output |
|---|---|---|---|
| Eulerian | `use_euler`, `calc_age` | Advection of the deposition time on the model grid (`calc_tracer_3D`), with `uz_star` | `t_dep_euler` |
| Lagrangian particles | `use_tracer` | Particle model ([tracer](https://github.com/fesmc/tracer)), with `ux`, `uy`, `uz`; namelist `tracer_nml` | `t_dep_trc` |
| Lagrangian layers | `use_elsa` | Layer model ([elsa](https://github.com/fesmc/elsa)), driven by the ice thickness, velocity and surface and basal mass balance; namelist `elsa_nml`, group `elsa_group` | `t_dep_elsa` |

All backends are off by default. Each one gives a deposition-time field on
the Yelmo vertical grid. The backend selected by `t_dep_source` (`"euler"`,
`"trc"` or `"elsa"`) is copied to `t_dep`, from which the depths of the
isochrones at the times `time_iso` [ka] are diagnosed (`depth_iso`; with
`calc_age = False`, only the first time is used). The Eulerian solver is
explicit (`tracer_method = "expl"`) or implicit in the vertical with an
artificial diffusion `tracer_impl_kappa` (`"impl"`); the horizontal advection
is explicit upwind in both cases. The elsa layer stack is sized by the
simulation end time `time_end`. On a restart, the state of the Lagrangian
backends is read from separate files next to the restart file
(`<restart>_tracer.nc`, and `<restart>_elsa.nc` with `elsa_restart = True`).

The tracers are updated after the material properties and before the
thermodynamics in each [time step](timestepping.md). The variables are listed
in the [tracer variable table](../yelmo-variables-ytrc.md).
