# Parameters

Yelmo reads its parameters from a namelist file, given to `yelmo_init` as
`filename`. The defaults are in `input/yelmo_defaults.nml`, which declares every
parameter together with its default value and a short description. A user
parameter file only needs to list the parameters it changes. Unknown parameters
stop the model at load time, as do values outside their allowed set or range.

The parameters are organised in groups, one per model component:

| Group | Component | Documentation |
|---|---|---|
| `&yelmo` | Domain, grid, restart, vertical grid, time stepping | [Time stepping](physics/timestepping.md) |
| `&ytopo` | Topography, mass conservation, front | [Mass conservation](physics/mass_conservation/index.md) |
| `&ycalv` | Calving | [Calving](physics/calving.md) |
| `&ydyn` | Momentum balance, basal friction | [Momentum balance](physics/momentum/diva.md), [Basal friction](physics/basal-friction.md) |
| `&ytill` | Till properties (`cb_ref`) | [Basal friction](physics/basal-friction.md#bed-coefficient) |
| `&ymat` | Rheology, enhancement | [Thermodynamics](physics/thermodynamics.md#coupling-to-the-flow) |
| `&ytrc` | Passive tracers (deposition time) | |
| `&ytherm` | Thermodynamics, bedrock | [Thermodynamics](physics/thermodynamics.md) |
| `&yhyd` | Basal hydrology (FastHydrology) | [Basal friction](physics/basal-friction.md#effective-pressure) |
| `&yelmo_masks`, `&yelmo_init_topo`, `&yelmo_data` | Input files for masks, initial topography and comparison data | |

Each group can be renamed in `&yelmo` (e.g. `nml_ydyn = "ydyn_north"` reads the
dynamics parameters from `&ydyn_north`), so that several Yelmo domains can share
one parameter file. Physical constants are read from
`input/yelmo_phys_const.nml`, in the group given by `yelmo.phys_const`
(default `"Earth"`; groups for the benchmark experiments, e.g. `"EISMINT"`, `"MISMIP3D"`, are also available).

Groups used only by a driver program, such as `&ctrl`, are not read by Yelmo.

::: {.callout-tip}
To browse, build, compare and validate parameter files from the command line, see
[`yelmo-config`](yelmo-config.md).
:::
