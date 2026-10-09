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

## Masks

`&yelmo_masks` makes Yelmo's masks once at init from the regions, zones and
basins of [FesmData](https://github.com/fesmc/FesmData) (fesm-utils
`regions`), with selection expressions such as
`"region:Greenland & ~zone:open_ocean"`:

| Key | Meaning |
|---|---|
| `regions_group` | Group of the regions (`path_regions`, `path_basins`, `basin_sets`, `masks`, `mask_<name>`); `"None"` = no regions |
| `basins` | Basin ids of `bnd%basins`: `"<set>"` or `"<set>.group"` (e.g. `"Zwally2012.group"`); `"None"` = none, `"domain"` = one basin (1) |
| `mask_ice_dynamic` | Where ice is dynamic; elsewhere no ice (`bnd%mask_ice`) |
| `mask_ice_fixed` | Where ice thickness is prescribed (over `mask_ice_dynamic`) |
| `relax`, `relax_tau` | Where ice relaxes to `H_ice_ref` and the timescale [yr] (`bnd%tau_relax`, used with `ytopo.topo_rel = -1`) |
| `mask_rmse` | Where the error metrics are computed |

Without regions the expressions can only be `"all"` or `"none"`. A driver may
pass the regions (`yelmo_init(..., reg=)`, on the Yelmo grid) instead of
`regions_group`, and `mask_ice` instead of the mask_ice expressions. The named
masks of the regions (`masks`, `mask_<name>`) become Yelmo's regional output
domains. `bnd%regions` holds the region codes of the deepest level.

::: {.callout-tip}
To browse, build, compare and validate parameter files from the command line, see
[`yelmo-config`](yelmo-config.md).
:::
