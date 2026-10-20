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
| `regions_group` | Group of the regions (`path_regions`, `path_basins`, `basin_sets`, `layers`, `masks`, ...; see fesm-utils `regions`); `"None"` = no layers |
| `basins` | Basin ids of `bnd%basins`: a layer of the regions, e.g. `"<set>"` or `"<set>.group"` (`"Zwally2012.group"`); `"None"` = none, `"domain"` = one basin (1) |
| `mask_ice_dynamic` | Where ice is dynamic; elsewhere no ice (`bnd%mask_ice`) |
| `mask_ice_fixed` | Where ice thickness is prescribed (over `mask_ice_dynamic`) |
| `relax`, `relax_tau` | Where ice relaxes to `H_ice_ref` and the timescale [yr] (`bnd%tau_relax`, used with `ytopo.topo_rel = -1`) |
| `border_width` | `N M`: a strip at the edge of the region where ice is allowed (`mask_ice_dynamic`, `mask_ice_fixed`): its outer `N` cells fixed (`H_ice = H_ice_ref`), the next `M` cells inward relaxing to `H_ice_ref` with `relax_tau` (needs `ytopo.topo_rel = -1`); `0 0` = none (below) |
| `mask_rmse` | Where the error metrics are computed |
| `dynamic_connected` | Keep only the dynamic cells connected to the ice sheet (the largest 4-connected region of present-day ice, else of the initial ice) through dynamic or fixed cells; the others get no ice (e.g. offshore islands) |

The strip of `border_width` is applied after the expressions, which stay for
hand-made zones. The distance of a cell to the edge of the region is counted in
cells in the Chebyshev (chessboard) metric, diagonal neighbours included, so the
strip has a constant width around corners and is rectangular on a rectangular
region; beyond the domain border is outside the region (except in periodic
directions). With `mask_ice_dynamic = "all"` the strip lies along the domain
border, e.g. for a domain whose border cuts through ice (GreenlandPaleo), with
`H_ice_ref` set in time by the driver. `yelmo.mask_border` is applied after the
masks group, on the outermost row only: `mask_border = "fixed"` is the special
case `border_width = 1 0` with `mask_ice_dynamic = "all"`; with a strip on the
domain border use `"fixed"` or `"dynamic"` (`"none"` would empty its outer row).
`"auto"` gives `"fixed"` for a domain it does not know.

With `experiment = "infinite"` and `mask_border = "dynamic"` the border cells
take the thickness of their inner neighbours. Ice flows out through the border
but not in: the border cells receive the outflow of their neighbours, and the
reset to their neighbours' thickness is booked in `mb_border` (the mass that
leaves the domain). The same holds for the non-periodic borders of
`"periodic-x"` and `"periodic-y"`. The border cells keep their mass balance
(smb, bmb, calving, ...), which the reset undoes; the totals of these terms in
`yelmo_ts.nc` leave the border cells out, and `mb_border_tot` is the total of
everything booked in them (their rates and the reset): the mass that leaves the
domain through its border, so that dV/dt = the terms + `mb_border_tot`.

Without layers the expressions can only be `"all"` or `"none"`. A driver may
pass the regions (`yelmo_init(..., reg=)`, on the Yelmo grid) instead of
`regions_group`, and `mask_ice` instead of the mask_ice expressions. The named
masks of the regions (`masks`, `mask_<name>`) become Yelmo's regional output
domains. `bnd%regions` holds the codes of the layer `region` (FesmData: the deepest
level), 0 without one.

::: {.callout-tip}
To browse, build, compare and validate parameter files from the command line, see
[`yelmo-config`](yelmo-config.md).
:::
