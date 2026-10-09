# C API

`src/yelmo_c_api.f90` exposes Yelmo to programs written in other languages
(C, C++, Julia, Python, ...) through C-bound routines (`bind(C)`). A host
program can initialize Yelmo, advance it in time, exchange 2D and 3D fields
with it and write restart files.

## Building

The C API is compiled into a shared library. All objects, including the static
dependencies and `fesm-utils`, must be position-independent, so build with
`pic=1` from a clean tree:

```bash
make clean
make yelmo-c pic=1
```

This builds `libyelmo.a` and then `libyelmo_c_api.so` in `libyelmo/include/`.
`make yelmo-c` without `pic=1` stops with an error. `fesm-utils` must also be
built with `pic=1` in its own checkout.

## Conventions

- **Instances.** The library holds two Yelmo instances. Every routine takes an
  `alias` string as its last argument: `"ylmo2"` selects the second instance,
  any other value the first.
- **Strings** are null-terminated C strings (at most 1028 characters; variable
  names at most 56).
- **Real numbers** are `double`. Yelmo converts them to its working precision.
- **Arrays** are in Fortran (column-major) order: a 2D field has shape
  `nx × ny` with $x$ varying fastest, a 3D field `nx × ny × nz` with the
  vertical index slowest. Staggered fields (e.g. `ux_bar` on the $x$ faces,
  acx-nodes) have the same `nx × ny` shape, with the value of cell $(i,j)$ on
  its right ($x$) or top ($y$) face. The caller allocates all buffers.
- **Integer arguments** are passed by value, except in `yelmo_init_grid` and
  `yelmo_get_grid_sizes`, which take pointers.
- **Units** are those of the [variable tables](yelmo-variables.md).

## Routines

C prototypes of all exported routines:

```c
/* Initialization */
void yelmo_init_grid(const char *grid_name, const int *nx, const int *ny,
                     const double *xc, const double *yc,
                     const double *lon, const double *lat, const double *area,
                     const char *alias);
void yelmo_init(const char *filename, const char *grid_def, double time,
                const char *alias);
void yelmo_init_state(double time, const char *thrm_method, const char *alias);

/* Time stepping and output */
void yelmo_step(double time, const char *alias);
void yelmo_restart_write(const char *filename, double time, const char *alias);

/* Grid */
void yelmo_get_grid_sizes(int *nx, int *ny, int *nz_aa, int *nz_ac,
                          int *nzr_aa, int *nzr_ac, const char *alias);
void yelmo_get_grid_info(double *xc, double *yc, double *zeta_aa, double *zeta_ac,
                         double *zeta_r_aa, double *zeta_r_ac, const char *alias);

/* Field exchange */
void yelmo_get_var2D(double *v, int nx, int ny, const char *name, const char *alias);
void yelmo_get_var3D(double *v, int nx, int ny, int nz, const char *name, const char *alias);
void yelmo_set_var2D(const double *v, int nx, int ny, const char *name, const char *alias);
void yelmo_set_var3D(const double *v, int nx, int ny, int nz, const char *name, const char *alias);

/* Effective-pressure callback */
typedef void (*yelmo_neff_cb)(int tag, const double *uxy_b, double *N_eff,
                              int nx, int ny);
void yelmo_set_neff_callback(yelmo_neff_cb cb, int tag, const char *alias);

/* Used internally by every routine above */
void yelmo_set_alias(const char *alias);
```

### Initialization

`yelmo_init(filename, grid_def, time, alias)` calls `yelmo_init` with the
parameter file `filename` (see [Parameters](parameters.md)), so the domain and
grid name are taken from the `&yelmo` group of that file. `grid_def` selects
how the grid is defined:

- `"file"`: read from `yelmo.grid_path`;
- `"name"`: generated from `yelmo.grid_name`;
- `"none"`: already defined by `yelmo_init_grid`.

`yelmo_init_grid(grid_name, nx, ny, xc, yc, lon, lat, area, alias)` defines
the grid from its axes `xc[nx]` and `yc[ny]` and the 2D fields `lon`, `lat`
and `area` (`nx × ny`). Call it before `yelmo_init` with `grid_def = "none"`.

`yelmo_init_state(time, thrm_method, alias)` initializes the model state
(`yelmo_init_state`), from the restart file if `yelmo.restart` is set,
otherwise with the temperature profile `thrm_method` (`"linear"`, `"robin"` or
`"robin-cold"`; see [Thermodynamics](physics/thermodynamics.md)). The
`"prescribed"` option needs a temperature field and is not available through
the C API. Set the boundary fields with `yelmo_set_var2D` before this call.

### Time stepping and output

`yelmo_step(time, alias)` advances the model to `time` [a] (`yelmo_update`,
see [Time stepping](physics/timestepping.md)). Set the boundary fields with
`yelmo_set_var2D` before each call.

`yelmo_restart_write(filename, time, alias)` writes a restart file
(`yelmo_restart_write`, see [Input and output](yelmo-io.md)).

There is no binding for `yelmo_end`: the instances live until the library is
unloaded.

### Grid

`yelmo_get_grid_sizes` returns the numbers of grid points: `nx`, `ny`, the ice
levels `nz_aa` (layer centres) and `nz_ac` (layer interfaces), and the bedrock
levels `nzr_aa` and `nzr_ac`. `yelmo_get_grid_info` then fills buffers of
these sizes with the axes `xc`, `yc` and the normalized heights `zeta_aa`,
`zeta_ac` (ice) and `zeta_r_aa`, `zeta_r_ac` (bedrock). A null pointer skips
an axis.

### Field exchange

`yelmo_get_var2D` / `yelmo_get_var3D` copy a field into the caller's buffer;
`yelmo_set_var2D` / `yelmo_set_var3D` copy the buffer into the model state.
Fields are addressed by a name made of a component prefix and the field name,
e.g. `"tpo_H_ice"`. The buffer dimensions must match the field: they are not
checked. 2D fields are `nx × ny`. 3D fields have `nz` = `nz_aa`, except
`dyn_uz` and `dyn_uz_star` (`nz_ac`), `thrm_T_rock` (`nzr_aa`) and
`trc_depth_iso` (the number of isochrones, `ytrc.time_iso`).

A getter called with an unknown name fills the buffer with the missing value
−9999. A setter called with an unknown name prints a message and leaves the
model unchanged. Setters only copy the field; nothing is recomputed until the
next `yelmo_step`. Integer and logical fields are exchanged as doubles: the
setter rounds `bnd_mask_ice` to the nearest integer, and `bnd_calv_mask` is
true where the value is non-zero. `bnd_smb_ref` is the surface mass balance
`bnd%smb`.

**2D getters** (`yelmo_get_var2D`):

| Prefix | Object | Names (after the prefix) |
|---|---|---|
| `bnd_` | `bnd` | `z_bed`, `z_bed_sd`, `z_sl`, `H_sed`, `smb_ref`, `T_srf`, `bmb_shlf`, `fmb_shlf`, `T_shlf`, `Q_geo`, `enh_srf`, `basins`, `basin_mask`, `regions`, `H_ice_ref`, `z_bed_ref`, `calv_mask`, `tau_relax`, `z_bed_corr`, `dzbdt_corr`, `mask_ice` |
| `dta_` | `dta%pd` | `pd_uxy_s`, `pd_H_grnd`, `pd_H_ice`, `pd_z_srf`, `pd_mask_bed` |
| `tpo_` | `tpo%now` | `H_ice`, `dHidt`, `dHidt_dyn`, `mb_net`, `mb_relax`, `mb_resid`, `mb_err`, `smb`, `bmb`, `fmb`, `dmb`, `cmb`, `bmb_ref`, `fmb_ref`, `dmb_ref`, `cmb_flt`, `cmb_flt_x`, `cmb_flt_y`, `cmb_grnd`, `cmb_grnd_x`, `cmb_grnd_y`, `cr_acx`, `cr_acy`, `calv_rate_flt`, `calv_rate_grnd`, `lsf`, `dlsfdt`, `z_srf`, `dzsdt`, `dzsdt_kin`, `dzbdt_kin`, `dHidt_vert`, `mask_kin`, `eps_eff`, `tau_eff`, `z_base`, `dzsdx`, `dzsdy`, `dHidx`, `dHidy`, `dzbdx`, `dzbdy`, `dzsdx_aa`, `dzsdy_aa`, `dHidx_aa`, `dHidy_aa`, `dzbdx_aa`, `dzbdy_aa`, `H_eff`, `H_grnd`, `H_calv`, `kt`, `z_bed_filt`, `f_grnd`, `f_grnd_acx`, `f_grnd_acy`, `f_grnd_ab`, `f_ice`, `f_grnd_bmb`, `f_grnd_pin`, `dist_margin`, `dist_grline`, `dHidt_dyn_raw_n`, `H_ice_n`, `z_srf_n`, `lsf_n`, `H_ice_dyn`, `f_ice_dyn`, `tau_relax`, `mask_bed`, `mask_grz`, `mask_frnt` |
| `dyn_` | `dyn%now` | `ux_bar`, `uy_bar`, `uxy_bar`, `ux_bar_prev`, `uy_bar_prev`, `ux_b`, `uy_b`, `uz_b`, `uxy_b`, `ux_s`, `uy_s`, `uz_s`, `uxy_s`, `ux_i_bar`, `uy_i_bar`, `uxy_i_bar`, `duxydt`, `duxdz_bar`, `duydz_bar`, `taud_acx`, `taud_acy`, `taud`, `taub_acx`, `taub_acy`, `taub`, `taul_int_acx`, `taul_int_acy`, `qq_gl_acx`, `qq_gl_acy`, `qq_acx`, `qq_acy`, `qq`, `visc_eff_int`, `N_eff`, `cb_tgt`, `cb_ref`, `c_bed`, `f_slide`, `beta_acx`, `beta_acy`, `beta`, `beta_eff`, `f_vbvs`, `ssa_err_acx`, `ssa_err_acy`, `ssa_mask_acx`, `ssa_mask_acy`, `strn2D_dxx`, `strn2D_dyy`, `strn2D_dxy`, `strn2D_dxz`, `strn2D_dyz`, `strn2D_de`, `strn2D_div`, `strn2D_f_shear`, `strn2D_eps_eig_1`, `strn2D_eps_eig_2` |
| `mat_` | `mat%now` | `enh_bar`, `ATT_bar`, `visc_bar`, `visc_int`, `f_shear_bar`, `strn2D_dxx`, `strn2D_dyy`, `strn2D_dxy`, `strn2D_dxz`, `strn2D_dyz`, `strn2D_de`, `strn2D_div`, `strn2D_f_shear`, `strn2D_eps_eig_1`, `strn2D_eps_eig_2`, `strs2D_txx`, `strs2D_tyy`, `strs2D_txy`, `strs2D_txz`, `strs2D_tyz`, `strs2D_te`, `strs2D_tau_eig_1`, `strs2D_tau_eig_2` |
| `thrm_` | `thrm%now` | `f_pmp`, `bmb_grnd`, `Q_b`, `Q_ice_b`, `T_prime_b`, `H_cts`, `bmb_grnd_star`, `bc_b`, `bmb_clamp`, `melt_int`, `Q_rock` |
| `hyd_` | `hyd%now` | `C_frz`, `Q_diss`, `Q_sens` |

The `dta_pd_` fields are present-day data, read-only.

**3D getters** (`yelmo_get_var3D`):

| Prefix | Object | Names (after the prefix) |
|---|---|---|
| `dyn_` | `dyn%now` | `ux`, `uy`, `uxy`, `uz`, `uz_star`, `ux_i`, `uy_i`, `duxdz`, `duydz`, `de_eff`, `visc_eff`, `strn_dxx`, `strn_dyy`, `strn_dxy`, `strn_dxz`, `strn_dyz`, `strn_de`, `strn_div`, `strn_f_shear`, `jvel_dxx`, `jvel_dxy`, `jvel_dxz`, `jvel_dyx`, `jvel_dyy`, `jvel_dyz`, `jvel_dzx`, `jvel_dzy`, `jvel_dzz` |
| `mat_` | `mat%now` | `enh`, `enh_bnd`, `ATT`, `visc`, `strn_dxx`, `strn_dyy`, `strn_dxy`, `strn_dxz`, `strn_dyz`, `strn_de`, `strn_div`, `strn_f_shear`, `strs_txx`, `strs_tyy`, `strs_txy`, `strs_txz`, `strs_tyz`, `strs_te` |
| `trc_` | `trc%now` | `t_dep`, `depth_iso` |
| `thrm_` | `thrm%now` | `enth`, `T_ice`, `omega`, `T_pmp`, `T_prime`, `Q_strn`, `dQsdT`, `cp`, `kt`, `advecxy`, `T_rock` |

**2D setters** (`yelmo_set_var2D`):

| Prefix | Object | Names (after the prefix) |
|---|---|---|
| `bnd_` | `bnd` | `z_bed`, `z_bed_sd`, `z_sl`, `H_sed`, `smb_ref`, `T_srf`, `bmb_shlf`, `fmb_shlf`, `T_shlf`, `Q_geo`, `enh_srf`, `basins`, `basin_mask`, `regions`, `H_ice_ref`, `z_bed_ref`, `calv_mask`, `tau_relax`, `z_bed_corr`, `dzbdt_corr`, `mask_ice` |
| `tpo_` | `tpo%now` | `H_ice` |
| `dyn_` | `dyn%now` | `N_eff`, `cb_tgt`, `cb_ref`, `c_bed` |
| `hyd_` | `hyd%now` | `N`, `W_til`, `C_frz`, `Q_diss`, `Q_sens` |

The `hyd_` setters let a host model drive the basal hydrology, e.g. with
`yhyd.bkt_N_closure = -1` for an externally set effective pressure `hyd_N`.

**3D setters** (`yelmo_set_var3D`):

| Prefix | Object | Names (after the prefix) |
|---|---|---|
| `dyn_` | `dyn%now` | `ux`, `uy`, `uz` |
| `thrm_` | `thrm%now` | `T_ice` |

### Effective-pressure callback

`yelmo_set_neff_callback(cb, tag, alias)` registers a host function that
computes the effective pressure from the basal velocity. With
`ydyn.solver = "diva"` or `"diva-noslip"`, Yelmo calls it in every Picard
iteration of the [velocity solution](physics/momentum/diva.md#picard-iteration)
as `cb(tag, uxy_b, N_eff, nx, ny)`, with

- `tag`: the integer passed at registration, handed back unchanged (e.g. to
  identify the host object);
- `uxy_b`: the basal speed [m a$^{-1}$] on the cell centres (`nx × ny`, input);
- `N_eff`: the effective pressure [Pa] (`nx × ny`), holding the current value
  on input and the new value on output.

Yelmo stores `N_eff` as the hydrology's effective pressure, recomputes the
friction coefficient `c_bed` from it, and requires `c_bed` to converge
together with the velocity. A steady hydrology owned by the host is thus
solved together with the sliding speed instead of lagging it by one step.
Set `yhyd.bkt_N_closure = -1` (externally set effective pressure). A null
`cb` unregisters the callback. The other momentum solvers do not call it.

## Example

```c
const char *a = "ylmo1";
int nx, ny, nz_aa, nz_ac, nzr_aa, nzr_ac;

yelmo_init("par/yelmo_initmip_grl.nml", "file", 0.0, a);
yelmo_get_grid_sizes(&nx, &ny, &nz_aa, &nz_ac, &nzr_aa, &nzr_ac, a);

double *smb  = malloc(sizeof(double) * nx * ny);
double *H    = malloc(sizeof(double) * nx * ny);
/* ... fill smb and the other boundary fields ... */
yelmo_set_var2D(smb, nx, ny, "bnd_smb_ref", a);

yelmo_init_state(0.0, "robin-cold", a);

for (double t = 10.0; t <= 1000.0; t += 10.0) {
    yelmo_set_var2D(smb, nx, ny, "bnd_smb_ref", a);
    yelmo_step(t, a);
}
yelmo_get_var2D(H, nx, ny, "tpo_H_ice", a);
yelmo_restart_write("yelmo_restart.nc", 1000.0, a);
```
