# Yelmo IO

## Writing output

Multiple generalized routines are available for writing the variables of a Yelmo instance (yelmo_class)
to a NetCDF file. The main public facing routines are the following:

```fortran
yelmo_write_init
yelmo_write_var
yelmo_write_step
yelmo_restart_write
```

These routines will be described briefly below.

### Standard file names

The Yelmo test programs and yelmox use the same file names as CLIMBER-X:

| File | Contents |
|---|---|
| `yelmo.nc` | Main output: 2D and 3D fields at the output interval |
| `yelmo_sm.nc` | Small output: a reduced set of 2D fields, written more often (initmip, yelmox) |
| `yelmo_ts.nc` | Time series of global diagnostics (`yelmo_regions_write`) |
| `yelmo_ts_<region>.nc` | Time series for sub-region `<region>` |
| `yelmo_restart.nc` | Restart snapshot (`yelmo_restart_write`) |
| `yelmo_killed.nc` | Snapshot written by `yelmo_check_kill` before the model is stopped |
| `yelmo_metrics.nc` | Numerical and speed metrics (`yelmo.write_metrics = True`) |
| `timesteps.nc` | Time steps and solver diagnostics (`yelmo.log_timestep = True`) |

The time-series names are the defaults set by `yelmo_region_init`. The killed,
metrics and timestep files are written by Yelmo in its output folder. For the
other files, the calling program chooses the name.

### yelmo_write_init

```fortran
subroutine yelmo_write_init(ylmo,filename,time_init,units,irange,jrange)
```

This routine can be used to initialize any file that will make use of one or more dimension axes of Yelmo variables. The dimension variables that will be written to the file are the following:

```fortran
xc, yc, month, zeta, zeta_ac, zeta_rock, time_iso, pd_time_iso, pc_steps, time [unlimited]
```

plus `depth_norm` when the tracer backend computes statistics. Some of the
dimension variables above are typically only needed for restart files
(`time_iso, pd_time_iso, pc_steps`), but are written as well to maintain
generality. The routine also writes the grid definition and the static fields
`mask_ice`, `basins`, `regions`, `z_bed_sd` and `H_sed`.

Importantly, `yelmo_write_init` can be used to initialize a regional output file by specifying the indices of the bounding box for the region of interest via the arguments `irange=[i1,i2], jrange=[j1,j2]`.

### yelmo_write_var

```fortran
subroutine yelmo_write_var(filename,varname,ylmo,n,ncid,irange,jrange)
```

This routine will write a variable to a given `filename` of an already existing NetCDF file, most likely but not necessarily initialized using `yelmo_write_init`. This routine will accept any variable `varname` that is listed in the [Yelmo variable tables](yelmo-variables.md), which will be written with the attributes specified in the table.

This routine can also be used to write regional output using the arguments `irange, jrange`.

### yelmo_write_step

```fortran
subroutine yelmo_write_step(ylmo,filename,time,nms,compare_pd,irange,jrange)
```

This routine will write several variables to a file for a given timestep. The variable names can be provided as a vector of strings via the `nms` argument (e.g., `nms=["H_ice","z_srf"]`). The routine calls `yelmo_write_var` for each variable listed. Optionally it is possible to write comparison fields with present-day data (`compare_pd=.TRUE.`), assuming it has been loaded into the `ylmo%dta` fields.

This routine can also be used to write regional output using the arguments `irange, jrange`.

Note that this routine can be challenging to use in Fortran, when custom variable names (`nms` argument) is used. This is because of the Fortran limitation on defining string arrays as inline arguments - namely, all strings in the array are required to have the same length.

Passing this argument would give an error:

```fortran
nms=["H_ice","z_srf","mask_bed"]
```

while this would be ok:

```fortran
nms=["H_ice   ","z_srf   ","mask_bed"]
```

For three variables this is not so cumbersome, but can be when many variables are listed.

If no argument is used, then a subset of useful variables is written:

```fortran
            names(1)  = "H_ice"
            names(2)  = "z_srf"
            names(3)  = "z_bed"
            names(4)  = "mask_bed"
            names(5)  = "uxy_b"
            names(6)  = "uxy_s"
            names(7)  = "uxy_bar"
            names(8)  = "ux_bar"
            names(9)  = "uy_bar"
            names(10) = "cb_ref"
            names(11) = "N_eff"
            names(12) = "beta"
            names(13) = "taub"
            names(14) = "taud"
            names(15) = "visc_bar"
            names(16) = "T_prime_b"
            names(17) = "hyd_W_til"
            names(18) = "mb_net"
            names(19) = "smb"
            names(20) = "bmb"
            names(21) = "cmb"
            names(22) = "z_sl"
```

### yelmo_restart_write

```fortran
subroutine yelmo_restart_write(ylmo,filename,time,init,irange,jrange)
```

This routine will save a snapshot of the Yelmo instance. Essentially the routine will loop over every field found in the [Yelmo variable tables](yelmo-variables.md) and write them to a NetCDF file. The tables are read at `yelmo_init` from `input/yelmo-variables-*.md`, so a field that is missing from a table is not written to the restart file. Optionally `init=.FALSE.` will allow writing of multiple timesteps to the same file (largely useful for diagnostic purposes, since the files can get very large).

This routine can also be used to write regional output using the arguments `irange, jrange`.

## Reading input

By specifying the parameter `yelmo.restart` to a restart file path, Yelmo will read the NetCDF file with a saved snapshot. The routines `yelmo_restart_read_topo_bnd` and `yelmo_restart_read` are generally used internally during `yelmo_init` and `yelmo_init_state`, respectively. So these routines will not typically be needed by a user externally.

If the `grid_name` attribute of the restart file differs from `yelmo.grid_name`,
the restart is interpolated to the model grid with conservative remapping
weights, which Yelmo computes itself. A simulation can therefore be continued
at another resolution by setting `yelmo.restart` to a restart file from the
other grid.

Two parameters set which topography is used:

- `yelmo.restart_H_ice` (default `False`): take the ice thickness from the
  restart file; otherwise it comes from the input topography file.
- `yelmo.restart_z_bed` (default `False`): use the bedrock elevation of the
  restart file. Otherwise, or when the restart comes from a coarser grid, the target
  bedrock is the present-day reference bedrock plus the isostatic displacement
  of the restart state. The model starts on the restart bedrock, and the
  difference to the target (e.g. the high-resolution detail missing in a coarse
  restart) is provided as the rate `bnd%dzbdt_corr` over `yelmo.restart_relax`
  years (default 1000), to be applied by the driver's isostasy model. With
  `restart_relax = 0`, the target bedrock is used at once.
