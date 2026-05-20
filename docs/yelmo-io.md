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

### yelmo_write_init

```fortran
subroutine yelmo_write_init(ylmo,filename,time_init,units,irange,jrange)
```

This routine can be used to initialize any file that will make use of one or more dimension axes of Yelmo variables. The dimension variables that will be written to the file are the following:

```fortran
xc, yc, month, zeta, zeta_ac, zeta_rock, age_iso, pd_age_iso, pc_steps, time [unlimited]
```

Some of the dimension variables above are typically only needed for restart files (`age_iso, pd_age_iso, pc_steps`), but are written as well to maintain generality.

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

This routine will write several variables to a file for a given timestep. The variable names can be provided as a vector of strings via the `nms` argument (e.g., `nms=["H_ice","z_srf"]`). The routine will write relevant model performance information and then individually call `yelmo_write_var` for each variable listed. Optionally it is possible to write comparison fields with present-day data (`compare_pd=.TRUE.`), assuming it has been loaded into the `ylmo%dta` fields.

This routine can also be used to write regional output using the arguements `irange, jrange`.

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
            names(8)  = "beta"
            names(9)  = "visc_bar"
            names(10) = "T_prime_b"
            names(11) = "H_w"
            names(12) = "mb_net"
            names(13) = "smb"
            names(14) = "bmb"
            names(15) = "cmb"
            names(16) = "z_sl"
```

### yelmo_write_restart

```fortran
subroutine yelmo_restart_write(ylmo,filename,time,init,irange,jrange)
```

This routine will save a snapshot of the Yelmo instance. Essentially the routine will loop over every field found in the [Yelmo variable tables](yelmo-variables.md) and write them to a NetCDF file. Optionally `init=.FALSE.` will allow writing of multiple timesteps to the same file (largely useful for diagnostic purposes, since the files can get very large).

This routine can also be used to write regional output using the arguements `irange, jrange`.

## Reading input

By specifying the parameter `yelmo.restart` to a restart file path, Yelmo will read the NetCDF file with a saved snapshot. The routines `yelmo_restart_read_topo_bnd` and `yelmo_restart_read` are generally used internally during `yelmo_init` and `yelmo_init_state`, respectively. So these routines will not typically be needed by a user externally.
