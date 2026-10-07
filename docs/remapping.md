# Remapping

Yelmo runs on a Cartesian (x/y) grid. Often input data comes in many formats, global lat/lon grids, projections and sets of points. It is important to have robust remapping tools.

Typically for a given domain, we define a Polar Stereographic projection to be able to convert lat/lon data points onto a Cartesian plane. For Antarctica, for example, the standard projection has the following parameters:

```bash
int polar_stereographic ;
    polar_stereographic:grid_mapping_name = "polar_stereographic" ;
    polar_stereographic:straight_vertical_longitude_from_pole = 0. ;
    polar_stereographic:latitude_of_projection_origin = -90. ;
    polar_stereographic:standard_parallel = -71. ;
    polar_stereographic:false_easting = 0. ;
    polar_stereographic:false_northing = 0. ;
```

## Naming files

For grids used by Yelmo, we generally use an abbreviation for the domain name followed by the resolution. So for Antarctica, we could have the grids `ANT-32KM` or `ANT-16KM` for a 32km or 16km grid, respectively. Data that have been projected onto these grids are saved with the grid name as a prefix followed by a general name that specifies the type of data, e.g., `CLIM` or `TOPO`, finally followed by more descriptive information about the specific dataset `IPSL-14Ma` or `IPSL-PD-CTRL`. For example, the RTopo-2.0.1 topography dataset processed onto the 32 km grid is called `ANT-32KM_TOPO-RTOPO-2.0.1.nc`.

## Fields Yelmo needs

To drive Yelmo with boundary conditions derived from a climate model, it needs the following fields to be defined on the Polar Stereographic grid:

- Climatological mean near-surface air temperature [monthly]
- Climatological mean precipitation [monthly]
- Surface elevation
- Sea level
- Ice thickness

- Climatological mean 3D ocean temperature [annual]
- Climatological mean 3D ocean salinity [annual]
- Oceanic bathymetry

Likely these would be processed into two or more separate files, e.g., one for climate `CLIM` variables and another for ocean `OCN` variables.

## Preprocessing data using `cdo`

As a first step, the Climate Data Operators `cdo` package is great for most preprocessing steps. It can handle averaging data over time and space, merging data files, extracting individual variables etc. See the extensive documentation and examples online.

For example, it is possible to use the command `cdo selvar` to extract specific variables from a file:

```bash
cdo selvar,t2m,precip diane_C14Ma_1_5PAL_SE_4750_4849_1M_histmth.nc ipsl_tmp1.nc
```

If you have several variables in individual files, you can then conveniently merge them into one file using `merge` (it's better if they have the same shape):

```bash
# Extract t2m to a temporary file
cdo selvar,t2m diane_C14Ma_1_5PAL_SE_4750_4849_1M_histmth.nc ipsl_tmp1.nc

# Extract precip to a temporary file
cdo selvar,precip diane_C14Ma_1_5PAL_SE_4750_4849_1M_histmth.nc ipsl_tmp2.nc

# Merge the two individual variable files into one convenient file
cdo merge ipsl_tmp1.nc ipsl_tmp2.nc ipsl_tmp3.nc
```

There are many other useful commands, particularly for getting monthly means `cdo monmean ...` and other statistics.

Resources:

CDO Documentation page:
[https://code.mpimet.mpg.de/projects/cdo/wiki/Cdo#Documentation](https://code.mpimet.mpg.de/projects/cdo/wiki/Cdo#Documentation)

CDO User guide:
[https://code.mpimet.mpg.de/projects/cdo/embedded/cdo.pdf](https://code.mpimet.mpg.de/projects/cdo/embedded/cdo.pdf)

CDO Reference card:
[https://code.mpimet.mpg.de/projects/cdo/embedded/cdo_refcard.pdf](https://code.mpimet.mpg.de/projects/cdo/embedded/cdo_refcard.pdf)

## Using `cdo` for remapping

To remap a data file from lat/lon coordinates to our projection, `cdo` needs a grid description file that describes the target Polar Stereographic projection grid. The grid description files of the standard Yelmo grids are in `maps/` (`maps/gengriddes.sh` generates one from a grid file). For example, for the 32 km Antarctic domain, `maps/grid_ANT-32KM.txt` is:

```bash
gridtype = projection
gridsize =      36481
xsize    =        191
ysize    =        191
xname    = xc
xunits   = km
yname    = yc
yunits   = km
xfirst   =    -3040.000000
xinc     =       32.000000
yfirst   =    -3040.000000
yinc     =       32.000000
grid_mapping = crs
grid_mapping_name = polar_stereographic
straight_vertical_longitude_from_pole =        0.000
latitude_of_projection_origin =      -90.000
standard_parallel =      -71.000
false_easting =        0.000
false_northing =        0.000
semi_major_axis =        6378137.000
inverse_flattening =         298.25722356
```

With this file defined, it's easy to perform projections using the `cdo remap*` commands. To perform a bicubic interpolation, call:

```bash
cdo remapbic,grid_ANT-32KM.txt diane_C14Ma_1_5PAL_SE_4750_4849_1M_histmth.nc ANT-32KM_test-bic.nc
```

Here, `remapbic` specifies bicubic interpolation and `grid_ANT-32KM.txt` defines the target grid as above. Then the source dataset is specified and the desired output file `ANT-32KM_test-bic.nc`.

To perform conservative interpolation, replace `remapbic` with `remapcon`:

```bash
cdo remapcon,grid_ANT-32KM.txt diane_C14Ma_1_5PAL_SE_4750_4849_1M_histmth.nc ANT-32KM_test-con.nc
```

Conservative interpolation is generally preferred, especially when going from a high resolution to a lower resolution, as it avoids unwanted interpolation artifacts and conserves the quantity being remapped. However, from low resolution to high resolution, conservative interpolation can result in more "blocky" fields with abrupt changes in values. Thus, in this case, bicubic interpolation, or conservative interpolation with additional Gaussian smoothing is better. The latter is not supported by `cdo`, but can be achieved with other tools.

One option for processing may be a conservative remapping, following by a smoothing step:

```bash
cdo remapcon,grid_ANT-32KM.txt diane_C14Ma_1_5PAL_SE_4750_4849_1M_histmth.nc ANT-32KM_test-con.nc
cdo smooth,radius=128km ANT-32KM_test-con.nc ANT-32KM_test-con-smooth.nc

```

The smoothing radius should be chosen such that it is the smallest value possible that removes blocky artifacts from the field.

## Summary

It can be tedious to process data from a climate model into the right format to drive Yelmo. Tools like `cdo` help to reduce this burden. Other tools like NetCDF Operator `NCO` and today numerous Python-based libraries and tools can also be used.

It is best to define a script or program with all the processing steps clearly defined. That way, when new data becomes available from the same model, it is easy to process it systematically (and reproducibly) in the same way without any trouble.

## Remapping restart file

Sometimes we may want to restart a simulation at a new resolution, e.g. perform
a spin-up at relatively low resolution and then continue at higher resolution.
Yelmo interpolates a restart file from another grid itself: when the
`grid_name` attribute of the restart file differs from `yelmo.grid_name`, it
computes conservative remapping weights between the two grids and interpolates
the restart state (see [Yelmo IO](yelmo-io.md#reading-input)). No separate
remapping step is needed.

For example, run a short 32 km Greenland simulation, which writes a restart file:

```bash
runme -r -e initmip -n par/yelmo_initmip.nml -o output/restarts/sim0-32km -p ctrl.time_end=100 ctrl.set_nm="set_grl_pd" yelmo.domain="Greenland" yelmo.grid_name="GRL-32KM"
```

and continue it at 16 km from this restart file:

```bash
runme -r -e initmip -n par/yelmo_initmip.nml -o output/restarts/sim1-16km -p ctrl.time_end=100 ctrl.set_nm="set_grl_pd" yelmo.domain="Greenland" yelmo.grid_name="GRL-16KM" yelmo.restart="../sim0-32km/yelmo_restart.nc"
```

To remap other fields between grids with `cdo`, `maps/genmap.sh` generates
conservative remapping weights from the grid description files, e.g.:

```bash
cdo gencon,grid_GRL-16KM.txt -setgrid,grid_GRL-32KM.txt GRL-32KM_REGIONS.nc scrip-con_GRL-32KM_GRL-16KM.nc
```

## Remapping fields inside a program

The module `yelmo_remapping` (re-exported by `use yelmo`) maps 2D and 3D fields
from another source, such as another model or another resolution, onto the
Yelmo grid:

- `yelmo_remap(var, var_in, mp, [name])` maps a 2D field with a coords map `mp`.
- `yelmo_remap(var, zeta, var_in, zeta_in, [mp], [name])` maps a 3D field (`name`:
  field name for messages, also in 2D). Each column
  is first interpolated linearly from the source levels `zeta_in` onto the Yelmo
  levels `zeta` (constant beyond the end levels), and then each level is
  remapped horizontally if `mp` is given. Levels are normalized heights (0 at
  the base, 1 at the top), as `zeta_aa`, `zeta_ac` and the bedrock axes.
- `yelmo_load_map(mp, grd, filename, src_grid_name, [method], [gen])` builds
  the map from the xc/yc axes of a NetCDF file onto the Yelmo grid, with the
  projection of the Yelmo grid (method `"con"` by default, as for restart
  files). With `gen = "cdo"`, a SCRIP map generated with `cdo` is loaded from
  `maps/` instead (default `"coords"`: weights computed by Yelmo).
- `yelmo_read_remap(var, grd, filename, varname, [zeta, zeta_name], [method])`
  reads a 2D or 3D field (the last time record, if there is a time dimension)
  and maps it onto the Yelmo grid. The field is remapped horizontally only if
  the file's axes differ from the Yelmo grid.

For example, to initialize the ice temperature from a field `T_ice` on the
levels `zeta` of another file:

```fortran
allocate(T_ice(nx,ny,size(yelmo1%thrm%par%z%zeta_aa)))
call yelmo_read_remap(T_ice,yelmo1%grd,"T_init.nc","T_ice",yelmo1%thrm%par%z%zeta_aa,"zeta")
call yelmo_init_state(yelmo1,time=time_init,thrm_method="prescribed",T_ice=T_ice)
```

With `thrm_method = "prescribed"`, Yelmo caps the temperature at the pressure
melting point, sets the water content to zero and computes the consistent
enthalpy.
