# Remapping

Yelmo runs on a Cartesian (x/y) grid. Often input data comes in many formats, global lat/lon grids, projections and sets of points. It is important to have robust remapping tools.

Typically for a given domain, we define a Polar Stereographic projection to be able to convert lat/lon data points onto a Cartesian plane. For Antarctica, for example, the standard projection has the following parameters:

```bash
int polar_stereographic ;
    polar_stereographic:grid_mapping_name = "polar_stereographic" ;
    polar_stereographic:straight_vertical_longitude_from_pole = 0. ;
    polar_stereographic:latitude_of_projection_origin = -71. ;
    polar_stereographic:angle_of_oblique_tangent = 19. ;
    polar_stereographic:scale_factor_at_projection_origin = 1. ;
    polar_stereographic:false_easting = 0. ;
    polar_stereographic:false_northing = 0. ;
```

## Naming files

For grids used by Yelmo, we generally use an abbreviation for the domain name followed by the resolution. So for Antarctica, we could have the grids `ANT-32KM` or `ANT-16KM` for a 32km or 16km grid, respectively. Data that have been projected onto these grids are saved with the grid name as a prefix followed by a general name that specifies the type of data, e.g., `CLIM` or `TOPO`, finally followed by more descriptive information about the specific dataset `IPSL-14Ma` or `IPSL-PD-CTRL`. For example, the latest topopgraphy dataset we use is called the RTopo2.0.1 dataset, so this is processed into a file called `ANT-32KM_TOPO-RTOPO-2.0.1.nc`.

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

If you have several variables in individual files, you can then conveniently merge them into one file usine `merge` (it's better if they have the same shape):

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

To remap a data file from lat/lon coordinates to our projection, `cdo` needs a grid description file that describes the target Polar Stereographic projection grid. For example, for a 32km resolution domain, we would use the following file named `grid_ANT-32KM.txt`:

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

Here, `remapbic` specifies bicubic interpolation and `grid_ANT-32KM.txt` defines the target grid as above. Then the source dataset is specified and the desired output file `ANT-32KM_test.nc`.

To perform conservative interpolation, replace `remapbic` with `remapcon`:

```bash
cdo remapcon,grid_ANT-32KM.txt diane_C14Ma_1_5PAL_SE_4750_4849_1M_histmth.nc ANT-32KM_test-con.nc
```

Conservative interpolation is generally preferred, especially when going from a high resolution to a lower resolution, as it avoids unwanted interpolation artifacts and conserves the quantity being remapped. However, from low resolution to high resolution, conservative interpolation can result in more "blocky" fields with abrupt changes in values. Thus, in this case, bicubic interpolation, or conservative interpolation with additional Gaussian smoothing is better. The latter is not supported by `cdo`, but can be acheived with other tools.

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

Sometimes we may want to restart a simulation at a new resolution - i.e., perform a spinup simulation at relatively low resolution and then continue the simulation at higher resolution.

1. Use `cdo` to remap the restart file based on the grid definition files.

```bash
# Define env variables as shortcuts to locations of grid files
grid_src=/Users/robinson/models/EURICE/gridding/maps/grid_GRL-32KM.txt
grid_tgt=/Users/robinson/models/EURICE/gridding/maps/grid_GRL-16KM.txt

# Call remapping
cdo remapcon,${grid_tgt} -setgrid,${grid_src} yelmo_restart.nc yelmo_restart_16km.nc
```

Let's do a test. First, run a short 32km Greenland simulation and generate a restart file:

```bash
./runme -r -e initmip -n par/yelmo_initmip.nml -o output/restarts/sim0-32km -p ctrl.time_end=100 ctrl.time_equil=0 ctrl.clim_nm="clim_pd_grl" yelmo.domain="Greenland" yelmo.grid_name="GRL-32KM"
```

That simulation should have produced a nice restart file. Let's test a normal 32km simulation that continues from this restart file.

```bash
./runme -r -e initmip -n par/yelmo_initmip.nml -o output/restarts/sim1-32km -p ctrl.time_end=100 ctrl.time_equil=0 ctrl.clim_nm="clim_pd_grl" yelmo.domain="Greenland" yelmo.grid_name="GRL-32KM" yelmo.restart="../sim0-32km/yelmo_restart.nc"
```

Ok, now generate scrip map file to interpolate from 32km down to 16km.

```bash
domain=Greenland
grid_name_src=GRL-32KM
grid_name_tgt=GRL-4KM
nc_src=../ice_data/${domain}/${grid_name_src}/${grid_name_src}_REGIONS.nc 

cdo gencon,grid_${grid_name_tgt}.txt -setgrid,grid_${grid_name_src}.txt ${nc_src} scrip-con_${grid_name_src}_${grid_name_tgt}.nc

```

Now let's try to run a simulation at 16km, loading the restart file from 32km

```bash
./runme -r -e initmip -n par/yelmo_initmip.nml -o output/restarts/sim2-16km -p ctrl.time_end=100 ctrl.time_equil=0 ctrl.clim_nm="clim_pd_grl" yelmo.domain="Greenland" yelmo.grid_name="GRL-16KM" yelmo.restart="../sim0-32km/yelmo_restart.nc"
```

The simulation is successful! (as of branch `alex-dev-2`, revision `1d9783fb`).
