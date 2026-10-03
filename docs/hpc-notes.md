# HPC Notes

`configme` supports the clusters below as machines (`configme list`) and detects
them from the hostname. Load the modules listed for each cluster first, so that
the compilers and the netCDF tools (`nf-config`, `nc-config`) are found, then
install as usual (see [Installation](getting-started.md)), e.g.

```bash
configme install yelmo -m dkrz_levante -c ifx
```

The input data for realistic domains are on each cluster in `<datapath>/ice_data`;
link them into the checkout with `ln -s <datapath>/ice_data ice_data`. To submit
jobs, set `hpc` and `account` in `.runme/config.toml` and use `runme -rs`.

## Running at PIK on HPC2024 (foote)

The following modules have to be loaded in order to compile and run the model.
For convenience you can also add those commands to your `.profile` file in your home directory.

```bash
module purge
module use /p/system/modulefiles/compiler \
            /p/system/modulefiles/gpu \
            /p/system/modulefiles/libraries \
            /p/system/modulefiles/parallel \
            /p/system/modulefiles/tools

module load intel/oneAPI/2024.0.0
module load netcdf-c/4.9.2
module load netcdf-fortran-intel/4.6.1
module load udunits/2.2.28
module load ncview/2.1.10
module load cdo/2.4.2
```

Install with `configme install yelmo -m pik_hpc2024 -c ifx`.

Data path:

```bash
datapath=/p/projects/megarun
```

## Running at AWI on albedo

Load the following modules in your `.bashrc` or `.bash_profile` file in your home directory.

```bash
module load intel-oneapi-compilers/2024.1.0
module load netcdf-c/4.8.1-openmpi4.1.3-oneapi2022.1.0
module load netcdf-fortran/4.5.4-oneapi2022.1.0
module load udunits/2.2.28
module load ncview/2.1.8
module load cdo/2.5.2
module load nco/5.0.1
module load python/3.11.7
module load git-lfs/3.1.2
```

Install with `configme install yelmo -m awi_albedo -c ifx`.

Data path:

```bash
datapath=/albedo/work/projects/p_forclima
```

## Running at DKRZ on levante

Load the following modules in your `.bashrc` file in your home directory.

```bash
# Tools
module load cdo/2.6.0-gcc-11.2.0
module load esmvaltool/2.5.0
module load ncview/2.1.8-gcc-11.2.0
module load git/2.43.7-gcc-11.2.0
module load python3/2023.01-gcc-11.2.0

# Compilers and libs
module load intel-oneapi-compilers/2023.2.1-gcc-11.2.0
module load netcdf-c/4.8.1-openmpi-4.1.2-intel-2021.5.0
module load netcdf-fortran/4.5.3-openmpi-4.1.2-intel-2021.5.0
```

Install with `configme install yelmo -m dkrz_levante -c ifx`.

Data path:

```bash
datapath=/work/ba1442
```
