# HPC Notes

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

When installing `fesm-utils` (see [Dependencies](dependencies.md)) use the `pik` script:

```bash
./install_pik.sh ifx
```

To link to data sources, use the following path:

```bash
datapath=/p/projects/megarun
```

## Running at AWI on albedo

Load the following modules in your `.bashrc` or `.bash_profile` file in your home directory.

```bash
    module load intel-oneapi-compilers/2024.0.0
    module load netcdf-c/4.8.1-openmpi4.1.3-oneapi2022.1.0
    module load netcdf-fortran/4.5.4-oneapi2022.1.0
    module load udunits/2.2.28
    module load ncview/2.1.8
    module load cdo/2.2.0
    module load python/3.10.4
```

When installing `fesm-utils` (see [Dependencies](dependencies.md)) use the `awi` script (which is a link to the `dkrz` script):

```bash
./install_awi.sh ifx
```

To link to data sources, use the following path:

```bash
datapath=/albedo/work/projects/p_forclima
```

## Running at DKRZ on levante

Load the following modules in your `.bashrc` file in your home directory.

```bash
# Tools
module load cdo/2.4.0-gcc-11.2.0
module load esmvaltool/2.5.0
module load ncview/2.1.8-gcc-11.2.0
module load git/2.43.3-gcc-11.2.0
module load python3/2023.01-gcc-11.2.0

# Compilers and libs
module load intel-oneapi-compilers/2023.2.1-gcc-11.2.0
module load netcdf-c/4.8.1-openmpi-4.1.2-intel-2021.5.0
module load netcdf-fortran/4.5.3-openmpi-4.1.2-intel-2021.5.0
```

When installing `fesm-utils` (see [Dependencies](dependencies.md)) use the `dkrz` script:

```bash
./install_dkrz.sh ifx
```

To link to data sources, use the following path:

```bash
datapath=/work/ba1442
```
