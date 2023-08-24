<!-- ========================================= -->
# Contents
<!-- ========================================= -->

1. [Repository description](#Description)
2. [Compilation](#Compilation)
3. [Usage](#Usage)
4. [Supplementary material](#Supplementary)

<!-- ========================================= -->
# 1. Repository description <a name="Description"></a>
<!-- ========================================= -->

This repository contains the tool to pre-process atmospheric fields for the atmospheric boundary layer model (ABL1D). This tool can read atmospheric fields from ECMWF products (IFS oper, ERA5, ERAI, ...). It allows to compute geostrophic pressure gradients, to interpolate the fields from the IFS vertical grid to the ABL1D vertical grid (which is defined in namelist) and to fill in the field values over land by extrapolating the values over sea.

This repository contains the following files:
* **src**: contains the Makefile and the source codes in Fortran 90 ;
* **ktest**: contains the files necessary to launch the tool and prepare the ABL1d forcing for CROCO BENGUELA_LR case study.

<!-- ========================================= -->
# 2. Compilation <a name="Compilation"></a>
<!-- ========================================= -->
 
To compile the source code, go to the src/ folder and run the command :
```bash
make
```

4 executables are created : 
* get_atm_LSfrc.exe: program to calculate the pressure gradients;
* vinterp_atm_frc.exe : program to vertically interpolate the IFS fields to the ABL1d;
* drown_atm_frc.exe: program to fill in the values on land by extrapolating the values on sea (avoid strong gradients near the coast);
* slp_smoothing.exe : program to test different methods to smooth the surface pressure (Gibbs effect). Using this program is not necessary if you only want to create the forcing for ABL1d.

To clean up the compiled code, use the command :
```bash
make clean
```

<!-- ========================================= -->
# 3. Usage <a name="Usage"></a>
<!-- ========================================= -->

To test the tool, go to the ktest/benguela_ifs_era5 folder. This folder contains all the files and scripts necessary for the
preprocessing tool. It contains the following folders: 
* **0_input_data**: contains the input files for the preprocessing tool. They correspond to the IFS files extracted from the ECMWF servers (in grib format).
* **1_run_preprocessing_tool_N_abl_50**: contains the scripts to run the pre-processing tool previously compiled.

To use the pre-processing tool on files containing only one date per file (0_input_data/ecmwf.20141106.00.nc, ...), the script :
```bash
./run_preprocessing.sh
```

**Note:**
* You must first convert the grib files to netcdf using the script (convert_grib_to_netcdf_with_cdo.sh) located in input_data directory.

To use the tool on a netcdf file containing all the dates (input_data/concatenate_all_ecmwf_files.nc), use the script :
```bash
./run_all_preprocessing.sh
```
**Note:**
* You must first concatenate the netcdf files using the script (concatenate_ecmwf_files.sh.sh) located in input_data.


<!-- ========================================= -->
# 4. Supplementary material <a name="Supplementary"></a>
<!-- ========================================= -->

To test the smoothing of the SLP, you can also test gcm-filters. Here are some commands to install python / gcm-filters if needed (with conda) and test some scripts :
```bash
https://gcm-filters.readthedocs.io/en/latest/how_to_contribute.html
conda create -n gcm-filters-env -c conda-forge --file requirements.txt --file requirements-dev.txt
conda activate gcm-filters-env
conda install ipython
conda install basemap
```

