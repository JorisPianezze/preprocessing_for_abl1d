<!-- ========================================= -->
# Contents
<!-- ========================================= -->

1. [Repository description](#Description)
2. [Compilation](#Compilation)
3. [Usage](#Usage)
4. [Supplementary material](#Supplementary)
5. [Reference](#Reference)

<!-- ========================================= -->
# 1. Repository description <a name="Description"></a>
<!-- ========================================= -->

This Git's repository contains the tools to pre-process atmospheric fields for the atmospheric boundary layer model (ABL1d ; Lemarié et al. 2021) which has been recently included in NEMO and CROCO oceanic models. This tool can read atmospheric fields from ECMWF products (IFS oper, ERA5, ERAI, ...). It allows to interpolate the fields from the IFS vertical grid to the ABL1d vertical grid (which is defined in namelist), to compute geostrophic pressure gradients, and to fill in the field values over land by extrapolating the values over sea.

This repository contains the following directories :

* **src**: contains the Makefile and the source codes in Fortran 90 ;
* **examples**: contains examples to extract and launch the tools.

<!-- ========================================= -->
# 2. Compilation <a name="Compilation"></a>
<!-- ========================================= -->

To compile the source code, go to the src/ folder and run the command :

```bash
make
```

4 executables are created : 

* **get_atm_LSfrc.exe**: program to calculate the large scale pressure gradients;
* **vinterp_atm_frc.exe** : program to vertically interpolate the IFS fields to the ABL1d;
* **drown_atm_frc.exe**: program to fill in the values on land by extrapolating the values on sea (avoid strong gradients near the coast);
* **slp_smoothing.exe** : program to test different methods to smooth the surface pressure (Gibbs effect). Using this program is not necessary if you only want to create the forcing for ABL1d.

To clean up the compiled code, use the command :

```bash
make clean
```

<!-- ========================================= -->
# 3. Usage <a name="Usage"></a>
<!-- ========================================= -->

To test the tool, go to the examples/benguela_ecmwf_era5/ folder corresponding to a test case for CROCO. This directory contains all the files and scripts necessary to run the preprocessing tool. It contains the following directories: 

* **0_input_data**: contains ane exemple of a script to extract ERA5 data. In this example, ERA5 data are extracted every 6 hours from 6 november 2014 at 00 UTC to 9 november 2014 at 18 UTC ;
* **1_run_preprocessing_tool_N_abl_50**: contains the scripts to run the pre-processing tools previously compiled ;
* **2_convert_to_croco**: contains the scripts to convert pre-processing outputs to a format readable by CROCO.

<!-- ----------------------------------------- -->
## 3.1 Extract ERA5 data 
<!-- ----------------------------------------- -->

First, go to **0_input_data** directory and extract ERA5 fields using :

```bash
python extract_era5_for_abl1d.py
```

Then you need to convert extracted grib files into netcdf using :

```bash
./convert_grib_to_netcdf.sh
```

At the end, you must have a file called era5_for_preprocessing_tools.nc into your **0_input_data** directory. This file contains all necessary informations for pre-proccessing tools.

<!-- ----------------------------------------- -->
## 3.2 Run pre-processing tools
<!-- ----------------------------------------- -->

Then, go to the **1_run_preprocessing_tool_N_abl_50** directory and run the pre-processing tool via :

```bash
./run_preprocessing.sh
```

At the end of pre-proccessing tools executions, you need to have new files called :

* **geos_era5_for_preprocessing_tools.nc**: contains large scale geostrophic winds or pressure

* **out_era5_for_preprocessing_tools.nc**: contains all 3d atmospheric fields without drowning over land

* **out_drown_era5_for_preprocessing_tools.nc**: contains all 3d atmospheric fields with drowning over land

* **vertical_grid_N_abl_50_H_2000m_DZ_10m.nc**: contains vertical grid for ABL1d

**Note**: Only **out_drown_era5_for_preprocessing_tools.nc** and **vertical_grid_N_abl_50_H_2000m_DZ_10m.nc** are used to prepare forcing for NEMO or CROCO oceanic models.

<!-- ----------------------------------------- -->
## 3.3 Convert for CROCO
<!-- ----------------------------------------- -->

Now, go to **2_convert_to_croco** directory to convert preproccess file to a readable format for CROCO using :

```bash
python create_croco_abl_from_ifs_era5.py
```

This script performs an horizontal interpolation of ERA5 fields on CROCO grid and rename variables.

The new file created is called **croco_abl_from_ifs_era5_N_abl_50_interp_linear.nc** and can be used directly to force CROCO using (#undef ONLINE).

<!-- ========================================= -->
# 4. Supplementary material <a name="Supplementary"></a>
<!-- ========================================= -->

To test the smoothing of the Sea Level Pressure, you can also test gcm-filters. Here are some commands to install python / gcm-filters if needed (with conda) and test some scripts :

```bash
https://gcm-filters.readthedocs.io/en/latest/how_to_contribute.html
conda create -n gcm-filters-env -c conda-forge --file requirements.txt --file requirements-dev.txt
conda activate gcm-filters-env
conda install ipython
conda install basemap
```

<!-- ========================================= -->
# 5. Reference <a name="Reference"></a>
<!-- ========================================= -->

Lemarié, F., Samson, G., Redelsperger, J.-L., Giordani, H., Brivoal, T., and Madec, G.: A simplified atmospheric boundary layer model for an improved representation of air–sea interactions in eddying oceanic models: implementation and first evaluation in NEMO (4.0), Geosci. Model Dev., 14, 543–572, https://doi.org/10.5194/gmd-14-543-2021, 2021.
