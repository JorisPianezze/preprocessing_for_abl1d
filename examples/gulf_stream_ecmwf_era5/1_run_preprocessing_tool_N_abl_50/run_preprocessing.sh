#!/bin/bash
#######################################################
#SBATCH -J prepabl 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 01:00:00      
#SBATCH -p normal256
#SBATCH  --exclusive
#SBATCH  --no-requeue
#######################################################

ulimit -s unlimited
ulimit -c 0

# ----------------------------------------
#   Load environment
# ----------------------------------------
if [ -e /home/cnrm_other/ge/erla/pianezzej/SAVE/env/env_tools.sh ]
then
  . /home/cnrm_other/ge/erla/pianezzej/SAVE/env/env_tools.sh
fi

# ----------------------------------------
#   Link executables
# ----------------------------------------
export dir_preprocessing_tools='../../../src'
ln -sf ${dir_preprocessing_tools}/*.exe .

# ----------------------------------------
#   Link input files
# ----------------------------------------
export dir_input_ecmwf='../0_input_data'
ln -sf ${dir_input_ecmwf}/concatenate_all_ecmwf_files.nc .

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ------ 
cp namelist_abl_tools_tmpl namelist_abl_tools 
./get_atm_LSfrc.exe   namelist_abl_tools
./vinterp_atm_frc.exe namelist_abl_tools
./drown_atm_frc.exe   namelist_abl_tools
# ------
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


