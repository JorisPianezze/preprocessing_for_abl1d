#!/bin/bash
# --------------------------------------------------------
#
#                 Author  (    date    ) :
#             J. Pianezze ( 21.08.2023 )
#
#                    ~~~~~~~~~~~~~~~
#       Script used to convert and concatenate ERA5 data
#             for ABL1d (preprocessing tools)
#                    ~~~~~~~~~~~~~~~
#
# --------------------------------------------------------

export list_of_files_to_convert=$(ls era5.*)

for fic in ${list_of_files_to_convert}
  do

  echo ''
  echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  echo '  Convert file:' ${fic}
  echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

  grib_set -s edition=2 ${fic} ${fic}.grb2
  cdo -f nc copy ${fic}.grb2 ${fic}.nc
  rm ${fic}.grb2

  done

export list_of_files_to_concatenate=$(ls era5.*.nc)

echo ''
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo '  List of files to concatenate:' ${list_of_files_to_concatenate}
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
time ncrcat -h ${list_of_files_to_concatenate} concatenate_all_ecmwf_files.nc
