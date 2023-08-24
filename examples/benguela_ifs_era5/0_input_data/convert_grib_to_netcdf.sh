#!/bin/bash
# --------------------------------------------------------
#            Auteur  (date de creation) :
#        J. Pianezze (   30.08.2022   )
#
# --------------------------------------------------------

list_files=`ls *.grib`
echo $list_files

for fic in ${list_files}
  do
  echo ''
  echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  echo '  Convert from grib to netcdf : ' ${fic}
  echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

  #grib_set -s edition=2 ${fic} ${fic}.grb2
  #cdo -f nc copy ${fic}.grb2 ${fic}.nc
  #rm ${fic}.grb2

  cdo -f nc copy ${fic} ${fic}.nc

  done

echo ''
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
echo '  Append all variables in one file                     '
echo '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

export file_uvtq=`ls *uvtq*nc`
export file_lnsp=`ls *lnsp*nc`
export file_surf=`ls surface*nc`

cp $file_uvtq era5_for_preprocessing_tools.nc

ncwa -C -O -a lev    $file_lnsp $file_lnsp
ncks -x -O -v lev    $file_lnsp $file_lnsp
ncks -A -v lnsp      $file_lnsp era5_for_preprocessing_tools.nc
ncks -A -v z,lsm,msl $file_surf era5_for_preprocessing_tools.nc

rm $file_lnsp $file_surf $file_uvtq
