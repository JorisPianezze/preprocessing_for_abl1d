#!/bin/python
# --------------------------------------------------------
#
#                 Author  (    date    ) :
#             J. Pianezze ( 14.08.2023 )
#
#                    ~~~~~~~~~~~~~~~
#       Script used to extract ERA5 instantaneous fields
#  for Meso-NH (PREP_REAL_CASE) or ABL1d (preprocessing tools)
#                    (1 time / file)
#                    ~~~~~~~~~~~~~~~
#
# --------------------------------------------------------
# https://cds.climate.copernicus.eu/api-how-to
# conda install cdsapi

import os, sys
import glob
import cdsapi
import datetime

cds = cdsapi.Client()

# #########################################################
# ###           To be defined by user                   ###
# #########################################################

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - -     First and last date to be extracted           - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
first_date_to_be_extracted = datetime.datetime(2005, 1, 1, 0, 0, 0)
last_date_to_be_extracted  = datetime.datetime(2005, 1, 1, 6, 0, 0)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - -          Type of data to be extracted             - -
# - -         analyses (an) or forecast (fc)            - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
type_data_to_be_extracted = 'an'

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - -     period_in_hr between two forcing files        - -
# - -          must be a multiple of 6                  - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
period_between_last_and_first_dates_in_hr = 6

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - -     Area to be extracted : 'North/West/South/East'- -
# - -     Benguela    : '-20.0/5.0/-40.0/25.0'          - -
# - -     Gulf Stream : '50.0/-90.0/20.0/-30.0'         - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
area_to_be_extracted = '50.0/-90.0/20.0/-30.0'

# #########################################################

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Define function to iterate over first and last dates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def range_for_date(start_date, end_date, period_in_hr):
  for n in range(int((end_date - start_date).total_seconds()/(3600.0*period_in_hr))+1):
    yield start_date + datetime.timedelta(seconds=n*3600.0*period_in_hr)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#   Loop over dates
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for date in range_for_date(first_date_to_be_extracted, last_date_to_be_extracted, period_between_last_and_first_dates_in_hr):

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Compute date and time variables
  #     date_an format is yyyy-mm-dd
  #     time_an format is hh
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  date_to_be_extracted   = str(date.year)+'-'+str(date.month).zfill(2)+'-'+str(date.day).zfill(2)
  time_to_be_extracted   = str(date.hour).zfill(2)
  name_of_extracted_file = str(date.year)+str(date.month).zfill(2)+str(date.day).zfill(2)+'.'+str(date.hour).zfill(2)

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Retrieve Model Level fields : u, v, t, q
  #     for Meso-NH and ABL1d
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cds.retrieve('reanalysis-era5-complete', {
        'date'     : date_to_be_extracted,
        'levelist' : '1/to/137',
        'levtype'  : 'ml',
        'param'    : 'u/v/t/q',
        'stream'   : 'oper',
        'time'     : time_to_be_extracted,
        'type'     : type_data_to_be_extracted,
        'area'     : area_to_be_extracted,
        'grid'     : '0.28125/0.28125',
    }, 'model_levels_uvtq_'+name_of_extracted_file+'.grib')

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Retrieve Model Level fields : lnsp
  #     for Meso-NH and ABL1d
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cds.retrieve('reanalysis-era5-complete', {
        'date'     : date_to_be_extracted,
        'levelist' : '1',
        'levtype'  : 'ml',
        'param'    : 'lnsp',
        'stream'   : 'oper',
        'time'     : time_to_be_extracted,
        'type'     : type_data_to_be_extracted,
        'area'     : area_to_be_extracted,
        'grid'     : '0.28125/0.28125',
    }, 'model_levels_lnsp_'+name_of_extracted_file+'.grib')

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Retrieve SurFaCe fields : z, lsm
  #     for Meso-NH and ABL1d
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cds.retrieve('reanalysis-era5-complete', {
        'date'     : date_to_be_extracted,
        'levtype'  : 'sfc',
        'param'    : 'z/lsm',
        'stream'   : 'oper',
        'time'     : time_to_be_extracted,
        'type'     : type_data_to_be_extracted,
        'area'     : area_to_be_extracted,
        'grid'     : '0.28125/0.28125',
    },  'surface_levels_'+name_of_extracted_file+'.grib')

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Retrieve SurFaCe fields : swlv1, ...
  #     for Meso-NH only
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cds.retrieve('reanalysis-era5-complete', {
      'date'     : date_to_be_extracted,
      'levtype'  : 'sfc',
      'param'    : '139/141/170/183/236/39/40/41/42',
      'stream'   : 'oper',
      'time'     : time_to_be_extracted,
      'type'     : type_data_to_be_extracted,
      'area'     : area_to_be_extracted,
      'grid'     : '0.28125/0.28125',
  }, 'surface_'+name_of_extracted_file+'.grib')

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Retrieve SurFaCe fields : msl
  #     for ABL1d only
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  cds.retrieve('reanalysis-era5-complete', {
        'date'     : date_to_be_extracted,
        'levtype'  : 'sfc',
        'param'    : 'msl',
        'stream'   : 'oper',
        'time'     : time_to_be_extracted,
        'type'     : type_data_to_be_extracted,
        'area'     : area_to_be_extracted,
        'grid'     : '0.28125/0.28125',
    },  'surface_msl_'+name_of_extracted_file+'.grib')

  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #   Concatenate & remove grib files
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  os.system('grib_copy surface_levels_'+name_of_extracted_file+'.grib    '+\
                      'surface_'+name_of_extracted_file+'.grib           '+\
                      'surface_msl_'+name_of_extracted_file+'.grib       '+\
                      'model_levels_uvtq_'+name_of_extracted_file+'.grib '+\
                      'model_levels_lnsp_'+name_of_extracted_file+'.grib '+\
                      'era5.'+name_of_extracted_file)

  for file in glob.glob('*.grib'):
    os.remove(file)
