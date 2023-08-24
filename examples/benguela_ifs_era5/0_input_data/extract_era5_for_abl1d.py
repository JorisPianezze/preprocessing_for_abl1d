#!/bin/python
# --------------------------------------------------------
#            Auteur  (date de creation) :
#        J. Pianezze (   20.08.2022   )
# --------------------------------------------------------
# https://cds.climate.copernicus.eu/api-how-to
# conda install cdsapi

import os
import cdsapi
import datetime
import numpy as np

c = cdsapi.Client()

# #########################################################
# ###           to be defined by user                   ###
# #########################################################
# - first_date      = first date to extract (has to be at 00 UTC)
# - last_date       = last  date to extract (has to be at 18 UTC)
# - period_in_hr    = period_in_hr between two forcing files
# - area_to_extract = 'North/West/South/East'
#
first_date          = datetime.datetime(2014,11, 6, 0, 0, 0)
last_date           = datetime.datetime(2014,11, 9,18, 0, 0)
period_in_hr        = 6
area_to_extract     = '-20.0/5.0/-40.0/25.0'
#
# #########################################################

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Compute date and time variables
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
date_an   = str(first_date.year)+'-'+str(first_date.month).zfill(2)+'-'+str(first_date.day).zfill(2)+'/to/'+\
            str( last_date.year)+'-'+str( last_date.month).zfill(2)+'-'+str( last_date.day).zfill(2)
time_an   = np.arange(0,24,period_in_hr)
time_an   = [format(x, '02d') for x in time_an]
file_an   = 'from_'+str(first_date.year)+'-'+str(first_date.month).zfill(2)+'-'+str(first_date.day).zfill(2)+'_at_00_UTC_to_'   +\
                    str( last_date.year)+'-'+str( last_date.month).zfill(2)+'-'+str( last_date.day).zfill(2)+'_at_18_UTC_every_'+\
                    str(   period_in_hr).zfill(2)+'hr'

first_date = first_date - datetime.timedelta(hours=24)
date_fc    = str(first_date.year)+'-'+str(first_date.month).zfill(2)+'-'+str(first_date.day).zfill(2)+'/to/'+\
             str( last_date.year)+'-'+str( last_date.month).zfill(2)+'-'+str( last_date.day).zfill(2)
step_fc    = np.arange(1,13,1)
step_fc    = [format(x, '02d') for x in step_fc]
last_date  = last_date  + datetime.timedelta(hours=24)
file_fc    = 'from_'+str(first_date.year)+'-'+str(first_date.month).zfill(2)+'-'+str(first_date.day).zfill(2)+'_at_07_UTC_to_'   +\
                     str( last_date.year)+'-'+str( last_date.month).zfill(2)+'-'+str( last_date.day).zfill(2)+'_at_06_UTC_every_'+\
                     str(              1).zfill(2)+'hr'

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract Model Level fields : u, v, t et q
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c.retrieve('reanalysis-era5-complete', {
      'date'     : date_an,
      'levelist' : '98/to/137',
      'levtype'  : 'ml',
      'param'    : 'u/v/t/q',
      'stream'   : 'oper',
      'time'     : time_an,
      'type'     : 'an',
      'area'     : area_to_extract,
      'grid'     : '0.28125/0.28125',
  }, 'model_levels_uvtq_'+file_an+'.grib')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract Model Level fields : lnsp
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c.retrieve('reanalysis-era5-complete', {
      'date'     : date_an,
      'levelist' : '1',
      'levtype'  : 'ml',
      'param'    : 'lnsp',
      'stream'   : 'oper',
      'time'     : time_an,
      'type'     : 'an',
      'area'     : area_to_extract,
      'grid'     : '0.28125/0.28125',
  }, 'model_levels_lnsp_'+file_an+'.grib')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract SurFaCe fields : z, lsm, msl
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c.retrieve('reanalysis-era5-complete', {
      'date'     : date_an,
      'levtype'  : 'sfc',
      'param'    : 'z/lsm/msl',
      'stream'   : 'oper',
      'time'     : time_an,
      'type'     : 'an',
      'area'     : area_to_extract,
      'grid'     : '0.28125/0.28125',
  },  'surface_levels_'+file_an+'.grib')
  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract SurFaCe fields : sw, lw and precip
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c.retrieve('reanalysis-era5-complete', {
      'date'     : date_fc,
      'levtype'  : 'sfc',
      'param'    : '228/176/177/175/169',
      'step'     : step_fc,
      'stream'   : 'oper',
      'time'     : ['06','18'],
      'type'     : 'fc',
      'area'     : area_to_extract,
      'grid'     : '0.28125/0.28125',
      }, 'flux_'+file_fc+'.grib')
