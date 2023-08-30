#!/usr/bin/python
# -*- coding: utf-8 -*-
# ----------------------------------------------------
#       Auteur (date de creation)
#  J. PIANEZZE (      31.08.2022)
# ----------------------------------------------------
#
print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
print('&&&                                           ')
print('&&&   Execution de :                          ')
print('&&&  ',__file__                                )
print('&&&                                           ')
print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
#
######################################################
###                                           
import os, sys
import netCDF4
import numpy as np
import xarray
import datetime
###                                           
######################################################

cfg_debug           = True

cfg_dir_input_files = '../1_run_preprocessing_tool_N_abl_50/'

# ~~~ Has to be the same as extract_era5_for_abl1d.py ~~~
first_date          = datetime.datetime(2014,11, 6, 0, 0, 0)
last_date           = datetime.datetime(2014,11, 9,18, 0, 0)
period_in_hr        = 6

if cfg_debug : print('++++++++++++++++++++++++++++++++++++++++++++++')
if cfg_debug : print('+++                                           ')
if cfg_debug : print('+++   0. Lecture des fichiers et des variables')
if cfg_debug : print('+++                                           ')
if cfg_debug : print('++++++++++++++++++++++++++++++++++++++++++++++')

# ---------------------------------------------------------------------
#  Get datasets
# ---------------------------------------------------------------------
file_croco_grd = netCDF4.Dataset('croco_grd.nc')
file_abl_grd   = netCDF4.Dataset(cfg_dir_input_files+'vertical_grid_N_abl_50_H_2000m_DZ_10m.nc')
file_abl_ifs   = netCDF4.Dataset(cfg_dir_input_files+'out_drown_era5_for_preprocessing_tools.nc')

# ---------------------------------------------------------------------
#  Read ABL vertical grid
# ---------------------------------------------------------------------
N_abl   = file_abl_grd.dimensions['jpka'] ; N_abl = N_abl.size - 1

zr_abl  = file_abl_grd.variables['ght'][:N_abl]
zw_abl  = file_abl_grd.variables['ghw'][:N_abl]
Hzr_abl = file_abl_grd.variables['e3t'][:N_abl]
Hzw_abl = file_abl_grd.variables['e3w'][:N_abl]

# ---------------------------------------------------------------------
#  Read IFS forcing time
# ---------------------------------------------------------------------
time_abl = file_abl_ifs.variables['time'] ; ntime_abl = time_abl.size

# ---------------------------------------------------------------------
#  Read CROCO dimensions
# ---------------------------------------------------------------------
xi_rho  = file_croco_grd.dimensions['xi_rho']  ; xi_rho  = xi_rho.size
eta_rho = file_croco_grd.dimensions['eta_rho'] ; eta_rho = eta_rho.size
xi_u    = file_croco_grd.dimensions['xi_u']    ; xi_u    = xi_u.size
eta_u   = file_croco_grd.dimensions['eta_u']   ; eta_u   = eta_u.size
xi_v    = file_croco_grd.dimensions['xi_v']    ; xi_v    = xi_v.size
eta_v   = file_croco_grd.dimensions['eta_v']   ; eta_v   = eta_v.size

if cfg_debug : print('++++++++++++++++++++++++++++++++++++++++++++++')
if cfg_debug : print('+++                                           ')
if cfg_debug : print('+++   1. Interpolation grille IFS a CROCO     ')
if cfg_debug : print('+++                                           ')
if cfg_debug : print('++++++++++++++++++++++++++++++++++++++++++++++')

# --------------------------------------------
#   Read CROCO grd
# ---------------------------------------------
ds_croco_grd = xarray.open_dataset('croco_grd.nc')
lon_rho      = ds_croco_grd.lon_rho # - 360.0
lat_rho      = ds_croco_grd.lat_rho
lon_u        = ds_croco_grd.lon_u # - 360.0
lat_u        = ds_croco_grd.lat_u
lon_v        = ds_croco_grd.lon_v # - 360.0
lat_v        = ds_croco_grd.lat_v

# --------------------------------------------
#   Read IFS grd
# ---------------------------------------------
ds_ifs       = xarray.open_dataset(cfg_dir_input_files+'../0_input_data/flux_from_2014-11-05_at_07_UTC_to_2014-11-10_at_06_UTC_every_01hr.grib.nc')
ds_abl_ifs   = xarray.open_dataset(cfg_dir_input_files+'out_drown_era5_for_preprocessing_tools.nc')

# --------------------------------------------
#   Read IFS variables
# ---------------------------------------------
tair_ifs     = ds_abl_ifs.tair
rhum_ifs     = ds_abl_ifs.humi
uwnd_ifs     = ds_abl_ifs.uwnd
vwnd_ifs     = ds_abl_ifs.vwnd
prate_ifs    = ds_ifs.tp
radlw_ifs    = ds_ifs.str
radlw_in_ifs = ds_ifs.strd
radsw_ifs    = ds_ifs.ssr

# --------------------------------------------
#   Get 10m variables
# ---------------------------------------------
tair_ifs     = tair_ifs.sel(jpka=1)
rhum_ifs     = rhum_ifs.sel(jpka=1)
uwnd_ifs     = uwnd_ifs.sel(jpka=1)
vwnd_ifs     = vwnd_ifs.sel(jpka=1)

# --------------------------------------------
#   Special treatment for fluxes and precipitation
# ---------------------------------------------
# --> 6h-average, period_in_hr between two forcing files
prate_ifs    = prate_ifs.resample   (time=str(period_in_hr)+'H',label='right').mean()
radlw_ifs    = radlw_ifs.resample   (time=str(period_in_hr)+'H',label='right').mean()
radlw_in_ifs = radlw_in_ifs.resample(time=str(period_in_hr)+'H',label='right').mean()
radsw_ifs    = radsw_ifs.resample   (time=str(period_in_hr)+'H',label='right').mean()

# --> date selection : date_first -> date_last
prate_ifs    = prate_ifs.sel   (time=slice(first_date, last_date))
radlw_ifs    = radlw_ifs.sel   (time=slice(first_date, last_date))
radlw_in_ifs = radlw_in_ifs.sel(time=slice(first_date, last_date))
radsw_ifs    = radsw_ifs.sel   (time=slice(first_date, last_date))

# --------------------------------------------
#   Interpolate at rho-points
# ---------------------------------------------
tair_abl     = tair_ifs.interp    (lat=lat_rho, lon=lon_rho, method='linear')
rhum_abl     = rhum_ifs.interp    (lat=lat_rho, lon=lon_rho, method='linear')
uwnd_abl     = uwnd_ifs.interp    (lat=lat_u  , lon=lon_u  , method='linear')
vwnd_abl     = vwnd_ifs.interp    (lat=lat_v  , lon=lon_v  , method='linear')
prate_abl    = prate_ifs.interp   (lat=lat_rho, lon=lon_rho, method='linear')
radlw_abl    = radlw_ifs.interp   (lat=lat_rho, lon=lon_rho, method='linear')
radlw_in_abl = radlw_in_ifs.interp(lat=lat_rho, lon=lon_rho, method='linear')
radsw_abl    = radsw_ifs.interp   (lat=lat_rho, lon=lon_rho, method='linear')

# ----------------------
#   remove negative prate
# -----------------------
prate_abl = np.where(prate_abl<0.0, 0.0, prate_abl)

if cfg_debug : print('++++++++++++++++++++++++++++++++++++++++++++++')
if cfg_debug : print('+++                                           ')
if cfg_debug : print('+++   2. Ecriture du fichier netcdf           ')
if cfg_debug : print('+++                                           ')
if cfg_debug : print('++++++++++++++++++++++++++++++++++++++++++++++')

file_croco_abl_ifs=netCDF4.Dataset(os.getcwd()+'/croco_blk_from_ifs_era5_interp_linear.nc','w')
file_croco_abl_ifs.Description='Atmospheric forcing for CROCO'

# ----------------------------------
# Create the dimensions of the files
# ----------------------------------
file_croco_abl_ifs.createDimension ('xi_rho'   , xi_rho )
file_croco_abl_ifs.createDimension ('eta_rho'  , eta_rho)
file_croco_abl_ifs.createDimension ('xi_u'     , xi_u   )
file_croco_abl_ifs.createDimension ('eta_u'    , eta_u  )
file_croco_abl_ifs.createDimension ('xi_v'     , xi_v   )
file_croco_abl_ifs.createDimension ('eta_v'    , eta_v  )
file_croco_abl_ifs.createDimension ('bulk_time', None   )

# ----------------------------------
# Create the variables of the files
# ----------------------------------
varout = file_croco_abl_ifs.createVariable('bulk_time','d',('bulk_time'                     )) ; varout.long_name = 'bulk formulation execution time' ; varout.units = 'days'    ; varout.cycle_length = 0.
varout = file_croco_abl_ifs.createVariable('tair'     ,'d',('bulk_time', 'eta_rho', 'xi_rho')) ; varout.long_name = 'Air temperature'                 ; varout.units = 'Celsius'
varout = file_croco_abl_ifs.createVariable('rhum'     ,'d',('bulk_time', 'eta_rho', 'xi_rho')) ; varout.long_name = 'Air specific humidity'           ; varout.units = 'g/kg'
varout = file_croco_abl_ifs.createVariable('prate'    ,'d',('bulk_time', 'eta_rho', 'xi_rho')) ; varout.long_name = 'Precipitation rate'              ; varout.units = 'cm/day'
varout = file_croco_abl_ifs.createVariable('radlw'    ,'d',('bulk_time', 'eta_rho', 'xi_rho')) ; varout.long_name = 'Net longwave radiation'          ; varout.units = 'W/m2'
varout = file_croco_abl_ifs.createVariable('radlw_in' ,'d',('bulk_time', 'eta_rho', 'xi_rho')) ; varout.long_name = 'Downward longwave radiation'     ; varout.units = 'W/m2'
varout = file_croco_abl_ifs.createVariable('radsw'    ,'d',('bulk_time', 'eta_rho', 'xi_rho')) ; varout.long_name = 'Net shortwave radiation'         ; varout.units = 'W/m2'
varout = file_croco_abl_ifs.createVariable('uwnd'     ,'d',('bulk_time', 'eta_u'  , 'xi_u'  )) ; varout.long_name = 'u-wind speed'                    ; varout.units = 'm/s'
varout = file_croco_abl_ifs.createVariable('vwnd'     ,'d',('bulk_time', 'eta_v'  , 'xi_v'  )) ; varout.long_name = 'v-wind speed'                    ; varout.units = 'm/s'

# ---------------------------------------
# Write out the data arrays into the file
# ---------------------------------------
file_croco_abl_ifs.variables['bulk_time'][:    ] = time_abl     [:    ]/24.0               # from hours to days
file_croco_abl_ifs.variables['tair']     [:,:,:] = tair_abl     [:,:,:]-273.15             # from K to C
file_croco_abl_ifs.variables['rhum']     [:,:,:] = rhum_abl     [:,:,:]*1E3                # from kg/kg to g/kg 
file_croco_abl_ifs.variables['prate']    [:,:,:] = prate_abl    [:,:,:]*1E2/(1.0/24.0)/1E3 # from kg/m2 to cm/day
file_croco_abl_ifs.variables['radlw']    [:,:,:] = radlw_abl    [:,:,:]/(1.0*3600.0)       # from J/m2 to W/m2
file_croco_abl_ifs.variables['radlw_in'] [:,:,:] = radlw_in_abl [:,:,:]/(1.0*3600.0)       # from J/m2 to W/m2
file_croco_abl_ifs.variables['radsw']    [:,:,:] = radsw_abl    [:,:,:]/(1.0*3600.0)       # from J/m2 to W/m2
file_croco_abl_ifs.variables['uwnd']     [:,:,:] = uwnd_abl     [:,:,:]
file_croco_abl_ifs.variables['vwnd']     [:,:,:] = vwnd_abl     [:,:,:]

# ---------------------------------------
# close the file
# ---------------------------------------
file_croco_abl_ifs.close()

print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
print('&&&                                           ')
print('&&&   Fin de :                                ')
print('&&&  ',__file__                                )
print('&&&                                           ')
print('&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&')
