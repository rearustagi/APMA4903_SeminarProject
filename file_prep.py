# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 12:01:29 2021

@author: rearu
"""
import numpy as np
import xarray as xr

months = ['JUN', 'JUL', 'AUG']
years1980s = np.arange(1980, 1990)
years2010s = np.arange(2010, 2020)

ds1=xr.open_dataset(path+'JUN1980.aijh12iWISO_20th_MERRA2_ANL.nc')
for y in years1980s:
    for m in months:
        if [m, y] != ["JUN", 1980]:
            fname = '{}{}.aijh12iWISO_20th_MERRA2_ANL.nc'.format(m, y)
            ds2=xr.open_dataset(path+fname)
            ds1 = xr.concat([ds1, ds2], dim="time")
        
#print(ds1.time.count().values)
ds1.to_netcdf(path=path+'1980s_prec.nc')

ds3=xr.open_dataset(path+'JUN2010.aijh12iWISO_20th_MERRA2_ANL.nc')
for y in years2010s:
    for m in months:
        if [m, y] != ["JUN", 2010]:
            fname = '{}{}.aijh12iWISO_20th_MERRA2_ANL.nc'.format(m, y)
            ds4=xr.open_dataset(path+fname)
            ds3 = xr.concat([ds3, ds4], dim="time")
        
#print(ds1.time.count().values)
ds1.to_netcdf(path=path+'2010s_prec.nc')