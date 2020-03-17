#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 16:12:00 2019

@author: deborahkhider

Cut GLDAS data for east africa
"""

import xarray as xr
import numpy as np
import os
import glob
import sys
import ast

bounding_box = ast.literal_eval(sys.argv[1]) #[10,60,-15,30]
path= sys.argv[2] #"./GLDASv2.1"
region = sys.argv[3] #Eastern Equatorial Africa
region = region.replace('_',' ')
output_name = sys.argv[4] #'GLDAS2.1_TP_2000_2018.nc'

folders = os.listdir(path)

years=[]
for folder in folders:
    try:
        year = int(folder)
        years.append(year)
    except:
        continue

years = np.sort(np.array(years))

for i in years:
    nc_files = (glob.glob(path+'/'+str(i)+'/*.nc4'))
    file_names=[]
    for file in nc_files:
        file_names.append(file)
    file_names.sort()
    #open
    data = xr.open_mfdataset(file_names)
    # get the values for the bounding box
    p_ = data.sel(lat=slice(bounding_box[2], bounding_box[3]),\
                  lon=slice(bounding_box[0],bounding_box[1]))
    try: valP =np.concatenate((valP,p_.Rainf_f_tavg.values),axis=0)
    except NameError: valP = p_.Rainf_f_tavg.values
    try: valT =np.concatenate((valT,p_.Tair_f_inst.values),axis=0)
    except NameError: valT = p_.Tair_f_inst.values
    try: time = np.concatenate((time,p_['time'].values),axis=0)
    except NameError: time = p_['time'].values

# Write the thresholds out as netcdf
lat = p_['lat'].values
lon = p_['lon'].values

da_P = xr.DataArray(valP,coords=[time,lat,lon],dims=['time','lat','lon'])
da_P.attrs = p_.Rainf_f_tavg.attrs

da_T = xr.DataArray(valT,coords=[time,lat,lon],dims=['time','lat','lon'])
da_T.attrs = p_.Tair_f_inst.attrs

ds = da_P.to_dataset(name='Rainf_f_tavg')
ds['Tair_f_inst'] = da_T
ds.attrs = p_.attrs
text = ds.attrs['comment']
comment = text + ', files have been cut to center around '+region+'. Only Precipitation and Temperature.'
ds.attrs['comment']=comment
ds.to_netcdf(output_name)