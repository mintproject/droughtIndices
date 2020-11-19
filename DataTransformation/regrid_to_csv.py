#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 17:24:20 2020

@author: deborahkhider

Regrid data for Ethiopia only and transform a netcdf file into csv
"""
import xarray as xr
import numpy as np
import pandas as pd
from spharm import Spharmt, regrid
import sys
import uuid
import ast

def openNetcdf(filename,bounding_box):

    ds=xr.open_dataset(filename)
    varname=list(ds.data_vars.keys())[0]
    if 'latitude' in ds.coords:
        p_ = ds.sel(latitude=slice(bounding_box[2], bounding_box[3]),\
                  longitude=slice(bounding_box[0],bounding_box[1]))
        da=p_[varname]
        lat = da['latitude'].values
        lon = da['longitude'].values
    elif 'lat' in ds.coords:
        p_ = ds.sel(lat=slice(bounding_box[2], bounding_box[3]),\
                  lon=slice(bounding_box[0],bounding_box[1]))
        da=p_[varname]
        lat = da['lat'].values
        lon = da['lon'].values
    elif 'Y' in ds.coords:
        p_ = ds.sel(Y=slice(bounding_box[2], bounding_box[3]),\
                  X=slice(bounding_box[0],bounding_box[1]))
        da=p_[varname]
        lat = da['Y'].values
        lon = da['X'].values

    return ds, da, lat, lon


def regrid_field(field, lat, lon, lat_new, lon_new):
    nlat_old, nlon_old = np.size(lat), np.size(lon)
    nlat_new, nlon_new = np.size(lat_new), np.size(lon_new)
    spec_old = Spharmt(nlon_old, nlat_old, gridtype='regular', legfunc='computed')
    spec_new = Spharmt(nlon_new, nlat_new, gridtype='regular', legfunc='computed')
    #remove nans
    field[np.isnan(field)] = 0
    field_new = []
    for field_old in field:
        regridded_field =  regrid(spec_old, spec_new, field_old, ntrunc=None, smooth=None)
        field_new.append(regridded_field)

    field_new = np.array(field_new)
    return field_new

def toNetcdf(ds,field_new,lat_new,lon_new):
    time=ds['time'].values
    varname=list(ds.data_vars.keys())[0]
    da_new=xr.DataArray(field_new,coords=[time,lat_new,lon_new],dims=['time','latitude','longitude'])
    ds2=da_new.to_dataset(name=varname)
    ds2.attrs=ds.attrs
    ds2[varname].attrs=ds[varname].attrs
    ds2.attrs['id'] = str(uuid.uuid4())
    ds2.to_netcdf(path = 'regrid_results.nc')

def toCsv(ds,field_new,lat_new,lon_new):
    time=ds['time'].values
    varname=list(ds.data_vars.keys())[0]

    latitude=[]
    longitude=[]
    time2=[]
    value=[]
    for idx,t in enumerate(time):
        for idx2, l in enumerate(lat_new):
            for idx3, lo in enumerate(lon_new):
                latitude.append(l)
                longitude.append(lo)
                time2.append(t)
                value.append(field_new[idx,idx2,idx3])
    #create a pandas dataframe
    d = {'time':time2,
         'latitude':latitude,
         'longitude':longitude,
         varname:value}
    df=pd.DataFrame(data=d)
    df.to_csv('regrid_results.csv')

if __name__ == "__main__":
    #bounding box set for ethiopia. TODO: make into an argument
    bounding_box= [32.5,45.8,2.5,15]
    ds, da, lat, lon = openNetcdf(sys.argv[1],bounding_box)
    # Enter new lats and lons
    
    #regrid
    if ast.literal_eval(sys.argv[2]) == True:
        lat_new=np.linspace(np.min(lat),np.max(lat),12)
        lon_new=np.linspace(np.min(lon),np.max(lon),12)
        field_new = regrid_field(da.values,
                              lat,
                              lon,
                              lat_new,
                              lon_new)
        print('ok')
    else:
        field_new=da.values
        lat_new = lat
        lon_new = lon
    toNetcdf(ds,field_new,lat_new,lon_new)
    toCsv(ds,field_new,lat_new,lon_new)
