#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 09:39:17 2019

@author: deborahkhider

Drought visualization - Make GIF for now
"""
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib.cm as cm
import imageio
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

#%% SPI and SPEI
#open the dataset
dataset_name = '/Users/deborahkhider/Documents/MINT/Drought/Oct7th indices/SPEI_MINT_2000-01-01_2017-12-01_20ff9f7a-7484-4525-a3a9-f594ea2ce8ab.nc'
dataset = xr.open_dataset(dataset_name)

#get the variable name
varname = list(dataset.data_vars.keys())[0]

#provenance
text = dataset.attrs['summary'].split()
model = text[text.index('from')+1][:-1]

#calib period
calib=  " ".join(text[text.index('calibration'):text.index('distribution')-3])

#get the data
val = dataset[varname].values
if 'latitude' in dataset.keys():
    lat = dataset['latitude'].values
elif 'lat' in dataset.keys():
    lat = dataset['lat'].values

if 'longitude' in dataset.keys():
    lon = dataset['longitude'].values
elif 'lon' in dataset.keys():
    lon = dataset['lon'].values

#Plot
proj = ccrs.PlateCarree(central_longitude = np.mean(lon))
idx = dataset['time'].values.size
count = list(np.arange(0,idx,1))
vmin = np.floor(np.nanmin(val))
vmax = np.ceil(np.nanmax(val))
filenames =[]

#Make a directory if it doesn't exit
if os.path.isdir('./figures') is False:
    os.makedirs('./figures')

for i in count:
    v = val[i,:,:]
    date = pd.to_datetime(dataset.time.values[i]).strftime("%B %Y")  
    fig,ax = plt.subplots(figsize=[15,10])
    ax = plt.axes(projection=proj)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.COASTLINE)
    if np.isnan(v).all() == False:
        img = plt.contourf(lon, lat, v, 60, 
                    transform = ccrs.PlateCarree(),
                    cmap=cm.BrBG,
                    vmin = vmin,
                    vmax = vmax) # need to return to img to make colorbar work
        m = plt.cm.ScalarMappable(cmap=cm.BrBG) # Following three lines necessary to lock ylim on colorbar
        m.set_array(v)
        m.set_clim(vmin, vmax) 
        tick_range = np.arange(-4,4.5,0.5)
        cbar = plt.colorbar(m, orientation='horizontal', ticks = tick_range, pad=0.15)
        long_name = dataset[varname].attrs['long_name']
        units = dataset[varname].attrs['units']
        if units == 'unitless':
            cbar.ax.set_title(long_name)
        else:
            cbar.ax.set_title(long_name+' ('+units+')')
    #cbar.ax.xaxis.set_ticks_position('top')
    #cb.ax.xaxis.set_label_position('top')
    ax.set_extent([23,48,3,15])
    #ax.set_xticks(np.linspace(23,48,4), crs=ccrs.PlateCarree())
    #ax.set_yticks(np.linspace(3,15,4), crs=ccrs.PlateCarree())
    #lon_formatter = LongitudeFormatter(zero_direction_label=True)
    #lat_formatter = LatitudeFormatter()
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlines = False
    gl.ylines = False
    gl.xlocator = mticker.FixedLocator(np.arange(23,49,2))
    gl.ylocator = mticker.FixedLocator(np.arange(3,14,2))
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 12, 'color': 'gray'}
    gl.ylabel_style = {'size': 12, 'color': 'gray'}
    #Add a title with time
    title = date+', inferred from '+model+ ', with a '+calib
    plt.title(title, fontsize=18, loc='left', pad=1)
    #save a jepg
    filename = './figures/'+varname+'_t'+str(i)+'.jpeg'
    filenames.append(filename)
    plt.savefig(filename)
    plt.close(fig)

#create a gif
images = []
for filename in filenames:
    images.append(imageio.imread(filename))

imageio.mimsave(model+'_'+calib+'_'+varname+'_movie.gif', images, duration=1)
    
#%%Pet  
#open the dataset
dataset_name = '/Users/deborahkhider/Documents/MINT/Drought/Oct7th indices/PET_MINT_2000-01-01_2017-12-01_976e5d5f-99f0-4804-90c5-84cc12b7108b.nc'
dataset = xr.open_dataset(dataset_name)

#get the variable name
varname = list(dataset.data_vars.keys())[0]

#provenance
text = dataset.attrs['summary'].split()
model = text[text.index('from')+1][:-1]

#get the data
val = dataset[varname].values
if 'latitude' in dataset.keys():
    lat = dataset['latitude'].values
elif 'lat' in dataset.keys():
    lat = dataset['lat'].values

if 'longitude' in dataset.keys():
    lon = dataset['longitude'].values
elif 'lon' in dataset.keys():
    lon = dataset['lon'].values

#Plot
proj = ccrs.PlateCarree(central_longitude = np.mean(lon))
idx = dataset['time'].values.size
count = list(np.arange(0,idx,1))
vmin = np.floor(np.nanmin(val))
vmax = np.ceil(np.nanmax(val))
filenames =[]

#Make a directory if it doesn't exit
if os.path.isdir('./figures') is False:
    os.makedirs('./figures')

for i in count:
    v = val[i,:,:]
    date = pd.to_datetime(dataset.time.values[i]).strftime("%B %Y")  
    fig,ax = plt.subplots(figsize=[15,10])
    ax = plt.axes(projection=proj)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.COASTLINE)
    if np.isnan(v).all() == False:
        img = plt.contourf(lon, lat, v, 60, 
                    transform = ccrs.PlateCarree(),
                    cmap=cm.cividis,
                    vmin = 0,
                    vmax = 500) # need to return to img to make colorbar work
        m = plt.cm.ScalarMappable(cmap=cm.cividis) # Following three lines necessary to lock ylim on colorbar
        m.set_array(v)
        m.set_clim(0, 500) 
        tick_range = np.arange(0,510,50)
        cbar = plt.colorbar(m, orientation='horizontal', ticks = tick_range, pad=0.15)
        long_name = dataset[varname].attrs['long_name']
        units = dataset[varname].attrs['units']
        if units == 'unitless':
            cbar.ax.set_title(long_name)
        else:
            cbar.ax.set_title(long_name+' ('+units+')')
    #cbar.ax.xaxis.set_ticks_position('top')
    #cb.ax.xaxis.set_label_position('top')
    ax.set_extent([23,48,3,15])
    #ax.set_xticks(np.linspace(23,48,4), crs=ccrs.PlateCarree())
    #ax.set_yticks(np.linspace(3,15,4), crs=ccrs.PlateCarree())
    #lon_formatter = LongitudeFormatter(zero_direction_label=True)
    #lat_formatter = LatitudeFormatter()
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                  linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlines = False
    gl.ylines = False
    gl.xlocator = mticker.FixedLocator(np.arange(23,49,2))
    gl.ylocator = mticker.FixedLocator(np.arange(3,14,2))
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 12, 'color': 'gray'}
    gl.ylabel_style = {'size': 12, 'color': 'gray'}
    #Add a title with time
    title = date+', inferred from '+model
    plt.title(title, fontsize=18, loc='left', pad=1)
    #save a jepg
    filename = './figures/'+varname+'_t'+str(i)+'.jpeg'
    filenames.append(filename)
    plt.savefig(filename)
    plt.close(fig)

#create a gif
images = []
for filename in filenames:
    images.append(imageio.imread(filename))

imageio.mimsave(model+'_'+varname+'_movie.gif', images, duration=1)
    
