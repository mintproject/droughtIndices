
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 17:15:34 2019

@author: deborahkhider
@description: Climate indices (drought) for DARPA World's modeler program
"""
import xarray as xr
import pandas as pd
import numpy as np
from climate_indices import compute,indices
import uuid
from datetime import date
import os
import glob
from calendar import monthrange
import sys
import ast
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import imageio
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.cm as cm


#%% Open files from various agencies
## CHIRPS
def openCHIRPS(dataset_name, bounding_box):
    """ Open CHIRPS dataset and returns the data

    Args:
        dataset_name (str): The name of the CHIRPS dataset
        bounding_box (list): lat/lon to cut to appropriate size

    Returns:
        da_precip (Xarra DataArray): A dataArray of precipitation

    """
    data = xr.open_dataset(dataset_name)
    p_ = data.sel(latitude=slice(bounding_box[2], bounding_box[3]),\
                  longitude=slice(bounding_box[0],bounding_box[1]))
    da_precip = p_.precip
    #da_precip_groupby = da_precip.stack(point=('latitude', 'longitude')).groupby('point')

    return da_precip

def openGLDAS(dataset_name, bounding_box, periodicity, netcdf):
    """Open GLDAS datasets and return precipitation and temperature

    Args:
        dataset_name (str): The name of the GLDAS folder
        bounding_box (list): lat/lon to cut to appropriate size
        periodicity (str): The temporal resolution of the input data. Useful to
            know how the data is organized.
        netcdf (bool): Whether the input data is in netCDF format.

    Returns:
        da_precip (Xarra DataArray): A dataArray of precipitation
        da_precip_groupby (Xarray DataArray): Xarray dataarray containing
            the precipitation data grouped by lat/lon
        da_temp (Xarra DataArray): A dataArray of temperature
        da_temp_groupby (Xarray DataArray): Xarray dataarray containing
            the temperature data grouped by lat/lon
    """

    # Loop over all datasets in the various folders to get the data
    if netcdf == False:
        file_names = []
        if periodicity == 'monthly':
            subdirs = [os.path.join(('./'+dataset_name),o) for o in os.listdir('./'+dataset_name)
                        if os.path.isdir(os.path.join(('./'+dataset_name),o))]
            for path in subdirs:
                nc_files = (glob.glob(path+'/*.nc4'))
                for file in nc_files:
                    file_names.append(file)
        file_names.sort()
        data = xr.open_mfdataset(file_names)
    if netcdf == True:
        data = xr.open_dataset(dataset_name)
    p_ = data.sel(lat=slice(bounding_box[2], bounding_box[3]),\
                  lon=slice(bounding_box[0],bounding_box[1]))
    ### PRECIPITATION
    da_precip = p_.Rainf_f_tavg
    ##unit transformation from kgm-2s-1 to mm
    #days in the month
    time = pd.to_datetime(da_precip.time.values)
    v= time.strftime("%Y,%m")
    days=[]
    for item in v:
        days.append(monthrange(int(item.split(',')[0]),int(item.split(',')[1]))[1])
    prcp = da_precip.values
    prcp = prcp*86400 #convert from kg/m2/s to mm/day
    for i in np.arange(0,len(days)):
        prcp[i,:,:] = prcp[i,:,:]*days[i]
    da_precip.values = prcp
    da_precip.attrs['units'] = 'mm'
    ### TEMPERATURE
    da_temp = p_.Tair_f_inst
    ## Unit conversion
    t = da_temp.values
    t = t-273.15
    da_temp.values = t
    da_temp.attrs['units'] = 'C'

    return da_precip, da_temp

def openFLDAS(dataset_name, bounding_box, periodicity, netcdf):
    """Open FLDAS datasets and return precipitation and temperature

    Args:
        dataset_name (str): The name of the GLDAS folder
        bounding_box (list): lat/lon to cut to appropriate size
        periodicity (str): The temporal resolution of the input data. Useful to
            know how the data is organized.
        netcdf (bool): Whether the input data is in netCDF format.

    Returns:
        da_precip (Xarra DataArray): A dataArray of precipitation
        da_precip_groupby (Xarray DataArray): Xarray dataarray containing
            the precipitation data grouped by lat/lon
        da_temp (Xarra DataArray): A dataArray of temperature
        da_temp_groupby (Xarray DataArray): Xarray dataarray containing
            the temperature data grouped by lat/lon
    """

    # Loop over all datasets in the various folders to get the data
    if netcdf == False:
        file_names = []
        if periodicity == 'monthly':
            subdirs = [os.path.join(('./'+dataset_name),o) for o in os.listdir('./'+dataset_name)
                        if os.path.isdir(os.path.join(('./'+dataset_name),o))]
            for path in subdirs:
                nc_files = (glob.glob(path+'/*.nc'))
                for file in nc_files:
                    file_names.append(file)
        file_names.sort()
        data = xr.open_mfdataset(file_names)
    if netcdf == True:
        data = xr.open_dataset(dataset_name)
    p_ = data.sel(Y=slice(bounding_box[2], bounding_box[3]),\
                  X=slice(bounding_box[0],bounding_box[1]))
    ### PRECIPITATION
    da_precip = p_.Rainf_f_tavg
    ##unit transformation from kgm-2s-1 to mm
    #days in the month
    time = pd.to_datetime(da_precip.time.values)
    v= time.strftime("%Y,%m")
    days=[]
    for item in v:
        days.append(monthrange(int(item.split(',')[0]),int(item.split(',')[1]))[1])
    prcp = da_precip.values
    prcp = prcp*86400 #convert from kg/m2/s to mm/day
    for i in np.arange(0,len(days)):
        prcp[i,:,:] = prcp[i,:,:]*days[i]
    da_precip.values = prcp
    da_precip.attrs['units'] = 'mm'
    ### TEMPERATURE
    da_temp = p_.Tair_f_tavg
    ## Unit conversion
    t = da_temp.values
    t = t-273.15
    da_temp.values = t
    da_temp.attrs['units'] = 'C'

    return da_precip, da_temp
#%% Compute indices from xarray-like dataset

##SPI
def SPI(da_precip, distribution = 'gamma', periodicity = 'monthly',\
        scales = 6, data_start_year = 'beginning', data_end_year = 'end',\
        calibration_start_year = 'beginning',\
        calibration_end_year = 'end'):
    """Calculate SPI from precipitation

    This function uses the SPI calculation from the climate indices package

    Args:
        da_precip (Xarray DataArray): A dataArray of precipitation
        distribution (str): The distribution used to fit the data. Default is
            'gamma'. To use a Pearson Type III distribution, enter 'pearson'
        periodicity (str): Either 'monthly' or 'daily'
        scales (int): The timescales on which the index is computed, either 6 or 12.
            Default is 6.
        data_start_year: Year to start computing  SPI - Default is first year in the data
        data_end_year: Year to stop computing  SPI - Default is first year in the data
        calibration_start_year: Start year for the calibration - Defauls is first year in the data
        calibration_end_year: End year for the calibration - Default is to set a 30-year  period,
            or the full dataset if shorter than 30 years

    Returns:
        ds_spi (Xarray DataArray): SPI index  DataArray
        info (dict): Dictionary containing relevant information about the calib period
    """

    #Perform some checks
    assert scales==6 or scales==12, 'Valid entries for timescales field should be 6 or 12'
    assert distribution == 'gamma' or distribution == 'pearson', "Valid entries for distribution field should be 'gamma' or pearson'"
    assert periodicity == 'monthly' or periodicity == 'daily', "Valid entries for periodicity field should be 'monthly' or 'daily'"

    if distribution == 'gamma':
        dist = indices.Distribution.gamma
    elif distribution == 'pearson':
        dist = indices.Distribution.pearson

    if periodicity == 'monthly':
        period = compute.Periodicity.monthly
    if periodicity == 'daily':
        period = compute.Periodicity.daily

    if data_start_year == 'beginning':
        data_start_year = int(np.min(da_precip.time.dt.year))
    elif data_start_year not in da_precip.time.dt.year.values:
        print('Start year not in dataset, using first available year')
        data_start_year = int(np.min(da_precip.time.dt.year))

    if data_end_year == 'end':
        data_end_year = int(np.max(da_precip.time.dt.year))
    elif data_end_year not in da_precip.time.dt.year.values:
        print('Start year not in dataset, using first available year')
        data_end_year = int(np.max(da_precip.time.dt.year))

    if calibration_start_year == 'beginning':
        calibration_start_year = int(np.min(da_precip.time.dt.year))
    elif calibration_start_year not in da_precip.time.dt.year.values:
        print('Calibration start year not in dataset, using first available year')
        calibration_start_year = int(np.min(da_precip.time.dt.year))

    if calibration_end_year == 'end':
        calibration_end_year = calibration_start_year + 29
        if calibration_end_year not in da_precip.time.dt.year.values:
            calibration_end_year = int(np.max(da_precip.time.dt.year))
    elif calibration_end_year not in da_precip.time.dt.year.values:
        print('Calibration end year not in dataset, using 30yr window or end of dataset')
        calibration_end_year = calibration_start_year + 29
        if calibration_end_year not in da_precip.time.dt.year.values:
            calibration_end_year = int(np.max(da_precip.time.dt.year))

    #Groupby
    if 'lat' in da_precip.coords:
        da_precip_groupby = da_precip.stack(point=('lat', 'lon')).groupby('point')
    elif 'latitude' in da_precip.coords:
        da_precip_groupby = da_precip_groupby = da_precip.stack(point=('latitude', 'longitude')).groupby('point')
    elif 'Y' in da_precip.coords:
        da_precip_groupby = da_precip.stack(point=('Y', 'X')).groupby('point')
    else:
        raise KeyError('latitude not found')

    #Perform calculation
    da_spi = xr.apply_ufunc(indices.spi,
                        da_precip_groupby,
                        scales,
                        dist,
                        data_start_year,
                        calibration_start_year,
                        calibration_end_year,
                        period)
    #Unstack
    da_spi = da_spi.unstack('point')
    #cut
    min_idx = np.where(da_spi.time.dt.year>=data_start_year)[0][0]
    max_idx = np.where(da_spi.time.dt.year<=data_end_year)[0][-1]
    r = np.arange(min_idx,max_idx+1,1)
    #Cut along the right dimensions
    da_spi_cut = np.take(da_spi,r,axis=da_spi.dims.index('time'))
    da_spi_cut['time'] = da_spi_cut.time[r]
    ds_spi=da_spi_cut.to_dataset(name='spi') #convert to xarrayDataset

    info = {'index':'SPI',
            'calibration_start': calibration_start_year,
            'calibration_end':  calibration_end_year,
            'distribution': distribution,
            'periodicity': periodicity,
            'timescales': scales}

    return ds_spi, info

## PET
def PET(da_temp, data_start_year = 'beginning', data_end_year = 'end'):
    """Calculate PET from temperature

    This function uses the PET calculation from the climate indices package.
    Monhtly only.
    Uses Thornthwaite equation.

    Args:
        da_temp (Xarray DataArray): A dataArray of temperature
        data_start_year: Year to start computing  SPI - Default is first year in the data
        data_end_year: Year to stop computing  SPI - Default is first year in the data

    Returns:
        ds_pet (Xarray DataSet): PET index (mm/month)
        da_pet (Xarray DataArray): PET index for full length of data (mm/month). To use with SPEI function
        info (dict): Dictionary containing relevant information
    """

    if data_start_year == 'beginning':
        data_start_year = int(np.min(da_precip.time.dt.year))
    elif data_start_year not in da_precip.time.dt.year.values:
        print('Start year not in dataset, using first available year')
        data_start_year = int(np.min(da_precip.time.dt.year))

    if data_end_year == 'end':
        data_end_year = int(np.max(da_precip.time.dt.year))
    elif data_end_year not in da_precip.time.dt.year.values:
        print('End year not in dataset, using last available year')
        data_end_year = int(np.max(da_precip.time.dt.year))

    #get latitude
    if 'lat' in da_temp.coords:
        lat = np.array(da_temp['lat'].values)
    elif 'latitude' in da_temp.coords:
        lat = np.array(da_temp['latitude'].values)
    elif 'Y' in da_temp.coords:
        lat = np.array(da_temp['Y'].values)
    else:
        raise KeyError('latitude not found')

    #Make sure that we have full years in the data or truncate
    if int(da_precip.time[0].dt.month) != 1:
        print("Full year not available for the beginning of the record, truncating...")
        start_year = int(da_precip.time.dt.year[1])
    else:
        start_year = int(da_precip.time.dt.year[0])

    if int(da_precip.time[-1].dt.month) != 12:
        print("Full year not available for the end of the record, truncating...")
        end_year = int(da_precip.time.dt.year[-2])
    else:
        end_year = int(da_precip.time.dt.year[-1])

    #cut the data if entire years are not available.
    min_idx = np.where(da_temp.time.dt.year>=start_year)[0][0]
    max_idx = np.where(da_temp.time.dt.year<=end_year)[0][-1]
    r = np.arange(min_idx,max_idx+1,1)
    #cut along the right dimensions
    da_temp_cut = np.take(da_temp,r,axis=da_temp.dims.index('time'))
    da_temp_cut['time'] = da_temp_cut.time[r]

    #Groupby
    if 'lat' in da_temp.coords:
        da_temp_groupby = da_temp.stack(point=('lat', 'lon')).groupby('point')
    elif 'latitude' in da_precip.coords:
        da_temp_groupby = da_temp.stack(point=('latitude', 'longitude')).groupby('point')
    elif 'Y' in da_precip.coords:
        da_temp_groupby = da_temp.stack(point=('Y', 'X')).groupby('point')
    else:
        raise KeyError('latitude not found')

    # perform calculation
    da_pet = xr.apply_ufunc(indices.pet,
                            da_temp_groupby,
                            lat,
                            start_year)
    da_pet = da_pet.unstack('point')
    #ds_pet_all = da_pet.to_dataset(name='pet')
    #cut to user define choices
    if data_start_year>start_year:
        min_idx = np.where(da_pet.time.dt.year>=data_start_year)[0][0]
    else:
        min_idx=0
    if data_end_year<end_year:
        max_idx = np.where(da_pet.time.dt.year<=data_end_year)[0][-1]
    else:
        max_idx = da_pet.time.size-1
    r = np.arange(min_idx,max_idx+1,1)
    #Cut along the right dimensions
    da_pet_cut = np.take(da_pet,r,axis=da_pet.dims.index('time'))
    da_pet_cut['time'] = da_pet_cut.time[r]
    ds_pet=da_pet_cut.to_dataset(name='pet')

    info={'index':'PET'}

    return ds_pet, da_pet, info

def SPEI(da_precip, da_temp, distribution = 'gamma', periodicity = 'monthly',\
        scales = 6, data_start_year = 'beginning', data_end_year = 'end',\
        calibration_start_year = 'beginning',\
        calibration_end_year = 'end'):
    """Calculate SPI from precipitation

    This function uses the SPI calculation from the climate indices package

    Args:
        da_precip (Xarray DataArray): A dataArray of precipitation
        da_temp (Xarray DataArray): A datarray of temperature values
        distribution (str): The distribution used to fit the data. Default is
            'gamma'. To use a Pearson Type III distribution, enter 'pearson'
        periodicity (str): Either 'monthly' or 'daily'
        scales (int): The timescales on which the index is computed, either 6 or 12.
            Default is 6.
        data_start_year: Year to start computing  SPI - Default is first year in the data
        data_end_year: Year to stop computing  SPI - Default is first year in the data
        calibration_start_year: Start year for the calibration - Defauls is first year in the data
        calibration_end_year: End year for the calibration - Default is to set a 30-year  period,
            or the full dataset if shorter than 30 years

    Returns:
        ds_spi (Xarray DataSet): SPEI index  DataArray
        info (dict): Dictionary containing relevant information about the calib period
    """

    #Perform some checks
    assert scales==6 or scales==12, 'Valid entries for timescales field should be 6 or 12'
    assert distribution == 'gamma' or distribution == 'pearson', "Valid entries for distribution field should be 'gamma' or pearson'"
    assert periodicity == 'monthly' or periodicity == 'daily', "Valid entries for periodicity field should be 'monthly' or 'daily'"

    if distribution == 'gamma':
        dist = indices.Distribution.gamma
    elif distribution == 'pearson':
        dist = indices.Distribution.pearson

    if periodicity == 'monthly':
        period = compute.Periodicity.monthly
    if periodicity == 'daily':
        period = compute.Periodicity.daily

    if data_start_year == 'beginning':
        data_start_year = int(np.min(da_precip.time.dt.year))
    elif data_start_year not in da_precip.time.dt.year.values:
        print('Start year not in dataset, using first available year')
        data_start_year = int(np.min(da_precip.time.dt.year))

    if data_end_year == 'end':
        data_end_year = int(np.max(da_precip.time.dt.year))
    elif data_end_year not in da_precip.time.dt.year.values:
        print('End year not in dataset, using last available year')
        data_end_year = int(np.max(da_precip.time.dt.year))

    if calibration_start_year == 'beginning':
        calibration_start_year = int(np.min(da_precip.time.dt.year))
    elif calibration_start_year not in da_precip.time.dt.year.values:
        print('Calibration start year not in dataset, using first available year')
        calibration_start_year = int(np.min(da_precip.time.dt.year))

    if calibration_end_year == 'end':
        calibration_end_year = calibration_start_year + 29
        if calibration_end_year not in da_precip.time.dt.year.values:
            calibration_end_year = int(np.max(da_precip.time.dt.year))
    elif calibration_end_year not in da_precip.time.dt.year.values:
        print('Calibration end year not in dataset, using 30yr window or end of dataset')
        calibration_end_year = calibration_start_year + 29
        if calibration_end_year not in da_precip.time.dt.year.values:
            calibration_end_year = int(np.max(da_precip.time.dt.year))

    #compute pet
    ds_pet,da_pet,info = PET(da_temp,data_start_year,\
                             data_end_year)

    #Resize the da_precip array as needed
    if np.min(da_precip.time.dt.year) != np.min(da_pet.time.dt.year):
        min_idx = np.where(da_precip.time.dt.year>=np.min(da_pet.time.dt.year))[0][0]
    else:
        min_idx=0
    if np.max(da_precip.time.dt.year) != np.max(da_pet.time.dt.year):
        max_idx = np.where(da_precip.time.dt.year<=np.max(da_pet.time.dt.year))[0][0]
    else:
        max_idx = da_pet.time.size-1
    r = np.arange(min_idx,max_idx+1,1)
    da_precip_cut = np.take(da_precip,r,axis=da_precip.dims.index('time'))

    #groupby
    if 'lat' in da_temp.coords:
        da_precip_groupby = da_precip_cut.stack(point=('lat', 'lon')).groupby('point')
        da_pet_groupby = da_pet.stack(point=('lat', 'lon')).groupby('point')
    elif 'latitude' in da_precip.coords:
        da_precip_groupby = da_precip_cut.stack(point=('latitude', 'longitude')).groupby('point')
        da_pet_groupby = da_pet.stack(point=('latitude', 'longitude')).groupby('point')
    elif 'Y' in da_precip.coords:
        da_precip_groupby = da_precip_cut.stack(point=('Y', 'X')).groupby('point')
        da_pet_groupby = da_pet.stack(point=('Y', 'X')).groupby('point')
    else:
        raise KeyError('latitude not found')

    #perform the calculation
    da_spei=xr.apply_ufunc(indices.spei,
                          da_precip_groupby,
                          da_pet_groupby,
                          scales,
                          dist,
                          period,
                          data_start_year,
                          calibration_start_year,
                          calibration_end_year)

    #Cut to the right years
    da_spei = da_spei.unstack('point')
    min_idx = np.where(da_spei.time.dt.year>=data_start_year)[0][0]
    max_idx = np.where(da_spei.time.dt.year<=data_end_year)[0][-1]
    r = np.arange(min_idx,max_idx+1,1)
    da_spei_cut = np.take(da_spei,r,axis=da_spei.dims.index('time'))
    da_spei_cut['time'] = da_spei_cut.time[r]
    ds_spei=da_spei_cut.to_dataset(name='spei') #convert to xarrayDataset

    info = {'index':'SPEI',
            'calibration_start': calibration_start_year,
            'calibration_end':  calibration_end_year,
            'distribution': distribution,
            'periodicity': periodicity,
            'timescales': scales}
    return ds_spei, info
#%% Return a netcdf using MINT conventions
def to_netcdfMint(ds, info, dataset_type, bounding_box, dir_out):
    """Returns a MINT-ready netcdf file with SPI values

    Args:
        ds (Xarray DataSet): A dataset of drought indices
        info (dict): Dictionary containing pertinent information frpm calculation
        dir_out (str): The out directory to write the netcdf files

    Returns:
        NetCDF ouput in MINT Format

    """
    #Adding the global attributes

    if os.path.isdir(dir_out+'/results') is False:
        os.makedirs(dir_out+'/results')

    ds.attrs['conventions'] = 'MINT-1.0'
    if info['index'] == 'SPI':
        long_name = 'Standardized Precipitation Index'
        ds.attrs['title'] = long_name
        ds.attrs['summary'] = info['periodicity']+' standardized precipitation index inferred from '\
        + dataset_type +', using a calibration period from '+\
        str(info['calibration_start']) + ' to '+ str(info['calibration_end'])+\
            ' using a '+ info['distribution']+ ' distribution with a '+\
            str(info['timescales'])+ '-month timescale.'
    elif info['index']== 'PET':
        long_name = 'Potential Evapotranspiration'
        ds.attrs['title'] = long_name
        ds.attrs['summary'] = 'Monthly potential evapotranspiration inferred from '\
        + dataset_type +'.'
        info['periodicity']='monthly'
    elif info['index'] == 'SPEI':
        long_name = 'Standardized Precipitation Evapotranspiration Index'
        ds.attrs['title'] = long_name
        ds.attrs['summary'] = info['periodicity']+' standardized precipitation  evapotranspiration index inferred from '\
        + dataset_type +', using a calibration period from '+\
        str(info['calibration_start']) + ' to '+ str(info['calibration_end'])+\
            ' using a '+ info['distribution']+ ' distribution with a '+\
            str(info['timescales'])+ '-month timescale.'
    ds.attrs['naming_authority'] = "MINT Workflow"
    ds.attrs['id'] = str(uuid.uuid4())
    ds.attrs['date_created'] = str(date.today())
    ds.attrs['date_modified']= str(date.today())
    ds.attrs['creator_name'] = 'Deborah Khider'
    ds.attrs['creator_email'] = 'khider@usc.edu'
    ds.attrs['institution'] = 'USC Information Sciences Institute'
    ds.attrs['project'] = 'MINT'
    ds.attrs['time_coverage_start'] = str(ds.time.values[0])
    ds.attrs['time_coverage_end'] = str(ds.time.values[-1])
    ds.attrs['time_coverage_resolution'] = info['periodicity']
    ds.attrs['geospatial_lat_min'] = bounding_box[2]
    ds.attrs['geospatial_lat_max'] = bounding_box[3]
    ds.attrs['geospatial_lon_min'] = bounding_box[0]
    ds.attrs['geospatial_lon_max'] = bounding_box[1]

    # Adding var attributes
    if info['index'] == 'SPI':
       ds.spi.attrs['title'] = 'Standardized Precipitation Index'
       ds.spi.attrs['standard_name'] = 'atmosphere_water__standardized_precipitation_wetness_index'
       ds.spi.attrs['long_name'] = 'Standardized Precipitation Index'
       ds.spi.attrs['units'] = 'unitless'
       ds.spi.attrs['valid_min'] = np.min(ds.spi.values)
       ds.spi.attrs['valid_max'] = np.max(ds.spi.values)
       ds.spi.attrs['valid_range'] = list((np.min(ds.spi.values), np.max(ds.spi.values)))
       ds.spi.attrs['missing_value'] = np.nan
    elif info['index'] == 'PET':
       ds.pet.attrs['title'] = 'Potential Evapotranspiration'
       ds.pet.attrs['standard_name'] = 'atmosphere_soil_water__thornthwaite_potential_evapotranspiration_volume'
       ds.pet.attrs['long_name'] = 'Potential Evapotranspiration'
       ds.pet.attrs['units'] = 'mm/month'
       ds.pet.attrs['valid_min'] = np.min(ds.pet.values)
       ds.pet.attrs['valid_max'] = np.max(ds.pet.values)
       ds.pet.attrs['valid_range'] = list((np.min(ds.pet.values), np.max(ds.pet.values)))
       ds.pet.attrs['missing_value'] = np.nan
    elif info['index'] == 'SPEI':
       ds.spei.attrs['title'] = 'Standardized Precipitation Evapotranspiration Index'
       ds.spei.attrs['standard_name'] = 'land_region_water__standardized_precipitation_evapotranspiration_drought_intensity_index'
       ds.spei.attrs['long_name'] = 'Standardized Precipitation Index'
       ds.spei.attrs['units'] = 'unitless'
       ds.spei.attrs['valid_min'] = np.min(ds.spei.values)
       ds.spei.attrs['valid_max'] = np.max(ds.spei.values)
       ds.spei.attrs['valid_range'] = list((np.min(ds.spei.values), np.max(ds.spei.values)))
       ds.spei.attrs['missing_value'] = np.nan

    #Write it out to file
    if dir_out[-1]=='/':
        path = dir_out+'results/'+info['index']+ '_' + ds.attrs['project'] + '_' +\
                dataset_type+'_'+ds.attrs['time_coverage_start'].split('T')[0]+\
                '_'+\
                ds.attrs['time_coverage_end'].split('T')[0]+\
                '_'+ds.attrs['id']+'.nc'
    else:
        path = dir_out+'/results/'+info['index']+ '_' + ds.attrs['project'] + '_' +\
                dataset_type+'_'+ds.attrs['time_coverage_start'].split('T')[0]+\
                '_'+\
                ds.attrs['time_coverage_end'].split('T')[0]+\
                '_'+ds.attrs['id']+'.nc'
    ds.to_netcdf(path = path)

def visualizeDroughtIndex(ds, dir_out, info, dataset_type):
    """ Visualization of drought index

    Args:
        ds (xarray dataset): the dataset containing the index
        dir_out (str): the output directory for the visualization
    """
    proj = ccrs.PlateCarree()
    idx = np.size(ds['time'])
    count = list(np.arange(0,idx,1))
    varname = list(ds.data_vars.keys())[0]
    filenames =[]

    #Make a directory for results/figures if it doesn't exit
    if dir_out[-1]=='/':
        if os.path.isdir(dir_out+'figures') is False:
            os.makedirs(dir_out+'figures')
        if os.path.isdir(dir_out+'results') is False:
            os.makedirs(dir_out+'results')
    else:
        if os.path.isdir(dir_out+'/figures') is False:
            os.makedirs(dir_out+'/figures')
        if os.path.isdir(dir_out+'/results') is False:
            os.makedirs(dir_out+'/results')

    if 'lat' in ds:
        lat = ds.lat.values
        lon = ds.lon.values
    elif 'latitude' in ds:
        lat = ds.latitude.values
        lon = ds.longitude.values
    elif 'Y' in ds:
        lat = ds.Y.values
        lon = ds.X.values
    else:
        raise KeyError('latitude not found')

    for i in count:
        v = ds[varname].values[i,:,:]
        date = pd.to_datetime(ds.time.values[i]).strftime("%B %Y")
        fig,ax = plt.subplots(figsize=[15,10])
        ax = plt.axes(projection=proj)
        ax.add_feature(cfeature.BORDERS)
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.RIVERS)
        levels = np.arange(-4,4.2,0.2)
        if np.isnan(v).all()==False:
            img = plt.contourf(lon, lat, v, levels=levels,
                cmap=cm.BrBG,
                transform=proj,
                vmin = -4,
                vmax =4)
            tick_range = np.arange(-4,4.5,0.5)
            cbar = plt.colorbar(img, orientation='horizontal',pad=0.1, ticks=tick_range)
            string = ', inferred from' + dataset_type +\
                    ' over a calibration period from '+ str(info['calibration_start']) +\
                    ' to '+ str(info['calibration_end']) + ', using a '+info['distribution']+\
                    ' distribution and '+str(info['timescales'])+'-month timescale.'
            if varname == 'spi':
                cbar.ax.set_title('Standardized Precipitation Index'+string)
            elif varname == 'pet':
                cbar.ax.set_title('Potential Evapotranspiration (mm/month)'+string)
            elif varname == 'spei':
                cbar.ax.set_title('Standardized Precipitation-Evapotranspiration Index'+string)
            else:
                cbar.ax.set_title(varname)
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xlines = False
        gl.ylines = False
        gl.xlocator = mticker.FixedLocator(np.linspace(np.round(np.min(lon)),np.round(np.max(lon)),5))
        gl.ylocator = mticker.FixedLocator(np.linspace(np.round(np.min(lat)),np.round(np.max(lat)),5))
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 12, 'color': 'gray'}
        gl.ylabel_style = {'size': 12, 'color': 'gray'}
        #Add a title with time
        plt.title(date, fontsize=18, loc='left', pad=1)
        #save a jepg
        if dir_out[-1]=='/':
            filename = dir_out+'figures/'+varname+'_t'+str(date)+'.png'
        else:
            filename = dir_out+'/figures/'+varname+'_t'+str(date)+'.png'
        filenames.append(filename)
        plt.savefig(filename)
        plt.close(fig)

    #create a gif
    if dir_out[-1]=='/':
        writer = imageio.get_writer(dir_out+'results/'+varname+'_'+ds.attrs['id']+'.mp4', fps=5)
    else:
        writer = imageio.get_writer(dir_out+'/results/'+varname+'_'+ds.attrs['id']+'.mp4', fps=5)

    for filename in filenames:
        writer.append_data(imageio.imread(filename))
    writer.close()



#%% Main
if __name__ == "__main__":
    dataset_type = sys.argv[1]
    index = sys.argv[4]
    dataset_name = sys.argv[2] #file name or directory name
    bounding_box = ast.literal_eval(sys.argv[5])
    distribution = sys.argv[6].lower()
    periodicity = sys.argv[7].lower()
    scales = int(sys.argv[8])
    try:
        data_start_year= int(sys.argv[9])
    except:
        data_start_year = sys.argv[9]
    try:
        data_end_year = int(sys.argv[10])
    except:
        data_end_year = sys.argv[10]
    try:
        calibration_start_year = int(sys.argv[11])
    except:
        calibration_start_year = sys.argv[11]
    try:
        calibration_end_year = int(sys.argv[12])
    except:
        calibration_end_year = sys.argv[12]
    dir_out = sys.argv[3]
    fig = ast.literal_eval(sys.argv[13])

    #Test
#    dataset_type = 'GLDAS'
#    dataset_name = 'GLDAS2.0_TP_1948_2010.nc'
#    dir_out = '/Users/deborahkhider/Desktop/'
#    index = 'SPI'
#    bounding_box = [23,48,3,15]
#    distribution = 'gamma'
#    periodicity = 'monthly'
#    scales = 6
#    data_start_year= 2000
#    data_end_year = 2010
#    calibration_start_year = 1981
#    calibration_end_year = 2010
#    fig = True

    ## Perform checks:
    #possible datatypes
    dataset_type_list = ['CHIRPS','GLDAS','FLDAS']
    if dataset_type not in dataset_type_list:
        raise ValueError("Dataset type not a valid entry. Use either 'CHIRPS', 'GLDAS', 'FLDAS'")
    #possible indices
    index_list = ['SPI', 'PET', 'SPEI']
    if index not in index_list:
        raise ValueError("Index is not a valid entry. Enter either 'SPI', 'PET' or 'SPEI'")

    #Some combinations are invalid
    if dataset_type == 'CHIRPS':
        if index in ['PET','SPEI']:
            raise ValueError("Index calculation not supported for this dataset.")
    #beginning and end years
    if calibration_start_year != 'beginning':
        calibration_start_year = int(calibration_start_year)
    if calibration_end_year != 'end':
        calibration_end_year = int(calibration_end_year)
    if data_start_year != 'beginning':
        data_start_year = int(data_start_year)
    if calibration_start_year>calibration_end_year:
        raise ValueError('The beginning of the calibration period should be set prior to the end')
    if data_start_year>data_end_year:
        raise ValueError('The beginning of the simulationperiod should be set prior to the end')
    #periodicity
    if periodicity !='monthly':
        raise ValueError("Periodicity parameter should be set to monthly")
    #Scales
    if scales not in [6,12]:
        raise ValueError('Scales should be 6 or 12')
    #distirbutions
    dist_list=['gamma','pearson']
    if distribution not in dist_list:
        raise ValueError("Valid distriubtion is 'gamma' or 'pearson'")
    ## Open datasets
    if dataset_type == 'CHIRPS':
        da_precip = openCHIRPS(dataset_name, bounding_box)
    elif dataset_type == 'GLDAS':
        if dataset_name.endswith('.nc')==True:
            da_precip,da_temp = openGLDAS(dataset_name, bounding_box, periodicity, True)
        else:
            da_precip,da_temp = openGLDAS(dataset_name, bounding_box, periodicity, False)
    elif dataset_type == 'FLDAS':
        if dataset_name.endswith('.nc')==True:
            da_precip,da_temp = openFLDAS(dataset_name, bounding_box, periodicity, True)
        else:
            da_precip,da_temp = openFLDAS(dataset_name, bounding_box, periodicity, False)
    ## Perform calcuculations
    if index == 'SPI':
        ds, info = SPI(da_precip, distribution, periodicity, scales,\
                           data_start_year, data_end_year, calibration_start_year, calibration_end_year)
    elif index == 'PET':
        ds, da_pet, info = PET(da_temp,data_start_year,data_end_year)
    elif index == 'SPEI':
        ds, info = SPEI(da_precip, da_temp,distribution, periodicity, scales,\
                            data_start_year, data_end_year, calibration_start_year, calibration_end_year)
    ## Write to file
    to_netcdfMint(ds, info, dataset_type, bounding_box, dir_out)
    ## Do vizualization if asked
    if fig == True:
        if index == 'PET':
            print("Visualization is not supported")
        else:
            visualizeDroughtIndex(ds, dir_out, info, dataset_type)
