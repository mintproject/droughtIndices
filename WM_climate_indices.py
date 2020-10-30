
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
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import json
import warnings
import matplotlib.cm as cm

#%% Open files from various agencies
## CHIRPS
def openCHIRPS(dataset_name, bounding_box, globe):
    """ Open CHIRPS dataset and returns the data

    Args:
        dataset_name (str): The name of the CHIRPS dataset
        bounding_box (list): lat/lon to cut to appropriate size
        globe (bool): If considering the full spatial coverage

    Returns:
        da_precip (Xarray DataArray): A dataArray of precipitation

    """
    
    data = xr.open_dataset(dataset_name)
    if globe == True:
        da_precip = data.precip
        warnings.warn('Global calculations require a few hours to complete',
                      UserWarning, stacklevel=2)
    else:    
        p_ = data.sel(latitude=slice(bounding_box[2], bounding_box[3]),\
                  longitude=slice(bounding_box[0],bounding_box[1]))
        da_precip = p_.precip

    return da_precip

def openGLDAS(dataset_name, bounding_box, globe, periodicity, netcdf):
    """Open GLDAS datasets and return precipitation and temperature

    Args:
        dataset_name (str): The name of the GLDAS folder
        bounding_box (list): lat/lon to cut to appropriate size
        globe (bool): If considering the full spatial coverage
        periodicity (str): The temporal resolution of the input data. Useful to
            know how the data is organized.
        netcdf (bool): Whether the input data is in netCDF format.

    Returns:
        da_precip (Xarra DataArray): A dataArray of precipitation
        da_temp (Xarra DataArray): A dataArray of temperature
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
    
    if globe == True:
        da_precip = data.Rainf_f_tavg
        da_temp = data.Tair_f_inst
        warnings.warn('Global calculations require a few hours to complete',
                      UserWarning, stacklevel=2)
    else:
        p_ = data.sel(lat=slice(bounding_box[2], bounding_box[3]),\
                  lon=slice(bounding_box[0],bounding_box[1]))
        da_precip = p_.Rainf_f_tavg
        da_temp = p_.Tair_f_inst
    ##unit transformation from kgm-2s-1 to mm
    #days in the month
    time = pd.to_datetime(da_precip.time.values)
    v = time.strftime("%Y,%m")
    days = []
    for item in v:
        days.append(monthrange(int(item.split(',')[0]),int(item.split(',')[1]))[1])
    da_precip = da_precip*86400
    da_precip = da_precip.T*days
    da_precip = da_precip.T
    da_precip.attrs['units'] = 'mm'
    ### TEMPERATURE
    ## Unit conversion
    da_temp = da_temp-273.15
    da_temp.attrs['units'] = 'C'
    
    return da_precip, da_temp

def openFLDAS(dataset_name, bounding_box, globe, periodicity, netcdf):
    """Open FLDAS datasets and return precipitation and temperature

    Args:
        dataset_name (str): The name of the GLDAS folder
        bounding_box (list): lat/lon to cut to appropriate size
        globe (bool): If considering the full spatial coverage
        periodicity (str): The temporal resolution of the input data. Useful to
            know how the data is organized.
        netcdf (bool): Whether the input data is in netCDF format.

    Returns:
        da_precip (Xarra DataArray): A dataArray of precipitation
        da_temp (Xarra DataArray): A dataArray of temperature
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
    if globe == True:
        da_precip = data.Rainf_f_tavg
        da_temp = data.Tair_f_inst
        warnings.warn('Global calculations require a few hours to complete',
                      UserWarning, stacklevel=2)
    else:
        p_ = data.sel(Y=slice(bounding_box[2], bounding_box[3]),\
                  X=slice(bounding_box[0],bounding_box[1]))
        da_precip = p_.Rainf_f_tavg
        da_temp = p_.Tair_f_tavg
    ##unit transformation from kgm-2s-1 to mm
    #days in the month
    time = pd.to_datetime(da_precip.time.values)
    v= time.strftime("%Y,%m")
    days=[]
    for item in v:
        days.append(monthrange(int(item.split(',')[0]),int(item.split(',')[1]))[1])
    da_precip = da_precip*86400
    da_precip = da_precip.T*days
    da_precip = da_precip.T
    da_precip.attrs['units'] = 'mm'
    ### TEMPERATURE
    ## Unit conversion
    da_temp = da_temp-273.15
    da_temp.attrs['units'] = 'C'

    return da_precip, da_temp

def openECMWF(dataset_name, bounding_box, globe):
    """Open ECMWF datasets and return precipitation and temperature

    Args:
        dataset_name (str): The name of the ECMWF file
        bounding_box (list): lat/lon to cut to appropriate size
        globe (bool): If considering the full spatial coverage
        
    Returns:
        da_precip (Xarra DataArray): A dataArray of precipitation
        da_temp (Xarra DataArray): A dataArray of temperature
    """

    # Loop over all datasets in the various folders to get the data
    data = xr.open_dataset(dataset_name)
    data = data.reindex(latitude=list(reversed(data.latitude)))
    if globe == True:
        da_precip = data.tp
        da_temp = data.t2m
        warnings.warn('Global calculations require a few hours to complete',
                      UserWarning, stacklevel=2)
    else:
        p_ = data.sel(latitude=slice(bounding_box[2], bounding_box[3]),\
                  longitude=slice(bounding_box[0],bounding_box[1]))
        da_precip = p_.tp
        da_temp = p_.t2m
    ##unit transformation m to mm
    da_precip = da_precip*1000
    da_precip.attrs['units'] = 'mm'
    ### TEMPERATURE
    ## Unit conversion
    da_temp = da_temp-273.15
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
    
    # Load the data
    if len(da_precip.dims) == 4: #Deal with ECMWF data
        da_precip = da_precip[:,0,:,:]
        da_precip.squeeze()
    da_precip.load()
    #Groupby
    if 'lat' in da_precip.coords:
        da_precip_groupby = da_precip.stack(point=('lat', 'lon')).groupby('point')
    elif 'latitude' in da_precip.coords:
        da_precip_groupby = da_precip.stack(point=('latitude', 'longitude')).groupby('point')
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
    t_start  = str(data_start_year)+'-01-01T00:00:00.000000000'
    t_end  = str(data_end_year)+'-12-01T00:00:00.000000000'
    da_spi_cut = da_spi.sel(time=slice(t_start,t_end))
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
        data_start_year = int(np.min(da_temp.time.dt.year))
    elif data_start_year not in da_temp.time.dt.year.values:
        print('Start year not in dataset, using first available year')
        data_start_year = int(np.min(da_temp.time.dt.year))

    if data_end_year == 'end':
        data_end_year = int(np.max(da_temp.time.dt.year))
    elif data_end_year not in da_temp.time.dt.year.values:
        print('End year not in dataset, using last available year')
        data_end_year = int(np.max(da_temp.time.dt.year))

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
    if int(da_temp.time[0].dt.month) != 1:
        print("Full year not available for the beginning of the record, truncating...")
        start_year = int(da_temp.time.dt.year[1])
    else:
        start_year = int(da_temp.time.dt.year[0])

    if int(da_temp.time[-1].dt.month) != 12:
        print("Full year not available for the end of the record, truncating...")
        end_year = int(da_temp.time.dt.year[-2])
    else:
        end_year = int(da_temp.time.dt.year[-1])

    #cut the data if entire years are not available.
    t_start  = str(start_year)+'-01-01T00:00:00.000000000'
    t_end =  str(end_year)+'-12-01T00:00:00.000000000'
    da_temp_cut = da_temp.sel(time=slice(t_start,t_end))

    if len(da_temp_cut.dims) == 4: #Deal with ECMWF data
        da_temp_cut = da_temp_cut[:,0,:,:]
    da_temp_cut.load()

    #Groupby
    if 'lat' in da_temp.coords:
        da_temp_groupby = da_temp_cut.stack(point=('lat', 'lon')).groupby('point')
    elif 'latitude' in da_precip.coords:
        da_temp_groupby = da_temp_cut.stack(point=('latitude', 'longitude')).groupby('point')
    elif 'Y' in da_precip.coords:
        da_temp_groupby = da_temp_cut.stack(point=('Y', 'X')).groupby('point')
    else:
        raise KeyError('latitude not found')
    # Load the data
    
    # perform calculation
    da_pet = xr.apply_ufunc(indices.pet,
                            da_temp_groupby,
                            lat,
                            start_year)
    da_pet_un = da_pet.unstack('point')
    if data_start_year>start_year:
        t_start= str(data_start_year)+'-01-01T00:00:00.000000000'
    if data_end_year<end_year:
        t_end  = str(data_end_year)+'-12-01T00:00:00.000000000'
    da_pet = da_pet_un.sel(time=slice(t_start,t_end))
    ds_pet=da_pet.to_dataset(name='pet')

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
    min_year = np.max([int(np.min(da_precip.time.dt.year)),int(np.min(da_pet.time.dt.year))])
    t_start = str(min_year)+'-01-01T00:00:00.000000000'
    max_year = np.min([int(np.max(da_precip.time.dt.year)),int(np.max(da_pet.time.dt.year))])
    t_end = str(max_year)+'-12-01T00:00:00.000000000'
    da_precip_cut = da_precip.sel(time=slice(t_start,t_end))

    #load the data
    if len(da_precip_cut.dims) == 4: #Deal with ECMWF data
        da_precip_cut = da_precip_cut[:,0,:,:]
    da_precip_cut.load()

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
    t_start  = str(data_start_year)+'-01-01T00:00:00.000000000'
    t_end  = str(data_end_year)+'-12-01T00:00:00.000000000'
    da_spei_cut = da_spei.sel(time=slice(t_start,t_end))
    ds_spei=da_spei_cut.to_dataset(name='spei') #convert to xarrayDataset

    info = {'index':'SPEI',
            'calibration_start': calibration_start_year,
            'calibration_end':  calibration_end_year,
            'distribution': distribution,
            'periodicity': periodicity,
            'timescales': scales}
    return ds_spei, info
#%% Return a netcdf using MINT conventions
def to_netcdfMint(ds, info, dataset_type, dir_out, dynamic_name):
    """Returns a MINT-ready netcdf file with SPI values

    Args:
        ds (Xarray DataSet): A dataset of drought indices
        info (dict): Dictionary containing pertinent information frpm calculation
        dir_out (str): The out directory to write the netcdf files
        dynamic_name (bool): Whether to generate a unique name dynamically using parameter settings. Otherwise returns, results.nc

    Returns:
        ds (xarray dataset): Metadata complete dataset (important for viz)
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
    if 'latitude' in ds.coords:
        ds.attrs['geospatial_lat_min'] = float(ds.latitude.min())
        ds.attrs['geospatial_lat_max'] = float(ds.latitude.max())
        ds.attrs['geospatial_lon_min'] = float(ds.longitude.min())
        ds.attrs['geospatial_lon_max'] = float(ds.longitude.max())
    elif 'lat' in ds.coords:
        ds.attrs['geospatial_lat_min'] = float(ds.lat.min())
        ds.attrs['geospatial_lat_max'] = float(ds.lat.max())
        ds.attrs['geospatial_lon_min'] = float(ds.lon.min())
        ds.attrs['geospatial_lon_max'] = float(ds.lon.max())
    elif 'Y' in ds.coords:
        ds.attrs['geospatial_lat_min'] = float(ds.Y.min())
        ds.attrs['geospatial_lat_max'] = float(ds.Y.max())
        ds.attrs['geospatial_lon_min'] = float(ds.X.min())
        ds.attrs['geospatial_lon_max'] = float(ds.X.max())
    

    # Adding var attributes
    if info['index'] == 'SPI':
       ds.spi.attrs['title'] = 'Standardized Precipitation Index'
       ds.spi.attrs['standard_name'] = 'atmosphere_water__standardized_precipitation_wetness_index'
       ds.spi.attrs['long_name'] = 'Standardized Precipitation Index'
       ds.spi.attrs['units'] = 'unitless'
       ds.spi.attrs['valid_min'] = float(ds.spi.min())
       ds.spi.attrs['valid_max'] = float(ds.spi.max())
       ds.spi.attrs['valid_range'] = list((ds.spi.attrs['valid_min'], ds.spi.attrs['valid_max']))
       ds.spi.attrs['missing_value'] = np.nan
    elif info['index'] == 'PET':
       ds.pet.attrs['title'] = 'Potential Evapotranspiration'
       ds.pet.attrs['standard_name'] = 'atmosphere_soil_water__thornthwaite_potential_evapotranspiration_volume'
       ds.pet.attrs['long_name'] = 'Potential Evapotranspiration'
       ds.pet.attrs['units'] = 'mm'
       ds.pet.attrs['valid_min'] = float(ds.pet.min())
       ds.pet.attrs['valid_max'] = float(ds.pet.max())
       ds.pet.attrs['valid_range'] = list((ds.pet.attrs['valid_min'], ds.pet.attrs['valid_max']))
       ds.pet.attrs['missing_value'] = np.nan
    elif info['index'] == 'SPEI':
       ds.spei.attrs['title'] = 'Standardized Precipitation Evapotranspiration Index'
       ds.spei.attrs['standard_name'] = 'land_region_water__standardized_precipitation_evapotranspiration_drought_intensity_index'
       ds.spei.attrs['long_name'] = 'Standardized Precipitation Index'
       ds.spei.attrs['units'] = 'unitless'
       ds.spei.attrs['valid_min'] = float(ds.spei.min())
       ds.spei.attrs['valid_max'] = float(ds.spei.max())
       ds.spei.attrs['valid_range'] = list((ds.spei.attrs['valid_min'], ds.spei.attrs['valid_max']))
       ds.spei.attrs['missing_value'] = np.nan

    #Write it out to file
    if dynamic_name == True:
        f = info['index']+ '_' + ds.attrs['project'] + '_' +\
                dataset_type+'_'+ds.attrs['time_coverage_start'].split('T')[0]+\
                '_'+\
                ds.attrs['time_coverage_end'].split('T')[0]+\
                '_'+ds.attrs['id']+'.nc'
    else:
        f = 'results.nc'
    if dir_out[-1]=='/':
        path = dir_out+'results/'+f
    else:
        path = dir_out+'/results/'+f
    ds.to_netcdf(path = path)
    
    return ds

#%% Visualization 
    
def visualizeDroughtIndex(ds, dir_out, dynamic_name):
    """ Visualization of drought index

    Args:
        ds (xarray dataset): the dataset containing the index
        dir_out (str): the output directory for the visualization
        dynamic_name (bool): Whether to generate a unique name dynamically using parameter settings. Otherwise returns, results.mp4
    """
    proj=ccrs.PlateCarree()
    idx = np.size(ds['time'])
    count = list(np.arange(0,idx,1))
    varname=list(ds.data_vars.keys())[0]
    filenames=[]

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
            
    # Get the levels for the countours
    if varname == 'spi' or varname == 'spei':
        levels = np.arange(-4,4.2,0.2)
    else:
        levels = np.linspace(float(ds[varname].min()),float(ds[varname].max()),40)
    for i in count:
        fig,ax = plt.subplots(figsize=[15,10])
        ax = plt.axes(projection=proj)
        if xr.apply_ufunc(np.isnan,ds[varname].isel(time=i)).all()==False:
            if varname == 'spi' or varname == 'spei':
                ds[varname].isel(time=i).plot.contourf(ax=ax,levels = levels, 
                  transform=ccrs.PlateCarree(), cmap=cm.BrBG, cbar_kwargs={'orientation':'horizontal'})
            else: 
                ds[varname].isel(time=i).plot.contourf(ax=ax,levels = levels, 
                  transform=ccrs.PlateCarree(), cmap=cm.viridis, cbar_kwargs={'orientation':'horizontal'})                
        ax.add_feature(cfeature.BORDERS)
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.RIVERS)
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
                      linewidth=2, color='gray', alpha=0.5, linestyle='--')
        gl.xlabels_top = False
        gl.ylabels_right = False
        gl.xlines = False
        gl.ylines = False
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.xlabel_style = {'size': 12, 'color': 'gray'}
        gl.ylabel_style = {'size': 12, 'color': 'gray'}
        #save a jepg
        
        if dir_out[-1]=='/':
            filename = dir_out+'figures/'+varname+'_t'+str(i)+'.jpeg'
        else:
            filename = dir_out+'/figures/'+varname+'_t'+str(i)+'.jpeg'
        filenames.append(filename)
        plt.savefig(filename)
        plt.close(fig)

    #create a gif
    if dynamic_name == True:
        f = varname+'_'+ds.attrs['id']+'.mp4'
    else:
        f = 'results.mp4'
    if dir_out[-1]=='/':
        writer = imageio.get_writer(dir_out+'results/'+f, fps=5)
    else:
        writer = imageio.get_writer(dir_out+'/results/'+f, fps=5)

    for filename in filenames:
        writer.append_data(imageio.imread(filename))
    writer.close()

#%% Main
if __name__ == "__main__":    
    #Open JSON file and get settings
    with open(sys.argv[1]) as json_file:
        config = json.load(json_file)
    if ast.literal_eval(config['debug']) == False:
        dataset_type = config['data']['dataset_type']
        index = config['index']['name']
        dataset_name = config['data']['dataset_name'] #file name or directory name
        globe = ast.literal_eval(config['spatial']['global'])
        if globe == False:
            bounding_box = ast.literal_eval(config['spatial']['bounding_box'])
        elif globe == True:
            bounding_box = []
        else:
            raise ValueError("globe option should be set as 'True' or 'False'")
        distribution = config['index']['distribution'].lower()
        periodicity = config['index']['periodicity'].lower()
        scales = int(config['index']['scales'].lower())
        try:
            data_start_year= int(config['index']['data_start_year'])
        except:
            data_start_year = config['index']['data_start_year']
        try:
            data_end_year = int(config['index']['data_end_year'])
        except:
            data_end_year = config['index']['data_end_year']
        try:
            calibration_start_year = int(config['index']['calibration_start_year'])
        except:
            calibration_start_year = config['index']['calibration_start_year']
        try:
            calibration_end_year = int(config['index']['calibration_end_year'])
        except:
            calibration_end_year = config['index']['calibration_end_year']
        dynamic_name = ast.literal_eval(config['output']['dynamic_name'])
        dir_out = config['output']['path']
        fig = ast.literal_eval(config['output']['fig'])
    elif ast.literal_eval(config['debug']) == 'True':
        dataset_type = 'GLDAS'
        dataset_name = 'GLDAS2.0_TP_1948_2010.nc'
        dynamic_name = 'True'
        dir_out = './'
        index = 'SPI'
        bounding_box = [23,48,3,15]
        globe = 'False'
        distribution = 'gamma'
        periodicity = 'monthly'
        scales = 6
        data_start_year= 2000
        data_end_year = 2010
        calibration_start_year = 1981
        calibration_end_year = 2010
        fig = True
    else:
        raise ValueError("debug options should be set as 'True' or 'False'")

    ## Perform checks:
    #possible datatypes
    dataset_type_list = ['CHIRPS','GLDAS','FLDAS','ECMWF']
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
    #distributions
    dist_list=['gamma','pearson']
    if distribution not in dist_list:
        raise ValueError("Valid distriubtion is 'gamma' or 'pearson'")
    ## Open datasets
    if dataset_type == 'CHIRPS':
        da_precip = openCHIRPS(dataset_name, bounding_box, globe)
    elif dataset_type == 'ECMWF':
        da_precip, da_temp = openECMWF(dataset_name, bounding_box, globe)
    elif dataset_type == 'GLDAS':
        if dataset_name.endswith('.nc')==True:
            da_precip,da_temp = openGLDAS(dataset_name, bounding_box, globe, periodicity, True)
        else:
            da_precip,da_temp = openGLDAS(dataset_name, bounding_box, globe, periodicity, False)
    elif dataset_type == 'FLDAS':
        if dataset_name.endswith('.nc')==True:
            da_precip,da_temp = openFLDAS(dataset_name, bounding_box, globe, periodicity, True)
        else:
            da_precip,da_temp = openFLDAS(dataset_name, bounding_box, globe, periodicity, False)
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
    ds = to_netcdfMint(ds, info, dataset_type, dir_out,dynamic_name)
    ## Do vizualization if asked
    if fig == True:
        visualizeDroughtIndex(ds, dir_out,dynamic_name)
