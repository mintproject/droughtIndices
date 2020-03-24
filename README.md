[![PyPI](https://img.shields.io/badge/python-3.7-yellow.svg)]()
[![license](https://img.shields.io/github/license/mintproject/droughtIndices.svg)]()

# droughtIndices
Model to calculate various drought indices from precipitation/temperature data.

**Table of contents**
* [What is it?](#what)
* [Version Information](#version)
* [User Guide](#quickstart)
* [Requirements](#req)
* [Files and folders in this repository](#files)
* [Contact](#contact)
* [License](#license)

## <a name = "what">What is it?</a>

This Python routine calculates three climate indices based on NetCDF input data from [CHIRPS](https://www.chc.ucsb.edu/data/chirps), [GLDAS](https://ldas.gsfc.nasa.gov/gldas), [FLDAS](https://ldas.gsfc.nasa.gov/fldas), [ECMWF](https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5):
* [SPI](https://climatedataguide.ucar.edu/climate-data/standardized-precipitation-index-spi): Standardized Precipitation Index, utilizing both gamma and Pearson Type III distributions. - CHIRPS, GLDAS, FLDAS, ECMWF
* [SPEI](https://www.researchgate.net/publication/252361460_The_Standardized_Precipitation-Evapotranspiration_Index_SPEI_a_multiscalar_drought_index): Standardized Precipitation Evapotranspiration Index, utilizing both gamma and Pearson Type III distributions. - GLDAS, FLDAS, ECMWF
* [PET](https://www.ncdc.noaa.gov/monitoring-references/dyk/potential-evapotranspiration): Potential Evapotranspiration, utilizing Thornthwaite equation. - GLDAS, FLDAS, ECMWF

This package is based upon [Climate Indices in Python](https://github.com/monocongo/climate_indices). The Climate Indices in Python project contains implementations of various algorithms which provide a geographical and temporal picture of the severity of precipitation and temperature anomalies for climate monitoring and research.

This project builds upon this software package to include additional sources of climate data.

## <a name = "version">Version information</a>
* v1.1.0: Support for datasets from ECMWF, use JSON for configuration, optimization using xarray
* v1.0.0: Support for datasets from FLDAS and combined netcdf file. Visualization support for SPI and SPEI. - First release
* v0.0.1: Support for three indices and datasets from CHIRPS and GLDAS

## <a name = "quickstart">User Guide</a>

Command line implementations:

`python WM_climate_indices.py config.json`

The config.json file contains the following fields:
* dataset_type: the type of NetCDF files. Valid entries are GLDAS, FLDAS or CHIRPS  
* dataset_name: Name of dataset file or parent folder
* path: The out directory to write the index files. Index files are in NetCDF format, visualization in mp4. All years contained in a single file. Two folders will be created: figures (which contains each frame for the movie) and results, which will include both the netcdf and mp4 files. Each dataset is given a unique ID, which can be used to match data and visualization.
* fig (bool): If True, will generate a movie.
* name: The name of the index to be calculated. Valid entries include:
  * SPI: Standardized Precipitation Index, utilizing both gamma and Pearson Type III distributions
  * PET: Potential Evapotranspiration, utilizing Thornthwaite equations
  * SPEI: Standardized Precipitation Evapotranspiration Index, utilizing both gamma and Pearson Type III distributions.
* distribution: Either 'gamma' or 'pearson'
* periodicity: 'monthly'. Daily calculations are not supported in the current version
* scales: Either '6' or '12' month
* data_start_year: The year for which to start calculating the Index
* data_end_year: The year for which to end the calculation
* calibration_start_year: The start year for the calibration
* calibration_end_year: The end year for the calibration.
* global: Bool to consider the full geographical extent of the dataset. bounding_box information will be ignored in this case.
* bounding_box: List of [min_lon,max_lon,min_lat,max_lat]

**Note**:
- A list of sample data and associated links can be in the data folder in this repository
- SPEI and PET require temperature inputs, not present in the CHIRPS dataset.
- The recommended length of time for calibration is 30 years
- PET calculation requires an entire year of data for  both  calibration and computation. If  not present, will adjust to last/first complete year on record. 
- An example config file is given in this repository

## <a name = "req">Requirements</a>
Tested under Python 3.7

Package requirements:
* xarray - 0.14.1
* dask - 2.2.0
* netcdf4 - 1.4.2
* bottleneck - 1.2.1
* pandas - 0.25.3
* numpy - 1.17.4
* climate-indices - 1.0.6 Install this package from source (not pip) -> https://github.com/monocongo/climate_indices
* cftime - 1.0.4.2
* cartopy - 0.17.0 Install using conda `conda install -c conda-forge cartopy`
* matplotlib - 3.1.0
* imageio - 2.6.1
* imageio-ffmpeg - 0.3.0

## <a name = "files">Files and folders in this directory</a>

### Files

* WM_climate_indices.py: Executable (climate indices calculation)
* config.json: JSON file containing configurable inputs data and parameters

### Folders

* Data: Link to relevant datasets for use with  the software
* Docs: Simulation Matrix for the mint mint project
* ExampleResults: Examples of outputs
* Utils: Utilities to cut and concatenate netcdf files from FLDAS and GLDAS

## <a name = "contact">Contact</a>

Please report issues to <khider@usc.edu>

## <a name ="license"> License </a>

The project is licensed under the Apache v2.0 License. Please refer to the file call license.
