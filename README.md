[![PyPI](https://img.shields.io/badge/python-3.7-yellow.svg)]()
[![license](https://img.shields.io/github/license/mintproject/droughtIndices.svg)]()

# droughtIndices
Model to calculate various drought indices from precipitation/temperature data.

**Table of contents**
* [What is it?](#what)
* [Version Information](#version)
* [User Guide](#quickstart)
* [Requirements](#req)
* [Files in this repository](#files)
* [Contact](#contact)
* [License](#license)

## <a name = "what">What is it?</a>

This Python routine calculates three climate indices based on NetCDF input data from [CHIRPS](https://www.chc.ucsb.edu/data/chirps) and [GLDAS](https://ldas.gsfc.nasa.gov/gldas):
* [SPI](https://climatedataguide.ucar.edu/climate-data/standardized-precipitation-index-spi): Standardized Precipitation Index, utilizing both gamma and Pearson Type III distributions. - CHIRPS and GLDAS
* [SPEI](https://www.researchgate.net/publication/252361460_The_Standardized_Precipitation-Evapotranspiration_Index_SPEI_a_multiscalar_drought_index): Standardized Precipitation Evapotranspiration Index, utilizing both gamma and Pearson Type III distributions. - GLDAS
* [PET](https://www.ncdc.noaa.gov/monitoring-references/dyk/potential-evapotranspiration): Potential Evapotranspiration, utilizing Thornthwaite equation. - GLDAS

This package is based upon [Climate Indices in Python](https://github.com/monocongo/climate_indices). The Climate Indices in Python project contains implementations of various algorithms which provide a geographical and temporal picture of the severity of precipitation and temperature anomalies for climate monitoring and research.

This project builds upon this software package to include additional sources of climate data, namely GLDAS and CHIRPS.

## <a name = "version">Version information</a>
* v0.0.1: Support for three indices and datasets from CHIRPS and GLDAS

## <a name = "quickstart">User Guide</a>

Command line implementations:

`python WM_climate_indices.py dataset_type index dataset_name bounding_box distribution periodicity scales data_start_year data_end_year calibration_start_year calibration_end_year dir_out`

where:
* dataset_type: the type of NetCDF files. Valid entries are GLDASv2.0, GLDASv2.1 or CHIRPS  
* Index: The index to be calculated.
  * SPI: Standardized Precipitation Index, utilizing both gamma and Pearson Type III distributions
  * PET: Potential Evapotranspiration, utilizing Thornthwaite equations
  * SPEI: Standardized Precipitation Evapotranspiration Index, utilizing both gamma and Pearson Type III distributions.
* dataset_name: Name of dataset file or parent folder
* bounding_box: List of [min_lon,max_lon,min_lat,max_lat]
* distribution: Either 'gamma' or 'pearson'
* periodicity: 'monthly'. Daily calculations are not supported in the current version
* scales: Either '6' or '12' month
* data_start_year: The year for which to start calculating the Index
* data_end_year: The year for which to end the calculation
* calibration_start_year: The start year for the calibration
* calibration_end_year: The end year for the calibration
* dir_out: The out directory to write the index files. Index files are in NetCDF format. All years contained in a single file.

**Note**:
- Assumes that the CHIRPS data is contained in a single NetCDF file. GLDAS data are organized in folders, with one month per file.
- SPEI and PET require temperature inputs, not present in the CHIRPS dataset.
- The recommended length of time for calibration is 30 years

Examples:
1. SPI from CHIRPS with a gamma distribution
`WM_climate_indices.py CHIRPS SPI chirps-v2.0.monthly.nc [23,48,3,15] gamma monthly 6 2000 2018 1981 2010 ./results/`

2. SPI from CHIRPS with a pearson distribution
`python WM_climate_indices.py CHIRPS SPI chirps-v2.0.monthly.nc [23,48,3,15] pearson monthly 6 2010 2012 beginning end ./results/`

3. SPI from GLDAS:
`python WM_climate_indices.py GLDASv2.0 SPI GLDASv2.0 [23,48,3,15] gamma monthly 6 2000 2018 1981 2010 ./results/`

4. PET from GLDAS
`python WM_climate_indices.py GLDASv2.0 PET GLDASv2.0 [23,48,3,15] gamma monthly 6 2000 2018 1981 2010 ./results/`

5. SPEI from GLDAS
`python WM_climate_indices.py GLDASv2.0 SPEI GLDASv2.0 [23,48,3,15] gamma monthly 6 2000 2018 1981 2010 ./results/`

## <a name = "req">Requirements</a>
Tested under Python 3.7

Package requirements:
* xarray - 0.12.1
* pandas - 0.25.0
* numpy - 1.16.4
* climate-indices - 1.0.6 Install this package from source (not pip)

## <a name = "files">Files in this directory</a>

* WM_climate_indices.py: Executable (climate indices calculation)
* DoughtIndexViz.py: Executable for the sample visualization
* spi_movie.gif: Sample visualization files

## <a name = "contact">Contact</a>

Please report issues to <khider@usc.edu>

## <a name ="license"> License </a>

The project is licensed under the Apache v2.0 License. Please refer to the file call license.
