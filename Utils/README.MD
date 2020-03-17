# Utilities

This folder contains executables to cut the global data for a certain region and concatenate the time axis.

One executable by type of data. Available data type:
*  GLDAS2

**Table of contents**

* [User Guide](#quickstart)
* [Requirements](#req)
* [Files in this repository](#files)

## <a name = "quickstart">User Guide</a>
Command line implementations:

`python DataTypeprep.py bounding_box path region output_name`

where:
* DataTypeprep.py is the python script for a particular dataset type
* bounding_box: List of [min_lon,max_lon,min_lat,max_lat]
* path: path to the file(s).
* Region: a string allowing to append the comment attribute to signal regional modification. Use `_` instead of space.
* dataset_name: Name of dataset file or parent folder
* output_name: The name of the concatenated netcdf file.

Example:
`python GLDASprep.py [10,60,-15,30] ./GLDASv2.1 Eastern_Equatorial_Africa trial.nc `

## <a name = "req">Requirements</a>
Tested under Python 3.7
Package requirements:
* xarray - 0.14.1
* dask - 2.2.0
* netcdf4 - 1.4.2
* bottleneck - 1.2.1
* numpy - 1.17.4

## <a name = "files">Files and in this directory</a>
* GLDASprep.py: Preparation script for GLDAS dataset.
* FLDASprep.py: Preparation script for FLDAS dataset.

**Note**: For FLDAS and GLDAS, assumes 1 file/month and the files are organized in yearly subfolder.
