# Regridding and data transformation

## Invocaation

`python regrid_to_csv.py path_to_netcdf regrid_flag`

## Options

path_to_netcdf: Path to the netcdf file to transform
regrid_flag : If set to True, will zoom onto Ethiopia on a 12x12 spatial grid

## Returns

a netcdf file with the regridded values. note that if regrid_flag is set to false, will be the same as the original file
a csv file with the regridded values

