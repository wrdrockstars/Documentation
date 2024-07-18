# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 22:48:24 2024

@author: rachit
"""
import imdlib as imd

start_yr = 2001
end_yr = 2022
variable = 'rain' # other options are ('tmin'/ 'tmax'/ 'rain')
file_dir = (r'\\192.168.9.63\Share\NWH\UK_DEM') #Path to save the files
# imd.get_data(variable, start_yr, end_yr, fn_format='yearwise', file_dir=file_dir)
# %pwd # present working directory

## grd to xarray
# CHANGE HERE ACC TO PARAMETER REQUIRED
# data = imd.open_data('rain', 2001, 2022,'yearwise', file_dir) 

# import xarray
#  #Remove NaN values
# ds['rain'].mean('time').plot()

import numpy as np
import xarray
import matplotlib.pyplot as plt
# Assuming you have your files named like 'rain_2000.grd', 'rain_2001.grd', etc.
# Modify the path pattern according to your file naming convention.
file_pattern = '%d.grd'

# Initialize an array to store the sums
rainfall_sums = []

# Iterate over the years
for year in range(2000, 2023):
    # Load the dataset for the current year
    file_path = file_pattern % year
    dataset = imd.open_data('rain', year, year,'yearwise', file_dir)
    ds = dataset.get_xarray()
    ds = ds.where(ds['rain'] != -999.)
    # Assuming 'rain' is the variable name in your dataset
    # Sum the rainfall values for the entire year
    yearly_sum = ds['rain'].sum('time')
    
    # Append the sum to the array
    rainfall_sums.append(yearly_sum)
    
    # Close the dataset
    ds.close()

rainfall_array = np.array(rainfall_sums)

rainfall = np.mean(rainfall_sums,axis=0)
plt.imshow(rainfall, origin='lower')

type(rainfall)
# Now rainfall_sums will contain the yearly sums of rainfall
print(rainfall)

import imdlib as imd
import geopandas as gpd
from shapely.geometry import mapping
import rioxarray as rio
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import numpy as np
#Define lat/long 
pr = ds.rio.set_spatial_dims('lon', 'lat')
pr = pr.rio.write_crs("epsg:4326")

# Get the shapefile to be extracted
sf = gpd.read_file(r"\\192.168.9.63\Share\NWH\UK_DEM\Uttarakhand.shp")

# Extract/mask from the shapefile
clipped = pr.rio.clip(sf.geometry.apply(mapping), sf.crs, all_touched=False)  # clipped is an xarray dataset
# all_touched (bool, optional) – If True, all pixels touched by geometries will be burned in. If false, only pixels whose center is within the polygon or that are selected by Bresenham’s line algorithm will be burned in.

# Check for the values by plotting
fig, axs1 = plt.subplots(figsize=(5, 2.7))
clipped['rain'].mean('time').plot()
fig, axs2 = plt.subplots(figsize=(10, 2.7))
clipped['rain'].mean(axis=(1, 2)).plot.line()

# Save to file in NetCDF format
save_nc = (r'\\192.168.9.63\Share\NWH\UK_DEM\IMD_2000-2022_UK.nc')
clipped.to_netcdf(save_nc)   # only execute to save new file


