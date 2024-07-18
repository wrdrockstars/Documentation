# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 16:15:55 2024

@author: rachit
"""
# from mpl_toolkits.basemap import Basemap


import xarray as xr
import matplotlib.pyplot as plt
import rioxarray
import geopandas
import pandas as pd
from shapely.geometry import mapping
import glob
import netCDF4 as nc
from datetime import datetime, timedelta
from tqdm import tqdm

# import warnings
# warnings.filterwarnings('ignore')

    
# =============================================================================
# # make a python list of dates starting from 1st Jan 2020 to 31st Dec 2020
# =============================================================================

# Define the start and end dates
start_date = datetime(2001, 1, 1)
end_date = datetime(2001, 12, 31)

# Initialize an empty list to store the dates
date_list = []

# Use a while loop to generate the dates
current_date = start_date
while current_date <= end_date:
    date_list.append(current_date)
    current_date += timedelta(days=1)

# # Print the list of dates
# date_list
# ============================================================================= 
path = glob.glob(rf"\\172.16.20.21\wrd_p_igbp\SCA_2000_2022\2002\*.NC")

# Adding time as 3rd dimension and saving them 
for idx, i in tqdm(enumerate(path)):     #enumerate(path) returns pairs of index and value from the path list
    data = xr.open_dataset(i)
    time_dim = [date_list[idx]]   # This var contains 
    ds_result = data.expand_dims(time=time_dim) # assign new dimension i.e. time to netcdf file
    ds_result.to_netcdf(rf"\\172.16.20.21\wrd_p_igbp\SCA_2000_2022\Output\2002\AWT_MODIS_FSC_{date_list[idx].strftime('%Y%m%d')}.NC")

data = xr.open_mfdataset(path, concat_dim='time', combine='nested')

# Clipping
mydata = data.rio.set_spatial_dims(x_dim="longitude", y_dim="latitude", inplace=True)
mydata = data.rio.write_crs("epsg:4326", inplace=True)

Himalayas = geopandas.read_file(r"\\172.16.20.21\\wrd_p_igbp\\IGBP_PROJECT_2023-24\\SnowCover_Himalayas\\Himalayas\\Himalayas.shp", crs="epsg:4326")

# Clip the data using the Himalayas geometry
geometries = [mapping(geom) for geom in Himalayas.geometry]
clipped = data.rio.clip(geometries, Himalayas.crs, all_touched=False)

# Plot the clipped data
snow = clipped['fSCA'].where(clipped['fSCA'] > 50)
# snow.plot()


Areas = []
for i in range(0,366):

    area = snow[i].count().values
    
    areasqkm = (int(area)*250000)/1000000
    
    Areas.append(areasqkm)
# #####################################

dataframe = pd.DataFrame(Areas) 
dataframe.to_csv(r'\\172.16.20.21\wrd_p_igbp\SCA_2000_2022\csv\2001.csv')

























for j in range(2001, 2020):
    
    # =============================================================================
    # # make a python list of dates starting from 1st Jan 2020 to 31st Dec 2020
    # =============================================================================
    
    # Define the start and end dates
    start_date = datetime(j, 1, 1)
    end_date = datetime(j, 12, 31)

    # Initialize an empty list to store the dates
    date_list = []

    # Use a while loop to generate the dates
    current_date = start_date
    while current_date <= end_date:
        date_list.append(current_date)
        current_date += timedelta(days=1)

    # # Print the list of dates
    # date_list
    # ============================================================================= 
    path = glob.glob(rf"\\172.16.20.21\wrd_p_igbp\SCA_2000_2022\{j}\*.NC")
    
    # Adding time as 3rd dimension and saving them 
    for idx, i in enumerate(path):     #enumerate(path) returns pairs of index and value from the path list
        data = xr.open_dataset(i)
        time_dim = [date_list[idx]]   # This var contains 
        ds_result = data.expand_dims(time=time_dim) # assign new dimension i.e. time to netcdf file
        ds_result.to_netcdf(rf"\\172.16.20.21\wrd_p_igbp\SCA_2000_2022\{j}\AWT_MODIS_FSC_{date_list[idx].strftime('%Y%m%d')}.NC")
    
    data = xr.open_mfdataset(glob.glob(rf"\\172.16.20.21\wrd_p_igbp\SCA_2000_2022\{j}\*.NC"))
       
    # Clipping
    mydata = data.rio.set_spatial_dims(x_dim="longitude", y_dim="latitude", inplace=True)
    mydata = data.rio.write_crs("epsg:4326", inplace=True)
    
    Himalayas = geopandas.read_file(r"\\172.16.20.21\\wrd_p_igbp\\IGBP_PROJECT_2023-24\\SnowCover_Himalayas\\Himalayas\\Himalayas.shp", crs="epsg:4326")
    
    # Clip the data using the Himalayas geometry
    geometries = [mapping(geom) for geom in Himalayas.geometry]
    clipped = data.rio.clip(geometries, Himalayas.crs, all_touched=False)
    
    # Plot the clipped data
    snow = clipped['fSCA'].where(clipped['fSCA'] > 50)
    # snow.plot()
    
    
    Areas = []
    for i in range(0,366):
    
        area = snow[i].count().values
        
        areasqkm = (int(area)*250000)/1000000
        
        Areas.append(areasqkm)
    # #####################################
    
    dataframe = pd.DataFrame(Areas) 
    dataframe.to_csv(r'\\172.16.20.21\wrd_p_igbp\SCA_2000_2022\csv\{j}.csv')


# plt.scatter(snow['time'], Areas)


