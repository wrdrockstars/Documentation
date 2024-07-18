# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 15:57:10 2024

@author: rachit
"""
import xarray as xr #  for better data slicing


data = xr.open_dataset(r"\\172.16.20.21\wrd_p_igbp\IGBP_PROJECT_2023-24\IGBP_Project_files\Data\Weather_data\APHRODITE\APHRO_V1801_R1\APHRO_MA_025deg_V1801R1.1998.nc")

subset = data.sel(lat=slice(21.375, 32.125), lon=slice(72.125, 89.875))
subset


import numpy as np
import xarray as xr
import pandas as pd
import glob
import xarray as xr

# =============================================================================
# Concatenating nc files
# =============================================================================
# List of filenames to be combined
file_paths = glob.glob(r'\\172.16.20.21\wrd_p_igbp\IGBP_PROJECT_2023-24\IGBP_Project_files\Data\Weather_data\APHRODITE\APHRO_V1801_R1\*.nc')             

# Open all files
datasets = [xr.open_dataset(file) for file in file_paths]

# Concatenate along the desired dimension (e.g., time, lat, lon)
combined_dataset = xr.concat(datasets, dim='time')  # Adjust 'time' with the appropriate dimension in your files

# Save the combined dataset to a new file
combined_dataset.to_netcdf(r'\\172.16.20.21\wrd_p_igbp\IGBP_PROJECT_2023-24\IGBP_Project_files\Data\Weather_data\APHRODITE\APHRO_V1801_R1\aph_rain_1998_2015.nc')


# =============================================================================
# Converting to SWAT text
# =============================================================================

# ds = xr.open_mfdataset(r"E:\Rachit\IGBP_Project\Data\Weather_data\ERA5\total_precip\CDSAPI\daily-precip*.nc", combine='by_coords', concat_dim="time", engine='netcdf4') # open multiple files as a single dataset
# ============================================================================
# opening ncfile with data slicing
ncfile = xr.open_dataset(r"\\172.16.20.21\wrd_p_igbp\IGBP_PROJECT_2023-24\IGBP_Project_files\Data\Weather_data\APHRODITE\APHRO_V1801_R1\aph_rain_1998_2015.nc").sel(lat=slice(21.375, 32.125), lon=slice(72.125, 89.875))
# check = ncfile.sel(lat=21.375, lon=89.375)['precip'].values

df = ncfile.to_dataframe().reset_index()           # convert netcdf to dataframe and to reset the indexes according to netcdf dimensions

# Replace very small non-zero values with zero
df['precip'] = df['precip'].apply(lambda x: 0 if abs(x) < 1e-10 else x)
# Replace nan values with zero
df['precip'] = df['precip'].apply(lambda x: 0 if pd.isna(x) else x) #
# non_null_rows = df.dropna(subset=['precip'])         # REMOVED NAN VALUES (values of coordinates within the extent of shapefile), but not used in this case to ensure same no. of files are created

group = df[['lat' , 'lon']].drop_duplicates() # To store lat lon pair in a separate dataframe and keep only unique combinations

# print(group.iloc[0]['lat']) # To access the value of a particular entry in a particular column of a DataFrame
# print(group)


for i in range(len(group)):
    filtered_df = df[(df['lat'] == group.iloc[i]['lat']) & (df['lon'] == group.iloc[i]['lon'])]
    rain = filtered_df["precip"]
    rain.to_csv(rf"\\172.16.20.21\wrd_p_igbp\IGBP_PROJECT_2023-24\IGBP_Project_files\Data\Weather_data\APHRODITE\\APHRO_V1801_R1\\rain{group.iloc[i]['lat']}_{group.iloc[i]['lon']}.txt".format(i), header=['19980101'], index=None, mode='w')
print('Converted to SWAT text format')

# =============================================================================
# Pcp station input file preparation for SWAT
# =============================================================================
import geopandas as gpd
import pandas as pd

geometry = gpd.points_from_xy(group.lon, group.lat)
gdfele = gpd.GeoDataFrame(group, geometry=geometry, crs='EPSG:4326')

# save the GeoDataFrame as a shapefile
gdfele.to_file(r'\\172.16.20.21\\wrd_p_igbp\\IGBP_PROJECT_2023-24\\IGBP_Project_files\\Data\\Weather_data\\APHRODITE\\APHRO_V1801_R1\\precip_points.shp', driver='ESRI Shapefile')

# =============================================================================
# Now Go to ArcMap--> Extract values to point using DEM and saved elevattemp.shp
# =============================================================================

# shp to dataframe
eledf = gpd.read_file(r"\\172.16.20.21\\wrd_p_igbp\\IGBP_PROJECT_2023-24\\IGBP_Project_files\\Data\\Weather_data\\APHRODITE\\APHRO_V1801_R1\\precip_elevation.shp")
eledf
data = []
for i in range(len(group)):
    data.append({
    'ID': i+1,
    'NAME': f"rain{group.iloc[i]['lat']}_{group.iloc[i]['lon']}",
    'LAT': group.iloc[i]['lat'],
    'LONG': group.iloc[i]['lon'],
    'ELEVATION': eledf['RASTERVALU'][i]
    })
    datadf = pd.DataFrame(data)
    datadf.to_csv(r"\\172.16.20.21\wrd_p_igbp\IGBP_PROJECT_2023-24\IGBP_Project_files\Data\Weather_data\APHRODITE\APHRO_V1801_R1\SWAT_Text_files_ASCII\\precip_stations.txt".format(i), header=True, index=None, mode='w')
