# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 19:27:39 2024

@author: rachit
"""

import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import numpy as np
import geopandas as gpd
from shapely.geometry import Point


file = pd.read_csv("D:\\Share\\rachit\\SWAT_data\\Rainfall_data\\rain_data_1993_2023.csv")

# df = ncfile.to_dataframe().reset_index()           # convert netcdf to dataframe and to reset the indexes according to netcdf dimensions

# non_null_rows = df.dropna(subset=['rain'])         # REMOVED NAN VALUES (values of coordinates within the extent of shapefile)
# non_null_rows

# -------------------------------------------------------------------------------------------
# # Convert DataFrame to a GeoDataFrame by creating a geometry column with Point geometries
# geometry = [Point(xy) for xy in zip(non_null_rows['lat'], non_null_rows['lon'])]
# gdf = gpd.GeoDataFrame(non_null_rows, geometry=geometry, crs="+proj=lcc +lat_1=12.472944444 +lat_2=35.172805555556 +lat_0=24 +lon_0=80 +x_0=4000000 +y_0=4000000 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
# gdf.crs
# gdf.set_crs("+proj=lcc +lat_1=12.472944444 +lat_2=35.172805555556 +lat_0=24 +lon_0=80 +x_0=4000000 +y_0=4000000 +ellps=WGS84 +datum=WGS84 +units=m")
# perform a spatial join between GeoDataFrame and the shapefile. This will assign each point in DataFrame to the corresponding polygon in the shapefile.
# sf = gpd.read_file(r'E:\\Rachit\\IGBP_Project\\Yamuna_test_shp\\Yamuna_test_shp.shp')
# gdf_joined = gpd.sjoin(gdf, sf, how='inner', predicate='within')
# gdf_joined
# sf
# -------------------------------------------------------------------------------------------
group = file[['lat' , 'lon']].drop_duplicates() # To store lat lon pair in a separate dataframe and keep only unique combinations

print(group.iloc[0]['lat']) # To access the value of a particular entry in a particular column of a DataFrame
# print(group)

# SAMPLE CHECK
# filtered_df = non_null_rows[(non_null_rows['lat'] == group.iloc[0]['lat']) & (non_null_rows['lon'] == group.iloc[0]['lon'])]
# filtered_df

for i in range(len(group)):
    filtered_df = file[(file['lat'] == group.iloc[i]['lat']) & (file['lon'] == group.iloc[i]['lon'])]
    rain = filtered_df["rain"]
    rain.to_csv(f"D://Share//rachit//SWAT_data//Rainfall_data//rain{group.iloc[i]['lat']}_{group.iloc[i]['lon']}.txt".format(i), header=['19930101'], index=None, mode='w')


#### Rain input file
##### --getting elevation values from DEM
### Converting dataframe to shapefile
import geopandas as gpd
import pandas as pd
# assuming the dataframe is named 'df' and has columns 'latitude' and 'longitude'
geometry = gpd.points_from_xy(group.lon, group.lat)
gdfele = gpd.GeoDataFrame(group, geometry=geometry, crs='EPSG:4326')

# save the GeoDataFrame as a shapefile
gdfele.to_file(r'D:\\Share\\rachit\\SWAT_data\\Rainfall_data\\rainelevation.shp', driver='ESRI Shapefile')

# =============================================================================
# Now Go to ArcMap--> Extract values to point using DEM and saved elevattemp.shp
# =============================================================================


eledf1 = gpd.read_file(r"D:\\Share\\rachit\\SWAT_data\\Rainfall_data\\rainelevation_values.shp") 
# eledf1
data = []
for i in range(len(group)):
    data.append({
    'ID': i+1,
    'NAME': f"rain{group.iloc[i]['lat']}_{group.iloc[i]['lon']}",
    'LAT': group.iloc[i]['lat'],
    'LONG': group.iloc[i]['lon'],
    'ELEVATION': eledf1['RASTERVALU'][i]
    })
    datadf = pd.DataFrame(data)
    datadf.to_csv(f"D:\\Share\\rachit\\SWAT_data\\Rainfall_data\\rainstations.txt".format(i), header=True, index=None, mode='w')



### Temperature
import rioxarray as rio
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import numpy as np

## Combine tmin & tmax

filetmax = pd.read_csv(r"D:\Share\rachit\SWAT_data\Temperature_data\tmax_1993_2023.csv")
filetmin = pd.read_csv(r"D:\Share\rachit\SWAT_data\Temperature_data\tmin_1993_2023.csv")

df3 = pd.concat([filetmax[['lat','lon','time','tmax']], filetmin['tmin']], axis=1)


# Rename the second column
# df = df.rename(columns={'A': 'C'})
group3 = df3[['lat' , 'lon']].drop_duplicates() # To store lat lon pair in a separate dataframe and keep only unique combinations

print(group3.iloc[0]['lat']) # To access the value of a particular entry in a particular column of a DataFrame



# # SAMPLE CHECK
# filtered_df3 = non_null_rows_t[(non_null_rows_t['lat'] == group3.iloc[0]['lat']) & (non_null_rows_t['lon'] == group3.iloc[0]['lon'])]
# filtered_df3
for i in range(len(group3)):
    filtered_df3 = df3[(df3['lat'] == group3.iloc[i]['lat']) & (df3['lon'] == group3.iloc[i]['lon'])]
    tmax = filtered_df3['tmax']
    tmin = filtered_df3['tmin']
    result_df = pd.concat([tmax, tmin], axis=1)
    result_df.to_csv(f"D:/Share/rachit/SWAT_data/Temperature_data/temp{group3.iloc[i]['lat']}_{group3.iloc[i]['lon']}.txt".format(i), header=['19930101',''], index=None, mode='w')

# =============================================================================
# Remove commas Source: ChatGPT
# =============================================================================
import os
import glob

# Specify the directory containing the text files
directory = r'D:\Share\rachit\SWAT_data\Temperature_data'  # Replace with the path to your directory

# Use glob to find all text files in the directory
file_paths = glob.glob(os.path.join(directory, '*.txt'))
a = file_paths
for file_path in file_paths:
    # Open the file for reading
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Remove the trailing comma from the first line if it exists
    if lines:
        lines[0] = lines[0].rstrip(',\n') + '\n'

    # Write the modified content back to the file
    with open(file_path, 'w') as file:
        file.writelines(lines)

print("Comma removed from the first line of all files.")
    
#### Temp station input file
##### --getting elevation values from DEM
### Converting dataframe to shapefile
import geopandas as gpd
import pandas as pd

geometry3 = gpd.points_from_xy(group3.lon, group3.lat)
gdfele3 = gpd.GeoDataFrame(group3, geometry=geometry3, crs='EPSG:4326')

# save the GeoDataFrame as a shapefile
gdfele3.to_file(r'D:\\Share\\rachit\\SWAT_data\\Temperature_data\\tempelevation.shp', driver='ESRI Shapefile')

# =============================================================================
# Now Go to ArcMap--> Extract values to point using DEM and saved elevattemp.shp
# =============================================================================

# shp to dataframe
eledf = gpd.read_file(r"D:\\Share\\rachit\\SWAT_data\\Temperature_data\\tempelevation_values.shp")
eledf
data = []
for i in range(len(group3)):
    data.append({
    'ID': i+1,
    'NAME': f"temp{group3.iloc[i]['lat']}_{group3.iloc[i]['lon']}",
    'LAT': group3.iloc[i]['lat'],
    'LONG': group3.iloc[i]['lon'],
    'ELEVATION': eledf['RASTERVALU'][i]
    })
    datadf = pd.DataFrame(data)
    datadf.to_csv(f"D://Share//rachit//SWAT_data//Temperature_data//tempstations.txt".format(i), header=True, index=None, mode='w')
