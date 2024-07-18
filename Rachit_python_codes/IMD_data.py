# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 10:14:49 2024

@author: rachit
"""

# =============================================================================
# Download IMD Rain, Temp Data & convert to SWAT Weather input Format
# =============================================================================

# pip install imdlib


# =============================================================================
# Download
# =============================================================================
import imdlib as imd

start_yr = 2023
end_yr = 2023
variable = 'tmin' # other options are ('tmin'/ 'tmax'/ 'rain')
file_dir = (r'E:\Rachit\IGBP_Project\Data\Weather_data\IMD_Ganga') #Path to save the files
# imd.get_data(variable, start_yr, end_yr, fn_format='yearwise', file_dir=file_dir)
# =============================================================================
%pwd # present working directory


# =============================================================================
# Changing imd grd format to xarray
# =============================================================================
# Opening grd data
data = imd.open_data('tmin', 2005, 2023,'yearwise', file_dir) # CHANGE HERE ACC TO PARAMETER REQUIRED

import xarray
ds = data.get_xarray()
print(ds)
# ds = ds.where(ds['rain'] != -999.) #Remove NaN values
# ds['rain'].mean('time').plot()  # plotting the mean rainfall over the years

# =============================================================================
# Clipping all grd files using shapefile and saving as nc file 
# =============================================================================
### (will keep all coordinates within the rect extent of shapefile)
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
sf = gpd.read_file(r"D:\Share\rachit\Stefy_Data\Watershed\watershed.shp")

# Extract/mask from the shapefile
clipped = pr.rio.clip(sf.geometry.apply(mapping), sf.crs, all_touched=True)         # clipped is an xarray dataset
# all_touched (bool, optional) – If True, all pixels touched by geometries will be burned in. If false, only pixels whose center is within the polygon or that are selected by Bresenham’s line algorithm will be burned in.

# Check for the values by plotting
# fig, axs1 = plt.subplots(figsize=(5, 2.7))
# clipped['rain'].mean('time').plot()
# fig, axs2 = plt.subplots(figsize=(10, 2.7))
# clipped['rain'].mean(axis=(1, 2)).plot.line()

# Save to file in NetCDF format
save_nc = (r'D:\Share\rachit\Stefy_Data\IMD_tmin_doyang.nc')
clipped.to_netcdf(save_nc)   # only execute to save new file

# =============================================================================
# Converting nc file to text files for SWAT weather input
# =============================================================================
# Rainfall
# import rioxarray as rio
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import numpy as np
import geopandas as gpd
from shapely.geometry import Point

ncfile = xr.open_dataset(r'D:\Share\rachit\Stefy_Data\IMD_tmax_doyang.nc')

df = ncfile.to_dataframe().reset_index()           # convert netcdf to dataframe and to reset the indexes according to netcdf dimensions

ncfile.close() # Closing the opened nc file

non_null_rows = df.dropna(subset=['rain'])         # REMOVED NAN VALUES (values of coordinates within the extent of shapefile)
non_null_rows

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
group = non_null_rows[['lat' , 'lon']].drop_duplicates() # To store lat lon pair in a separate dataframe and keep only unique combinations

print(group.iloc[0]['lat']) # To access the value of a particular entry in a particular column of a DataFrame
# print(group)

# SAMPLE CHECK
# filtered_df = non_null_rows[(non_null_rows['lat'] == group.iloc[0]['lat']) & (non_null_rows['lon'] == group.iloc[0]['lon'])]
# filtered_df

for i in range(len(group)):
    filtered_df = non_null_rows[(non_null_rows['lat'] == group.iloc[i]['lat']) & (non_null_rows['lon'] == group.iloc[i]['lon'])]
    rain = filtered_df["rain"]
    rain.to_csv(f"D://Share//rachit//Stefy_Data//rain//rain{group.iloc[i]['lat']}_{group.iloc[i]['lon']}.txt".format(i), header=['20050101'], index=None, mode='w')


#### Rain input file
##### --getting elevation values from DEM
### Converting dataframe to shapefile
import geopandas as gpd
import pandas as pd
# assuming the dataframe is named 'df' and has columns 'latitude' and 'longitude'
geometry = gpd.points_from_xy(group.lon, group.lat)
gdfele = gpd.GeoDataFrame(group, geometry=geometry, crs='EPSG:4326')

# save the GeoDataFrame as a shapefile
gdfele.to_file(r'D://Share//rachit//Stefy_Data\\rain_elevation.shp', driver='ESRI Shapefile')

# =============================================================================
# Now Go to ArcMap--> Extract values to points using DEM and saved elevattemp.shp
# =============================================================================


eledf1 = gpd.read_file(r"D://Share//rachit//Stefy_Data//rain_points.shp") 
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
    datadf.to_csv(f"D://Share//rachit//Stefy_Data//rain//rainstations.txt".format(i), header=True, index=None, mode='w')



## Combine tmin & tmax

ncfile1 = xr.open_dataset(r"D:\Share\rachit\Stefy_Data\IMD_tmax_doyang.nc")
ncfile2 = xr.open_dataset(r"D:\Share\rachit\Stefy_Data\IMD_tmin_doyang.nc")

df1 = ncfile1.to_dataframe().reset_index()
df2 = ncfile2.to_dataframe().reset_index()           # convert netcdf to dataframe and to reset the indexes according to netcdf dimensions

df3 = pd.concat([df1[['lat','lon','time','tmax']], df2['tmin']], axis=1)
df3
non_null_rows_t = df3.dropna(subset=['tmin','tmax'])       # REMOVED NAN VALUES (values of coordinates within the extent of shapefile)
non_null_rows_t 
# Rename the second column
# df = df.rename(columns={'A': 'C'})
group3 = non_null_rows_t[['lat' , 'lon']].drop_duplicates() # To store lat lon pair in a separate dataframe and keep only unique combinations
group3
print(group3.iloc[0]['lat']) # To access the value of a particular entry in a particular column of a DataFrame
print(group3)


# # SAMPLE CHECK
# filtered_df3 = non_null_rows_t[(non_null_rows_t['lat'] == group3.iloc[0]['lat']) & (non_null_rows_t['lon'] == group3.iloc[0]['lon'])]
# filtered_df3
for i in range(len(group3)):
    filtered_df3 = non_null_rows_t[(non_null_rows_t['lat'] == group3.iloc[i]['lat']) & (non_null_rows_t['lon'] == group3.iloc[i]['lon'])]
    tmax = filtered_df3['tmax']
    tmin = filtered_df3['tmin']
    result_df = pd.concat([tmax, tmin], axis=1)
    result_df.to_csv(f"D:\\Share\\rachit\\Stefy_Data/temp/temp{group3.iloc[i]['lat']}_{group3.iloc[i]['lon']}.txt".format(i), header=['20050101',''], index=None, mode='w')


    
#### Temp station input file
##### --getting elevation values from DEM
### Converting dataframe to shapefile
import geopandas as gpd
import pandas as pd

geometry3 = gpd.points_from_xy(group3.lon, group3.lat)
gdfele3 = gpd.GeoDataFrame(group3, geometry=geometry3, crs='EPSG:4326')

# save the GeoDataFrame as a shapefile
gdfele3.to_file(r'D:\\Share\\rachit\\Stefy_Data\\tempelevation.shp', driver='ESRI Shapefile')

# =============================================================================
# Now Go to ArcMap--> Extract values to point using DEM and saved elevattemp.shp
# =============================================================================

# shp to dataframe
eledf = gpd.read_file(r"D:\\Share\\rachit\\Stefy_Data\\temp_points.shp")
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
    datadf.to_csv(f"D:\\Share\\rachit\\Stefy_Data\\temp/tempstations.txt".format(i), header=True, index=None, mode='w')



### Extras
# group = non_null_rows.groupby(['lat', 'lon'])['rain']
# print(group)
# lat = non_null_rows['lat'].tolist()
# lon = non_null_rows['lon'].tolist()
# lat
# unique_lat_lon_df = non_null_rows[['time','rain']].drop_duplicates()
# unique_lat_lon_df


## Script to convert grd to csv
variable = 'rain' # change the variable name according to requirement. Other options are ('tmin'/ 'tmax')
if variable == 'rain':
    grid_size = 0.25 # grid spacing in deg
    y_count = 129 # no of grids in y direction
    x_count = 135 # no of grids in x direction
    x = 66.5 # starting longitude taken from control file (.ctl)
    y = 6.5 # starting latitude taken from control file (.ctl)
elif variable == 'tmax' or variable == 'tmin':
    grid_size = 1 # grid spacing in deg
    y_count = 31 # no of grids in y direction
    x_count = 31 # no of grids in x direction
    x = 67.5 # starting longitude taken from control file (.ctl)
    y = 7.5 # starting latitude taken from control file (.ctl)

#print(grid_size,x_count, y_count, x, y)
data
data.shape
np_array = data.data
#print(np_array[0,0,0])
#xr_objecct = data.get_xarray()
#type(xr_objecct)
#xr_objecct.mean('time').plot()
years_no = (end_yr - start_yr) + 1
#print(years_no)
day = 0
for yr in range(0,years_no):
    f = open("E:/Rachit/IGBP_Project/Data/Weather_data/IMD/"+str(start_yr+yr)+"_"+str(variable)+".csv",'w') # just change the path where you want to save csv file
    if ((start_yr+yr) % 4 == 0) and ((start_yr+yr) % 100 != 0):  # check for leap year
        days = 366
        count = yr + days
    elif ((start_yr+yr) % 4 == 0) and ((start_yr+yr) % 100 == 0) and ((start_yr+yr) % 400 == 0):
        days = 366
        count = yr + days
    else:
        days = 365
        count = yr + days

    day = day + days

    f.write("X,Y,")
    for d in range(0, days):
        f.write(str(d+1))
        f.write(",")
    f.write("\n")
    #print(np_array[364,0,0])
    for j in range(0, y_count):

        for i in range(0, x_count):

            f.write(str((i * grid_size) + x))
            f.write(",")
            f.write(str((j * grid_size) + y))
            f.write(",")
            time = 0
            for k in range(day-days, day):

                val = np_array[k,i,j]
                if val == 99.9000015258789 or val == -999:
                    f.write(str(-9999))
                    f.write(",")
                else:
                    f.write(str(val))
                    f.write(",")


            f.write("\n")
    print("File for " + str(start_yr + yr) + "_" + str(variable) + " is saved")
print("CSV conversion successful !")
## Getting elevations
import rasterio
from rasterio.transform import from_origin
import pandas as pd

# Load your DEM file
dem_path = "E:\Rachit\IGBP_Project\Data\DEM\MERIT_Ganga.tif"
dem_dataset = rasterio.open(dem_path)

# Read your DataFrame with lat and lon columns
# Assuming your DataFrame is named df
# You can replace this with your actual DataFrame
df = group

# Function to get elevation at a specific lat, lon
def get_elevation(lat, lon):
    row, col = dem_dataset.index(lon, lat)
    return dem_dataset.read(1)[row, col]

# Add a new column 'elevation' to your DataFrame
df['elevation'] = df.apply(lambda col: get_elevation(row['lat'], row['lon']), axis=1)

# Print the resulting DataFrame with elevation values
print(df)

dem_dataset
import rasterio

# Load your DEM file
dem_path = "E:\Rachit\IGBP_Project\Data\DEM\MERIT_Ganga_gcs.tif"
dem_dataset = rasterio.open(dem_path)

# Iterate over rows and columns to print latitude and longitude
for row in range(dem_dataset.height):
    for col in range(dem_dataset.width):
        lon, lat = dem_dataset.xy(row, col)
        print(f"Latitude: {lat}, Longitude: {lon}, Elevation: {dem_dataset.read(1)[row, col]}")