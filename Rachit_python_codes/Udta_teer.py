from matplotlib import pyplot as plt

from netCDF4 import Dataset
import numpy as np
import glob

from tqdm import tqdm


# =============================================================================
# Part 1 - Merging all .nc files into a single .nc file by adding time dimension
# =============================================================================


file_list = sorted(glob.glob(r"\\172.16.20.21\wrd_p_igbp\SCA_2000_2022\TEST\Data\*.nc"))

first_data = Dataset(file_list[0],'r')
lat = first_data.variables['latitude'][:]
lon = first_data.variables['longitude'][:]
lat_dim = len(lat)
lon_dim = len(lon)
var_name = 'fSCA'  # Change this to your variable name
# var_units = src.variables[var_name].units
first_data = None

# Create the output file
output_file = r'\\172.16.20.21\wrd_p_igbp\SCA_2000_2022\TEST\merged_data.nc'

merged_out_file = Dataset(output_file, 'w', format='NETCDF4')    

merged_out_file.createDimension('time', None)
merged_out_file.createDimension('lat', lat_dim)
merged_out_file.createDimension('lon', lon_dim)    

# Create variables
times = merged_out_file.createVariable('time', 'f4', ('time',))
latitudes = merged_out_file.createVariable('lat', 'f4', ('lat',))
longitudes = merged_out_file.createVariable('lon', 'f4', ('lon',))
fSCA = merged_out_file.createVariable(var_name, 'f4', ('time', 'lat', 'lon',),zlib=True, complevel=9)

# Assign attributes
latitudes.units = 'degrees north'
longitudes.units = 'degrees east'
# temperatures.units = var_units

# Assign latitude and longitude data
latitudes[:] = lat
longitudes[:] = lon

# Initialize the time variable
times.units = 'days since 2001-01-01'  # Change the reference date if necessary
times.calendar = 'gregorian'

# Iterate through the files and add data to the new file
time_index = 0
for file in tqdm(file_list):
    with Dataset(file, 'r') as src:
        fSCA_data = src.variables[var_name][:]
        fSCA[time_index, :, :] = fSCA_data
        # Assuming each file corresponds to one day, you can increment the time index
        times[time_index] = time_index
        time_index += 1

merged_out_file = None

print(f"Data successfully merged into {output_file}")



# =============================================================================
# Part 2 - Opening merged .nc file and extracting daily snow cover area
# =============================================================================

all_data = Dataset(output_file)

import geopandas as gpd
from shapely.geometry import Point
from shapely.geometry import mapping

# Load the shapefile
shapefile = r"\\172.16.20.21\\wrd_p_igbp\\IGBP_PROJECT_2023-24\\SnowCover_Himalayas\\Himalayas\\Himalayas.shp"
gdf = gpd.read_file(shapefile)

# Assuming your netCDF data has latitude and longitude variables
lat = all_data.variables['lat'][:]
lon = all_data.variables['lon'][:]
var_name = 'fSCA'  # Change to your variable name
data = all_data.variables[var_name][:]

# Create a mask for the region within the shapefile
def point_in_shapefile(lat, lon, gdf):
    point = Point(lon, lat)
    return any(g.contains(point) for g in gdf.geometry)

# =============================================================================
# Experiment Starts
# =============================================================================

for geom2 in gdf.geometry:
    continue

geom2_points = mapping(geom2)

# =============================================================================
# geom2_points = geom2_points['coordinates']
# =============================================================================

# Create an empty mask with the same shape as the data
mask = np.zeros((len(lat), len(lon)), dtype=bool)


for i in tqdm(range(len(lat))):
    for j in range(len(lon)):
        point = Point(lat[i],lon[j])        
        mask[i, j] = geom2.contains(point)

# =============================================================================
# Experiment Ends
# =============================================================================

# Apply the mask
for i in tqdm(range(len(lat))):
    for j in range(len(lon)):
        mask[i, j] = point_in_shapefile(lat[i], lon[j], gdf)     

# Use the mask to clip the data
clipped_data = np.where(mask, data, np.nan)

snow_cover_clipped_data = np.where((clipped_data > 50) & (clipped_data <=100), clipped_data, np.nan)

# =============================================================================
# =============================================================================
# # # Save the clipped data to a new netCDF file
# # output_file = 'clipped_data.nc'
# # with nc.Dataset(output_file, 'w', format='NETCDF4') as dst:
# #     # Create dimensions
# #     dst.createDimension('lat', len(lat))
# #     dst.createDimension('lon', len(lon))
# #     
# #     # Create variables
# #     latitudes = dst.createVariable('lat', 'f4', ('lat',))
# #     longitudes = dst.createVariable('lon', 'f4', ('lon',))
# #     clipped_var = dst.createVariable(var_name, 'f4', ('lat', 'lon',), fill_value=np.nan)
# #     
# #     # Assign data to variables
# #     latitudes[:] = lat
# #     longitudes[:] = lon
# #     clipped_var[:, :] = clipped_data
# # 
# #     # Copy attributes
# #     for name in ds.variables[var_name].ncattrs():
# #         clipped_var.setncattr(name, ds.variables[var_name].getncattr(name))
# =============================================================================
# =============================================================================

# =============================================================================
# print(f"Clipped data saved to {output_file}")
# =============================================================================
