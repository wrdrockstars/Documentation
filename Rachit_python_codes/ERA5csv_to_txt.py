# -*- coding: utf-8 -*-
"""
Created on Thu May  9 15:25:18 2024

@author: rachit
"""

import os
import pandas as pd

def combine_csv_files(folder_path, output_file):
    # Get a list of all CSV files in the folder
    csv_files = [file for file in os.listdir(folder_path) if file.endswith('.csv')]
    
    if not csv_files:
        print("No CSV files found in the folder.")
        return
    
    # Read the first CSV file to get the header
    first_csv = pd.read_csv(os.path.join(folder_path, csv_files[0]))
    header = list(first_csv.columns)
    
    # Check if all CSV files have the same header
    for file in csv_files[1:]:
        df = pd.read_csv(os.path.join(folder_path, file))
        if list(df.columns) != header:
            print(f"CSV file {file} has a different header.")
            return
    
    # Combine all CSV files into a single DataFrame
    combined_df = pd.concat([pd.read_csv(os.path.join(folder_path, file)) for file in csv_files], ignore_index=True)
    
    # Write the combined DataFrame to a new CSV file
    combined_df.to_csv(output_file, index=False)
    print(f"Combined CSV files saved to {output_file}.")

# Usage example
folder_path = r"\\172.16.20.21\wrd_p_igbp\IGBP_PROJECT_2023-24\IGBP_Project_files\Data\Weather_data\ERA5\csv_files"
output_file = r"\\172.16.20.21\wrd_p_igbp\IGBP_PROJECT_2023-24\IGBP_Project_files\Data\Weather_data\ERA5\csv_files\era_1998_2019.csv"
combine_csv_files(folder_path, output_file)


# =============================================================================
# CSV to text 
# =============================================================================

import pandas

csvFile = pandas.read_csv(r"\\172.16.20.21\\wrd_p_igbp\\IGBP_PROJECT_2023-24\\IGBP_Project_files\\Data\\Weather_data\\ERA5\\csv_files\\era_1998_2019.csv")


group1 = csvFile[(csvFile['latitude'] >= 21.375) & (csvFile['latitude'] <= 32.125) & (csvFile['longitude'] >= 72.125) & (csvFile['longitude'] <= 89.875)]
group = group1[['latitude' , 'longitude']].drop_duplicates()  # To store lat lon pair in a separate dataframe and keep only unique combinations

# print(group.iloc[0]['latitude']) # To access the value of a particular entry in a particular column of a DataFrame
# print(group)


for i in range(len(group)):
    filtered_csvFile = csvFile[(csvFile['latitude'] == group.iloc[i]['latitude']) & (csvFile['longitude'] == group.iloc[i]['longitude'])]
    temp = filtered_csvFile[["maximum_2m_air_temperature","minimum_2m_air_temperature"]] - 273.15
    temp.to_csv(f"E:/Rachit/IGBP_Project/Data/Weather_data/ERA5/Temp_2m/CDSAPI/temp_max_min_1998_2019/celsius/temp{group.iloc[i]['latitude']}_{group.iloc[i]['longitude']}.txt".format(i), header=['19980101',''], index=None, mode='w')
print('Converted to SWAT text format')


# =============================================================================
# Making points from data file
# =============================================================================
import geopandas as gpd
import pandas as pd

geometry = gpd.points_from_xy(group.longitude, group.latitude)
gdfele = gpd.GeoDataFrame(group, geometry=geometry, crs='EPSG:4326')

# save the GeoDataFrame as a shapefile
gdfele.to_file(r'E:\\Rachit\\IGBP_Project\\Data\\Weather_data\\ERA5\\Temp_2m\\CDSAPI\\temp_points.shp', driver='ESRI Shapefile')

# =============================================================================
# Now Go to ArcMap--> Extract values to point using DEM and saved elevattemp.shp
# =============================================================================
eledf = gpd.read_file(r"E:\Rachit\IGBP_Project\Data\Weather_data\ERA5\Temp_2m\CDSAPI\temp_max_min_1998_2019\temp_elevation1.shp")
eledf
data = []
for i in range(len(group)):
    data.append({
    'ID': i+1,
    'NAME': f"temp{group.iloc[i]['latitude']}_{group.iloc[i]['longitude']}",
    'LAT': group.iloc[i]['latitude'],
    'LONG': group.iloc[i]['longitude'],
    'ELEVATION': eledf['RASTERVALU'][i]
    })
datacsvFile = pd.DataFrame(data)
# (To match with APHRODITE station files) sort first by longitude and then by latitude, ensuring that longitude is sorted in ascending order and latitude is sorted within each longitude group
datacsvFile.sort_values(by=['LONG', 'LAT'], ascending=[True, True], inplace=True) 
datacsvFile.to_csv(f"E:\\Rachit\\IGBP_Project\\Data\\Weather_data\\ERA5\\Temp_2m\\CDSAPI\\temp_stations.txt".format(i), header=True, index=None, mode='w')
