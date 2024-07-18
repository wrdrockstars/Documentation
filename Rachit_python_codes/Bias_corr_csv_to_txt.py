# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 12:29:13 2024

@author: rachit
"""

# =============================================================================
# CSV to text 
# =============================================================================

import pandas

csvFile = pandas.read_csv(r"\\172.16.20.21\wrd_p_igbp\IGBP_PROJECT_2023-24\IGBP_Project_files\Data\Weather_data\Bias_corrected\p_ra_LS.csv")


group1 = csvFile[(csvFile['latitude'] >= 21.375) & (csvFile['latitude'] <= 32.125) & (csvFile['longitude'] >= 72.125) & (csvFile['longitude'] <= 89.875)]
group = group1[['latitude' , 'longitude']].drop_duplicates()  # To store lat lon pair in a separate dataframe and keep only unique combinations

# print(group.iloc[0]['latitude']) # To access the value of a particular entry in a particular column of a DataFrame
print(len(group))


for i in range(len(group)):
    filtered_csvFile = csvFile[(csvFile['latitude'] == group.iloc[i]['latitude']) & (csvFile['longitude'] == group.iloc[i]['longitude'])]
    pcp = filtered_csvFile.transpose().iloc[2:]  #transpose dataframe row to column and select rows from third to all using slicing
    pcp.to_csv(rf"\\172.16.20.21/wrd_p_igbp/IGBP_PROJECT_2023-24/IGBP_Project_files/Data/Weather_data/Bias_corrected/pcp{group.iloc[i]['latitude']}_{group.iloc[i]['longitude']}.txt".format(i), header=['19980101'], index=None, mode='w')
print('Converted to SWAT text format')


# =============================================================================
# Making points from data file
# =============================================================================
import geopandas as gpd
import pandas as pd

geometry = gpd.points_from_xy(group.longitude, group.latitude)
gdfele = gpd.GeoDataFrame(group, geometry=geometry, crs='EPSG:4326')

# save the GeoDataFrame as a shapefile
gdfele.to_file(r'\\172.16.20.21/wrd_p_igbp/IGBP_PROJECT_2023-24/IGBP_Project_files/Data/Weather_data/Bias_corrected/pcp_points.shp', driver='ESRI Shapefile')

# =============================================================================
# Now Go to ArcMap--> Extract values to point using DEM and saved elevattemp.shp
# =============================================================================
eledf = gpd.read_file(r"\\172.16.20.21/wrd_p_igbp/IGBP_PROJECT_2023-24/IGBP_Project_files/Data/Weather_data/Bias_corrected/elevation.shp")
eledf
data = []
for i in range(len(group)):
    data.append({
    'ID': i+1,
    'NAME': f"pcp{group.iloc[i]['latitude']}_{group.iloc[i]['longitude']}",
    'LAT': group.iloc[i]['latitude'],
    'LONG': group.iloc[i]['longitude'],
    'ELEVATION': eledf['RASTERVALU'][i]
    })
datacsvFile = pd.DataFrame(data)
# (To match with APHRODITE station files) sort first by longitude and then by latitude, ensuring that longitude is sorted in ascending order and latitude is sorted within each longitude group
datacsvFile.sort_values(by=['LONG', 'LAT'], ascending=[True, True], inplace=True) 
datacsvFile.to_csv(r"\\172.16.20.21/wrd_p_igbp/IGBP_PROJECT_2023-24/IGBP_Project_files/Data/Weather_data/Bias_corrected/pcp_stations.txt".format(i), header=True, index=None, mode='w')
