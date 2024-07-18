# -*- coding: utf-8 -*-
"""
Created on Thu May  9 22:31:01 2024

@author: rachit
"""
import ee
ee.Initialize()
import geemap
Map = geemap.Map()
# Define the region of interest (ROI) as a feature collection (e.g., from a table)
roi = ee.FeatureCollection("projects/ee-rachitknl/assets/Rectangle")
# Define the latitude and longitude bounds for your area of interest
lat_min = 21.375
lat_max = 32.125
lon_min = 72.125
lon_max = 89.875

# Create a rectangle geometry for the area of interest
Map.addLayer(roi)
# Print the area of interest
print("Area of Interest:", roi)

# Define the ImageCollection with the specified ID and filter it based on time and region.
collection = ee.ImageCollection("ECMWF/ERA5/DAILY").select(['maximum_2m_air_temperature',
  # 'mean_2m_air_temperature',
  'minimum_2m_air_temperature'
  # 'total_precipitation',
  # 'u_component_of_wind_10m',
  # 'v_component_of_wind_10m'
  ]).filterBounds(roi).filterDate('1990-01-01', '1998-01-01')

# Create a function to extract values for all pixels in the ROI for each image
def extractValues(image): 
    date = image.date().format('yyyy-MM-dd')

    # Extract values for all pixels in the ROI
    values = image.sample(
    region= roi,
    scale= 27830,
    geometries= True
    )
    # Map over the values to create features
    def process_features(features, date):
        processed_features = features.map(lambda feature: (
            ee.Feature(
                None,
                feature.toDictionary().combine(
                    ee.Dictionary({
                        'date': date,
                        'latitude': ee.List(feature.geometry().coordinates().reverse()).get(0),
                        'longitude': ee.List(feature.geometry().coordinates().reverse()).get(1),
                        'maximum_2m_air_temperature': feature.getNumber('maximum_2m_air_temperature'),
                        # 'mean_2m_air_temperature': feature.getNumber('mean_2m_air_temperature'),
                        'minimum_2m_air_temperature': feature.getNumber('minimum_2m_air_temperature')
                        # 'total_precipitation': feature.getNumber('total_precipitation'),
                        # 'u_component_of_wind_10m': feature.getNumber('u_component_of_wind_10m'),
                        # 'v_component_of_wind_10m': feature.getNumber('v_component_of_wind_10m')
                    })
                )
            )
        ))
        return processed_features
    return values

# // Map the extractValues function over the ImageCollection
extractedData = collection.map(extractValues).flatten();

# // Filter out features with null values
extractedData = extractedData.filter(ee.Filter.notNull([
  'latitude',
  'longitude',
  'maximum_2m_air_temperature',
  # // 'mean_2m_air_temperature',
  'minimum_2m_air_temperature'
  # // 'total_precipitation',
  # // 'u_component_of_wind_10m',
  # // 'v_component_of_wind_10m'
]))

# // Print the first feature in the extractedData collection to inspect its properties
print('Extracted Data:', extractedData.first())


import os
os.chdir(r'E:\\Rachit\\IGBP_Project\\Data\\Weather_data\\ERA5\\Temp_2m')
geemap.ee_to_csv(extractedData, filename='era5_data_1990_1997.csv')

# Print the properties of the first feature in the extractedData collection
first_feature = extractedData.first()
properties = first_feature.propertyNames()
print('Properties of the first feature:', properties.getInfo())
