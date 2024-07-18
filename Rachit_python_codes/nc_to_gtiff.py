# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 11:59:00 2024

@author: rachit
"""

# =============================================================================
# NC to Tiff
# =============================================================================
import glob
from netCDF4 import Dataset
from osgeo import gdal
from osgeo import osr

ds = gdal.Open(r"E:\Rachit\Extras\Practical_Mtech\3B-DAY.MS.MRG.3IMERG.20230601-S000000-E235959.V07B.tif")
dt = ds.ReadAsArray()

# Replace path with your actual folder path
folder_path = glob.glob(r"E:\\Rachit\\Extras\\Practical_Mtech\\Daily_IMERGv7nc\\*.nc4")

# Loop through each file in the folder
for file_path in folder_path:
    variable_name = 'precipitation'
    dataset = Dataset(file_path)
    data = dataset.variables[variable_name][0,:,:] # array slicing 1st argument is time (here it is only one day file but depends on netcdf structure) 
    
    dst_filename= rf"E:\Rachit\Extras\Practical_Mtech\Daily_IMERGv7tiff\New_folder\{file_path[-37:-4]}.tif" # output location
    driver = gdal.GetDriverByName("GTiff")
    # Specify the dimensions of the image
    rows, cols = dt.shape
    # Create the output raster dataset
    dst_dataset = driver.Create(dst_filename,
                                cols,
                                rows,
                                bands=1,
                                eType = gdal.GDT_Float64)

    # Set required metadata and georeferencing information
    
    crs = osr.SpatialReference()
    crs.ImportFromEPSG(4326)   #  EPSG code for WGS84
    crs.ExportToWkt()
    transform = ds.GetGeoTransform()
    # transform =[74.067552, 0.1000000229141223274, 0, 35.436552, 0, -0.1000000229141223274]

    dst_dataset.SetProjection(crs.ExportToWkt())
    dst_dataset.SetGeoTransform(transform)

    # Write the array data to the raster band
    band = dst_dataset.GetRasterBand(1)
    band.WriteArray(data)
    #close the dataset to save changes
    dst_dataset.FlushCache()
    dst_dataset=None

#############################################################################################
