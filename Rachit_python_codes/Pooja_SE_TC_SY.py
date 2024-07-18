# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 22:48:40 2024

@author: rachit
"""

from osgeo import gdal
import numpy as np

# Open the raster files
src_se = gdal.Open(r"\\192.168.9.65\share_supriya\SEDIMENT_YIELD_DEPOSITION\SE\SE.tif")
src_tc = gdal.Open(r"\\192.168.9.65\share_supriya\SEDIMENT_YIELD_DEPOSITION\TC\TC.tif")

# Read raster data
se = src_se.GetRasterBand(1).ReadAsArray()
tc = src_tc.GetRasterBand(1).ReadAsArray()

# Compute SY based on conditions
sy = np.where(se == tc, tc, np.where(se < tc, se, tc))

# Compute Deposition based on conditions
deposition = np.where(se == tc, 0, np.where(se > tc, se - tc, 0))

#Checking
se_check = sy + deposition
print(se_check)

# Write SY and Deposition to output raster files
driver = gdal.GetDriverByName('GTiff')
profile = src_se.GetGeoTransform()  # Copy metadata from SE.tif

dst_sy = driver.Create(r'\\192.168.9.65\share_supriya\SEDIMENT_YIELD_DEPOSITION\SY.tif', src_se.RasterXSize, src_se.RasterYSize, 1, gdal.GDT_Float32)
dst_sy.SetGeoTransform(profile)
dst_sy.GetRasterBand(1).WriteArray(sy)
dst_sy.FlushCache()

dst_dep = driver.Create(r'\\192.168.9.65\share_supriya\SEDIMENT_YIELD_DEPOSITION\Deposition.tif', src_se.RasterXSize, src_se.RasterYSize, 1, gdal.GDT_Float32)
dst_dep.SetGeoTransform(profile)
dst_dep.GetRasterBand(1).WriteArray(deposition)
dst_dep.FlushCache()

dst_se_check = driver.Create(r'\\192.168.9.65\share_supriya\SEDIMENT_YIELD_DEPOSITION\SE_check.tif', src_se.RasterXSize, src_se.RasterYSize, 1, gdal.GDT_Float32)
dst_se_check.SetGeoTransform(profile)
dst_se_check.GetRasterBand(1).WriteArray(se_check)
dst_se_check.FlushCache()

# Close the raster files
src_se = None
src_tc = None
dst_sy = None
dst_dep = None
se_check = None
