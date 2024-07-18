# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 18:54:30 2024

@author: rachit
"""

# Pcp Long Term Annual Average
import glob
from osgeo import gdal
import numpy as np

from matplotlib import pyplot as plt

path = r"E:\\Rachit\\FLDAS_forcodetest\\wgetdata\\Pcp"

counter = 0

yearly_sum_of_yearly_sum = []

for year in range(2001,2021):    

    
    daily_file_list = glob.glob(path+"\*"+"A"+str(year)+"*.tif")


    if daily_file_list == []:            
        continue

    yearly_files = []

    for file in daily_file_list:                            
        data = gdal.Open(file)
        arr = data.ReadAsArray()
        arr[arr == -9999] = np.nan  # Replace with the actual NoData value                
        yearly_files.append(arr)
        # print(np.min(arr))
        # print(np.max(arr))
    
# =============================================================================
#         monthly_sum_array = np.array(monthly_sum)
#         monthly_mean = np.mean(monthly_sum_array,axis=0)
#         
#         plt.imshow(monthly_mean)
#         plt.show()
#         counter=counter+1
#         print(counter)
#         print("Month",month)
#         print("Year",year)
# =============================================================================


    yearly_sum_array = np.array(yearly_files)
    yearly_sum = np.sum(yearly_sum_array,axis=0)
    yearly_sum_of_yearly_sum.append(yearly_sum)
    yearly_sum = None

    if yearly_sum_of_yearly_sum == []:
        continue

yearly_sum_of_yearly_sum_array = np.array(yearly_sum_of_yearly_sum)
yearly_mean_of_yearly_sum = np.mean(yearly_sum_of_yearly_sum_array,axis=0)*86400  
# yearly_mean_of_yearly_means = np.where(yearly_mean_of_yearly_means<0,
#                                         0,
#                                         yearly_mean_of_yearly_means)

print("Max - ",np.max(yearly_mean_of_yearly_sum))
print("Min - ",np.min(yearly_mean_of_yearly_sum))
print("maxlocation - ",np.where(yearly_mean_of_yearly_sum == np.max(yearly_mean_of_yearly_sum)))
plt.imshow(yearly_mean_of_yearly_sum)
plt.colorbar()
plt.show()

########## Saving as tiff file
dst_filename=(f"E:\Rachit\FLDAS_forcodetest\output\pcp\pcp_yearly_mean.tif")
driver = gdal.GetDriverByName("GTiff")
# Specify the dimensions of the image
rows, cols = yearly_mean_of_yearly_sum.shape
# Create the output raster dataset
dst_dataset = driver.Create(dst_filename, cols, rows, bands=1, eType = gdal.GDT_Float64)

# Set required metadata and georeferencing information
from osgeo import osr
crs = osr.SpatialReference()
crs.ImportFromEPSG(4326)   #  EPSG code for WGS84
crs.ExportToWkt()
transform =data.GetGeoTransform()  #  Identity transform

dst_dataset.SetProjection(crs.ExportToWkt())
dst_dataset.SetGeoTransform(transform)

# Write the array data to the raster band
band = dst_dataset.GetRasterBand(1)
band.WriteArray(yearly_mean_of_yearly_sum)
#close the dataset to save changes
dst_dataset.FlushCache()
dst_dataset=None
# Clipping pcp annual avg
def ClipAll(infol,shpin,outfol):
    
    """Parameters are infol = path to input folder, shpin = path to input shape file, outfol = path to folder for clipped rasters

    Usage Example :-
    
    ClipAll(infol = r"D:\\Landsat\\Study Area", shpin = r"D:\\Study Area\\Study.shp", outfol = r"D:\\Study Area\\Clipped Output")
    
    use help(ClipAll) for clearer details
    
    """
    
    from osgeo import gdal
    import glob
    import datetime
    from tqdm import tqdm
    
    start_time = datetime.datetime.now()
    infiles = glob.glob(infol + "\*_yearly_mean.tif") 
    
    for i in tqdm(range(len(infiles))):
        split = infiles[i].split("\\")
        outfilename = split[-1][:-4]
        dsclip=gdal.Warp(destNameOrDestDS = outfol + '\\' + outfilename + '_clip' + '.tif',
                         srcDSOrSrcDSTab = infiles[i],
                         cutlineDSName=shpin,
                         cropToCutline=True,
                         dstNodata= np.nan)  #dstNodata to define no data value
        dsclip = None
        
    end_time = datetime.datetime.now()
    
    time_taken = end_time - start_time
    
    output = str(len(infiles))+" Files Successfully Clipped !!" + "\nTime Taken was " + str(time_taken)
    
    return(output)
names=['Asan','Aglar','Gangotri','Hanval','Jhelum','Parbati','Suketi','Tawi','Beas_Manali'] #
for name in names:
    result = ClipAll(r'E:\\Rachit\\FLDAS_forcodetest\\output\\pcp',r"E:\\Rachit\\NWH\\shapefiles_watershed\\"+name+".shp",r"E:\Rachit\FLDAS_forcodetest\output\pcp\\"+name)
    print(result)