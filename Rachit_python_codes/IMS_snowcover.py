# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 11:53:25 2014

@author: rachit
"""

# =============================================================================
# Snow Cover Time Series Analysis
# =============================================================================
import pandas as pd
from osgeo import gdal
import numpy as np
import matplotlib.pyplot as plt
import tempfile
import glob
import os
import zipfile
import gzip  # to unzip .gz files
import shutil

files = glob.glob(r"\\172.16.20.21\wrd_p_igbp\LSP_Project\SnowCover\GIS_1km\2014\*.gz")

# =============================================================================
# Unzipping the files
# =============================================================================

# Create temporary directory for unzipped files
os.chdir(r"\\172.16.20.21\wrd_p_igbp\LSP_Project\SnowCover\GIS_1km")
temp_dir = tempfile.mkdtemp(dir=r'\\172.16.20.21\wrd_p_igbp\LSP_Project\SnowCover\GIS_1km')

for file in files:
    # print(file)
    
    with gzip.open(file, 'rb') as f_in:
        with open(temp_dir + f"\{file[-30:]}.tif", 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        # Get the folder name
      # foldername = [
      #     info.filename for info in zf.infolist()
      #     if info.is_dir()][0]
      # Extract all the data
      # zf.extractall(temp_dir)

# =============================================================================
# Clipping all files
# =============================================================================
def ClipAll(infol,shpin,outfol):
        
    from osgeo import gdal
    import glob
    import datetime
    from tqdm import tqdm
    
    start_time = datetime.datetime.now()
    infiles = glob.glob(infol + "\*.tif") 
    
    for i in tqdm(range(len(infiles))):
        split = infiles[i].split("\\")
        outfilename = split[-1][:-4]
        dsclip=gdal.Warp(destNameOrDestDS = outfol + '\\' + outfilename + '_clip' + '.tif',
                         srcDSOrSrcDSTab = infiles[i],
                         cutlineDSName=shpin,
                         cropToCutline=True,
                         dstNodata = 0)
        dsclip = None
        
    end_time = datetime.datetime.now()
    
    time_taken = end_time - start_time
    
    output = str(len(infiles))+" Files Successfully Clipped !!" + "\nTime Taken was " + str(time_taken)
    
    return(output)

result = ClipAll(temp_dir,r"\\172.16.20.21\wrd_p_igbp\LSP_Project\SnowCover\UK1.shp", r'\\172.16.20.21\wrd_p_igbp\LSP_Project\SnowCover\GIS_1km\UK_clip\2014')



shutil.rmtree(temp_dir) # to delete the temporary unzipped files

snowarea= []
tsnowarea = []
def snowcoverarea(path):
    tiff = gdal.Open(path)
#     lulc.GetProjection() 
    tiffilearr = tiff.ReadAsArray()
    allpixels = (tiffilearr == 1) + (tiffilearr == 2) + (tiffilearr == 3) + (tiffilearr == 4)
    snowpixels = tiffilearr == 4        
    
    snow = np.extract(snowpixels, tiffilearr)
    totalarea = np.extract(allpixels, tiffilearr)
    totalsnowarea = len(snow)*1*1 #since the spatial resolution of ims snow data is 1m hence multiplied by 1
    area = (len(snow)*1*1 ) / len(totalarea) 

    print(f'relative snow covered area in UK: {area} km\u00b2, \ntotal area covered with snow {(len(snow)*1*1 )} km\u00b2') 
          
    snowarea.append(area)
    tsnowarea.append(totalsnowarea)

folder = glob.glob(r'\\172.16.20.21\\wrd_p_igbp\\LSP_Project\\SnowCover\\GIS_1km\\UK_clip\\2014\\*.tif')


for i in folder:
    print(f'For {i}:')
    print(snowcoverarea(i))
    
# =============================================================================
# List to Dataframe to Excel
# =============================================================================
filename = [file_name[-36:-29] for file_name in folder]
# Convert list to DataFrame
df = pd.DataFrame({'Date': filename,'Relative Snow Covered Area (ratio)': snowarea, 'total area covered with snow (km²)': tsnowarea })

# Specify the Excel file path
excel_file_path = r'\\172.16.20.21\wrd_p_igbp\LSP_Project\SnowCover\GIS_1km\UK_clip\2014\snow_area.xlsx'

# Save DataFrame to Excel
df.to_excel(excel_file_path, index=False)

print(f"Data saved to {excel_file_path}")


# =============================================================================
# For Creating Binary TIFF files
# =============================================================================

# Define the output directories
output_directory = r"\\172.16.20.21\wrd_p_igbp\LSP_Project\SnowCover\GIS_1km\UK_clip\\2014"

# Create the output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

# Iterate through each TIFF file in the input directory
for i in folder:
    dataset = gdal.Open(i)
        
    # Read the data as an array
    tiffile = dataset.ReadAsArray()
    
    # Create the binary mask
    binary_mask = np.zeros_like(tiffile, dtype=np.uint32)
    binary_mask[(tiffile == 1) | (tiffile == 2) | (tiffile == 3)] = 0
    binary_mask[(tiffile == 4)] = 1
    binary_mask[(tiffile == 0)] = -9999 #nodata
    
    # Define the output file path
    output_path = os.path.join(output_directory, f"{os.path.splitext(i)[0]}_binary.tif")
    
    # Write the binary mask to a new TIFF file
    driver = gdal.GetDriverByName("GTiff")
    output_dataset = driver.Create(output_path, dataset.RasterXSize, dataset.RasterYSize, 1, gdal.GDT_Byte)
    output_dataset.SetGeoTransform(dataset.GetGeoTransform())
    output_dataset.SetProjection(dataset.GetProjection())
    output_dataset.GetRasterBand(1).WriteArray(binary_mask)
    output_dataset.FlushCache()
    output_dataset = None
    
    print(f"Binary mask saved to {output_path}")






