# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 17:29:41 2024

@author: rachit
"""
# =============================================================================
# Snowcover data for MHAM basins
# =============================================================================
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

result = ClipAll(r'\\172.16.20.21\wrd_p_igbp\LSP_Project\SnowCover\GIS_1km\UK_clip\2024\binary',r"E:/Rachit/Extras/MHAM/mham_boundary/alaknanda.shp", r'\\172.16.20.21\wrd_p_igbp\LSP_Project\SnowCover\GIS_1km\UK_clip\2024\binary\Alaknanda')



snowarea= []
tsnowarea = []
def snowcoverarea(path):
    tiff = gdal.Open(path)
#     lulc.GetProjection() 
    tiffilearr = tiff.ReadAsArray()
    allpixels = (tiffilearr == 1) + (tiffilearr == 0)
    snowpixels = tiffilearr == 1        
    
    snow = np.extract(snowpixels, tiffilearr)
    totalarea = np.extract(allpixels, tiffilearr)
    totalsnowarea = len(snow)*1*1 #since the spatial resolution of ims snow data is 1m hence multiplied by 1
    area = (len(snow)*1*1 ) / len(totalarea) 

    print(f'relative snow covered area in Alaknanda: {area} km\u00b2, \ntotal area covered with snow {(len(snow)*1*1 )} km\u00b2') 
          
    snowarea.append(area)
    tsnowarea.append(totalsnowarea)

folder = glob.glob(r'\\172.16.20.21\wrd_p_igbp\LSP_Project\SnowCover\GIS_1km\UK_clip\2024\binary\Alaknanda\*.tif')


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
excel_file_path = r'\\172.16.20.21\wrd_p_igbp\LSP_Project\SnowCover\GIS_1km\UK_clip\2024\binary\Alaknanda\snow_area_Alk_2024.xlsx'

# Save DataFrame to Excel
df.to_excel(excel_file_path, index=False)

print(f"Data saved to {excel_file_path}")