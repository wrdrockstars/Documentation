# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 19:46:05 2024

@author: rachit
"""

# =============================================================================
# # Snow Cover Time Series Analysis
# =============================================================================
## Clipping all files
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
    infiles = glob.glob(infol + "\\*.tif") 
    
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

ClipAll(r'E:\Rachit\Snow\2014',r'E:\\Rachit\\Snow\\Snow_boundary\\NWH_states_dissolve.shp',r'E:\Rachit\Snow\2014\clipped')

# =============================================================================
# # plotting
# =============================================================================
import pandas as pd
from osgeo import gdal
import numpy as np
import matplotlib.pyplot as plt
import glob
folder1 = glob.glob('E:\\Rachit\\Snow\\DecJan2015-24\\2023-24\\*.tif')
folder2 = glob.glob('E:\\Rachit\\Snow\\DecJan2015-24\\2022-23\\*.tif')
folder3 = glob.glob('E:\\Rachit\\Snow\\DecJan2015-24\\2021-22\\*.tif')
folder4 = glob.glob('E:\\Rachit\\Snow\\DecJan2015-24\\2020-21\\*.tif')
folder5 = glob.glob('E:\\Rachit\\Snow\\DecJan2015-24\\2019-20\\*.tif')
folder6 = glob.glob('E:\\Rachit\\Snow\\DecJan2015-24\\2018-19\\*.tif')
folder7 = glob.glob('E:\\Rachit\\Snow\\DecJan2015-24\\2017-18\\*.tif')
folder8 = glob.glob('E:\\Rachit\\Snow\\DecJan2015-24\\2016-17\\*.tif')
folder9 = glob.glob('E:\\Rachit\\Snow\\DecJan2015-24\\2015-16\\*.tif')
folder10 = glob.glob('E:\\Rachit\\Snow\\DecJan2015-24\\2014-15\\*.tif')
dates = ['Dec1', 'Dec2', 'Dec3', 'Dec4', 'Dec5',
 
 'Dec6',

 'Dec7',

 'Dec8',

 'Dec9',

 'Dec10',

 'Dec11',

 'Dec12',

 'Dec13',

 'Dec14',

 'Dec15',

 'Dec16',

 'Dec17',

 'Dec18',

 'Dec19',

 'Dec20',

 'Dec21',

 'Dec22',

 'Dec23',

 'Dec24',

 'Dec25',

 'Dec26',

 'Dec27',
 
 'Dec28',

 'Dec29',

 'Dec30',
 
 'Dec31', 'Jan1', 'Jan2', 'Jan3', 'Jan4', 'Jan5', 
 'Jan6',
 
 'Jan7',
 
 'Jan8',
 
 'Jan9',
 
 'Jan10',
 
 'Jan11',
 
 'Jan12',

 'Jan13',

 'Jan14',
 
 'Jan15',

 'Jan16',
 
 'Jan17',

 'Jan18',

 'Jan19',

 'Jan20',
 
 'Jan21',
 
 'Jan22',

 'Jan23',

 'Jan24',
 
 'Jan25',

 'Jan26',
 
 'Jan27',

 'Jan28',

 'Jan29',

 'Jan30',
 
 'Jan31']
# for i in range(1,32):
#     decday = f'Dec{i}:'
#     dates.append(decday)
#     janday = f'Jan{i}'
#     dates.append(janday)
# dates
snowarea1= []
def snowcoverarea1(path):
    tiff1 = gdal.Open(path)
#     lulc.GetProjection() 
    tiffilearr1 = tiff1.ReadAsArray()
    allpixels1 = (tiffilearr1 == 1) + (tiffilearr1 == 2) + (tiffilearr1 == 3) + (tiffilearr1 == 4)
    snowpixels1 = tiffilearr1 == 4      
    
    # 7 value is allocated to pixels in built area class in Sentinel-2 10m Land Use/Land Cover timeseries  
    
    snow1 = np.extract(snowpixels1, tiffilearr1)
    totalarea1 = np.extract(allpixels1, tiffilearr1)
    area1 = (len(snow1)*1*1 ) / len(totalarea1) #since the spatial resolution of ims snow data is 1m
    
    print(f'relative snow covered area in the NWH: {area1} km\u00b2') 
    snowarea1.append(area1)

snowarea2= []
def snowcoverarea2(path):
    tiff2 = gdal.Open(path)
#     lulc.GetProjection() 
    tiffilearr2 = tiff2.ReadAsArray()
    allpixels2 = (tiffilearr2 == 1) + (tiffilearr2 == 2) + (tiffilearr2 == 3) + (tiffilearr2 == 4)
    snowpixels2 = tiffilearr2 == 4      
    
    # 7 value is allocated to pixels in built area class in Sentinel-2 10m Land Use/Land Cover timeseries  
    
    snow2 = np.extract(snowpixels2, tiffilearr2)
    totalarea2= np.extract(allpixels2, tiffilearr2)
    area2 = (len(snow2)*1*1 ) / len(totalarea2) #since the spatial resolution of ims snow data is 1m
    
    print(f'relative snow covered area in the NWH: {area2} km\u00b2') 
    snowarea2.append(area2)
snowarea3= []
def snowcoverarea3(path):
    tiff3 = gdal.Open(path)
#     lulc.GetProjection() 
    tiffilearr3 = tiff3.ReadAsArray()
    allpixels3 = (tiffilearr3 == 1) + (tiffilearr3 == 2) + (tiffilearr3 == 3) + (tiffilearr3 == 4)
    snowpixels3 = tiffilearr3 == 4      
    
    # 7 value is allocated to pixels in built area class in Sentinel-2 10m Land Use/Land Cover timeseries  
    
    snow3 = np.extract(snowpixels3, tiffilearr3)
    totalarea3 = np.extract(allpixels3, tiffilearr3)
    area3 = (len(snow3)*1*1 ) / len(totalarea3) #since the spatial resolution of ims snow data is 1m
    
    print(f'relative snow covered area in the NWH: {area3}') 
    snowarea3.append(area3)
snowarea4= []
def snowcoverarea4(path):
    tiff4 = gdal.Open(path)
#     lulc.GetProjection() 
    tiffilearr4 = tiff4.ReadAsArray()
    allpixels4 = (tiffilearr4 == 1) + (tiffilearr4 == 2) + (tiffilearr4 == 3) + (tiffilearr4 == 4)
    snowpixels4 = tiffilearr4 == 4      
    
    # 7 value is allocated to pixels in built area class in Sentinel-2 10m Land Use/Land Cover timeseries  
    
    snow4 = np.extract(snowpixels4, tiffilearr4)
    totalarea4 = np.extract(allpixels4, tiffilearr4)
    area4 = (len(snow4)*1*1 ) / len(totalarea4) #since the spatial resolution of ims snow data is 1m
    
    print(f'relative snow covered area in the NWH: {area4}') 
    snowarea4.append(area4)
snowarea5= []
def snowcoverarea5(path):
    tiff5 = gdal.Open(path)
#     lulc.GetProjection() 
    tiffilearr5 = tiff5.ReadAsArray()
    allpixels5 = (tiffilearr5 == 1) + (tiffilearr5 == 2) + (tiffilearr5 == 3) + (tiffilearr5 == 4)
    snowpixels5 = tiffilearr5 == 4      
    
    # 7 value is allocated to pixels in built area class in Sentinel-2 10m Land Use/Land Cover timeseries  
    
    snow5 = np.extract(snowpixels5, tiffilearr5)
    totalarea5 = np.extract(allpixels5, tiffilearr5)
    area5 = (len(snow5)*1*1 ) / len(totalarea5) #since the spatial resolution of ims snow data is 1m
    
    print(f'relative snow covered area in the NWH: {area5} ') 
    snowarea5.append(area5)
snowarea6= []
def snowcoverarea6(path):
    tiff6 = gdal.Open(path)
#     lulc.GetProjection() 
    tiffilearr6 = tiff6.ReadAsArray()
    allpixels6 = (tiffilearr6 == 1) + (tiffilearr6 == 2) + (tiffilearr6 == 3) + (tiffilearr6 == 4)
    snowpixels6 = tiffilearr6 == 4      
    
    # 7 value is allocated to pixels in built area class in Sentinel-2 10m Land Use/Land Cover timeseries  
    
    snow6 = np.extract(snowpixels6, tiffilearr6)
    totalarea6 = np.extract(allpixels6, tiffilearr6)
    area6 = (len(snow6)*1*1 ) / len(totalarea6) #since the spatial resolution of ims snow data is 1m
    
    print(f'relative snow covered area in the NWH : {area6}') 
    snowarea6.append(area6)
snowarea7= []
def snowcoverarea7(path):
    tiff7 = gdal.Open(path)
#     lulc.GetProjection() 
    tiffilearr7 = tiff7.ReadAsArray()
    allpixels7 = (tiffilearr7 == 1) + (tiffilearr7 == 2) + (tiffilearr7 == 3) + (tiffilearr7 == 4)
    snowpixels7 = tiffilearr7 == 4      
    
    # 7 value is allocated to pixels in built area class in Sentinel-2 10m Land Use/Land Cover timeseries  
    
    snow7 = np.extract(snowpixels7, tiffilearr7)
    totalarea7 = np.extract(allpixels7, tiffilearr7)
    area7 = (len(snow7)*1*1 ) / len(totalarea7) #since the spatial resolution of ims snow data is 1m
    
    print(f'relative snow covered area in the NWH: {area7}') 
    snowarea7.append(area7)
snowarea8= []
def snowcoverarea8(path):
    tiff8 = gdal.Open(path)
#     lulc.GetProjection() 
    tiffilearr8 = tiff8.ReadAsArray()
    allpixels8 = (tiffilearr8 == 1) + (tiffilearr8 == 2) + (tiffilearr8== 3) + (tiffilearr8 == 4)
    snowpixels8 = tiffilearr8 == 4      
    
    # 7 value is allocated to pixels in built area class in Sentinel-2 10m Land Use/Land Cover timeseries  
    
    snow8 = np.extract(snowpixels8, tiffilearr8)
    totalarea8 = np.extract(allpixels8, tiffilearr8)
    area8 = (len(snow8)*1*1 ) / len(totalarea8) #since the spatial resolution of ims snow data is 1m
    
    print(f'relative snow covered area in the NWH: {area8} ') 
    snowarea8.append(area8)
snowarea9= []
def snowcoverarea9(path):
    tiff9 = gdal.Open(path)
#     lulc.GetProjection() 
    tiffilearr9 = tiff9.ReadAsArray()
    allpixels9 = (tiffilearr9 == 1) + (tiffilearr9 == 2) + (tiffilearr9 == 3) + (tiffilearr9 == 4)
    snowpixels9 = tiffilearr9 == 4      
    
    # 7 value is allocated to pixels in built area class in Sentinel-2 10m Land Use/Land Cover timeseries  
    
    snow9 = np.extract(snowpixels9, tiffilearr9)
    totalarea9 = np.extract(allpixels9, tiffilearr9)
    area9 = (len(snow9)*1*1 ) / len(totalarea9) #since the spatial resolution of ims snow data is 1m
    
    print(f'relative snow covered area in the NWH : {area9}') 
    snowarea9.append(area9)
snowarea10= []
def snowcoverarea10(path):
    tiff10 = gdal.Open(path)
#     lulc.GetProjection() 
    tiffilearr10 = tiff10.ReadAsArray()
    allpixels10 = (tiffilearr10 == 1) + (tiffilearr10 == 2) + (tiffilearr10 == 3) + (tiffilearr10 == 4)
    snowpixels10 = tiffilearr10 == 4      
    
    # 7 value is allocated to pixels in built area class in Sentinel-2 10m Land Use/Land Cover timeseries  
    
    snow10 = np.extract(snowpixels10, tiffilearr10)
    totalarea10 = np.extract(allpixels10, tiffilearr10)
    area10 = (len(snow10)*1*1 ) / len(totalarea10) #since the spatial resolution of ims snow data is 1m
    
    print(f'relative snow covered area in the NWH: {area10} ') 
    snowarea10.append(area10)
for i in folder1:
    print(f'For {i}:')
    print(snowcoverarea1(i))
for i in folder2:
    print(f'For {i}:')
    print(snowcoverarea2(i))
for i in folder3:
    print(f'For {i}:')
    print(snowcoverarea3(i))
for i in folder4:
    print(f'For {i}:')
    print(snowcoverarea4(i))
for i in folder5:
    print(f'For {i}:')
    print(snowcoverarea5(i))
for i in folder6:
    print(f'For {i}:')
    print(snowcoverarea6(i))
for i in folder7:
    print(f'For {i}:')
    print(snowcoverarea7(i))
for i in folder8:
    print(f'For {i}:')
    print(snowcoverarea8(i))
for i in folder9:
    print(f'For {i}:')
    print(snowcoverarea9(i))
for i in folder10:
    print(f'For {i}:')
    print(snowcoverarea10(i))
    
for i in range (0,9):
    snowcoverarea
# **Plotting Snow cover Area**
years = ['2023-24','2022-23', '2021-22','2020-21','2019-20','2018-19','2017-18','2016-17','2015-16', '2014-15' ]
plt.figure(figsize=(15,5))
plt.plot(dates,snowarea1)
plt.plot(dates,snowarea2)
plt.plot(dates,snowarea3)
plt.plot(dates,snowarea4)
plt.plot(dates,snowarea5)
plt.plot(dates,snowarea6)
plt.plot(dates,snowarea7)
plt.plot(dates,snowarea8)
plt.plot(dates,snowarea9)
plt.tick_params(axis='x', rotation=55, labelsize=7)
plt.title('Relative Snow cover area of NWH in Dec & Jan from 2015 to 2024')
plt.xlabel('Dates', fontsize=10)
plt.xticks(rotation=45)
plt.xticks(dates[::2])    # to plot 
plt.ylabel('Relative Snow Cover Area')
plt.legend(years,loc = 'upper center', fontsize='small')
# plt.savefig('F:\IIRS\Course material\Mod 4 - Scientific Geocomputing\Project\Data\Output\_area.png')
plt.show()

# =============================================================================
# Write to csv (organise/transpose in excel)
# =============================================================================

import csv
 
# field names 
fields = ['Date'] + years 
   
# data rows of csv file
rows=[snowarea1,snowarea2,snowarea3,snowarea4,snowarea5,snowarea6,snowarea7,snowarea8,snowarea9,snowarea10]
    
with open(r'E:\Rachit\Snow\SnowNWH.csv', 'w', newline='') as f:
         # using csv.writer method from CSV package
    write = csv.writer(f)
    write.writerow(fields)
    
    for i in range(0,10):
        write.writerow(rows[i])