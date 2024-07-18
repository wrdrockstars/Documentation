# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 16:11:14 2023

@author: Gaurish Singhal
"""
import ee
import numpy as np
import matplotlib.pyplot as plt
ee.Initialize()
#copy the asset id
geometry=ee.FeatureCollection('projects/ee-gsinghal572/assets/watershed').geometry()

def funsum(m) :
  filtered = dataset.filter(ee.Filter.calendarRange(m, m, 'month'));
  sum1 = filtered.reduce(ee.Reducer.sum());
  return sum1.set('month', m);


start = ee.Date('2023-01-01');
end = ee.Date('2024-01-01'); #select end date according to dataset
dataset = ee.ImageCollection('NASA/GPM_L3/IMERG_V06').filterDate(start,end);
#for loading CHIRPS or ERA5 or GPM data just add: "UCSB-CHG/CHIRPS/DAILY", 'ECMWF/ERA5/MONTHLY',"NASA/GPM_L3/IMERG_V06" respt..
band='precipitationCal'# for CHIRPS or ERA5 or GPM use band: precipitation,total_precipitation ,precipitationCal respt..

monthlysum = ee.ImageCollection(
  ee.List.sequence(1, 12).map(funsum));

no=end.difference(start,'year');
no=ee.Number(no).round();
print(no.getInfo())

suminfo=monthlysum.getInfo()

monthlymean = monthlysum.map(lambda image : image.divide(ee.Image(no)));
  
meaninfo= monthlymean.getInfo()

pptsum=monthlysum.toList(999)
evi2weekList=monthlymean.toList(999)

band=band+'_sum'

'''
# use this section for getting data for specific month
precipitation=np.array(ee.Image(evi2weekList.get(7)).select(band).sampleRectangle(region=geometry).getInfo())
print(precipitation.shape)
plt.imshow(precipitation)
p=precipitation.item()

'''
for i in range (0,12):
# precipitationsum=ee.Image(pptsum.get(i))# //note index 0 is the first image
 precipitation=ee.Image(evi2weekList.get(i))
# precipitation=precipitation.multiply(1000) #//note index 0 is the first image
 month=i+1
 #
 task = ee.batch.Export.image.toDrive(image = precipitation.select(band),                                     
                                     region=geometry,  # an ee.Geometry object.
                                     description= str(month)+'_GPM_Mean_Rainfall_mm',
                                     folder='Beesalpur',
                                     scale=500,# 11132 for GPM, 5566 for CHIRPS, 27830 for ERA5
                                     crs='EPSG:4326',
                                     maxPixels=1e13)
 task.start()
#task.status()  



