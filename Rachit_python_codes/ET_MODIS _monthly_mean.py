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
geometry=ee.FeatureCollection('projects/ee-gsinghal572/assets/Tawi').geometry()

def funsum(m) :
  filtered = dataset.filter(ee.Filter.calendarRange(m, m, 'month'));
  sum1 = filtered.reduce(ee.Reducer.sum());
  return sum1.set('month', m);


start = ee.Date('2000-01-01');
end = ee.Date('2023-01-01'); #select end date according to dataset

# =============================================================================
# "NASA/FLDAS/NOAH01/C/GL/M/V001"; ET: Evap_tavg; kg/m^2/s; SM0-10:SoilMoi00_10cm_tavg;SoilMoi10_40cm_tavg; Snow depth:SnowDepth_inst; 	
#Snow water equivalent: SWE_inst;kg/m^2; Near surface air temperature:Tair_f_tavg; K Scale:11132
# =============================================================================
dataset = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR').filterDate(start,end);#
#for loading CHIRPS or ERA5 or GPM data just add: "UCSB-CHG/CHIRPS/DAILY", 'ECMWF/ERA5_LAND/DAILY_AGGR',"NASA/GPM_L3/IMERG_V06","" respt..
band='runoff_sum'# for CHIRPS or ERA5 or GPM use band: precipitation,total_precipitation_sum ,precipitationCal, , respt..
#snow_depth_water_equivalent,snow_depth,snow_cover,potential_evaporation_sum,runoff_sum,total_precipitation_sum, temperature_2m,volumetric_soil_water_layer_1, ERA5
#NDSI_Snow_Cover: MODIS/061/MOD10A1" Terra MODIS
#ET : MODIS/061/MOD16A2GF
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
precipitation=np.array(ee.Image(evi2weekList.get(gsinghal572)).select(band).sampleRectangle(region=geometry).getInfo())
print(precipitation.shape)
plt.imshow(precipitation)
p=precipitation.item()

'''
for i in range (0,12):
# precipitationsum=ee.Image(pptsum.get(i))# //note index 0 is the first image
 precipitation=ee.Image(evi2weekList.get(i)) #//note index 0 is the first image
 month=i+1
 print(month)
 #
 task = ee.batch.Export.image.toDrive(image = precipitation.select(band).multiply(1000).clip(geometry), #   
                                     region=geometry,  # an ee.Geometry object.
                                     description= str(month)+'_ERA5_Runoff_'+band+'_Mean_4326_mm',
                                     folder='Tawi_ERA_Runoff',
                                     scale=500,# 11132 for GPM, 5566 for CHIRPS, 11132 for ERA5
                                     crs='EPSG:4326',
                                     maxPixels=1e13)
 task.start()
#task.status()  



