# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 14:36:11 2024

@author: gaurishs
"""

import glob
import numpy as np
from osgeo import gdal

cpath=r"\\192.168.9.65\share_supriya\soil erosion\C"
rpath=r"\\192.168.9.65\share_supriya\soil erosion\R"
dst=r"\\192.168.9.65\share_supriya\soil erosion\RC"
lspath=r"\\192.168.9.65\share_supriya\soil erosion\LS\LS .tif"
kpath=r"\\192.168.9.65\share_supriya\soil erosion\K\K.tif"
ppath=r"\\192.168.9.65\share_supriya\soil erosion\P\P.tif"

ls=gdal.Open(lspath)
ls=ls.GetRasterBand(1).ReadAsArray()

k=gdal.Open(kpath)
k=k.GetRasterBand(1).ReadAsArray()

p=gdal.Open(ppath)
p=p.GetRasterBand(1).ReadAsArray()


for filnamc in glob.glob(cpath + '/*.tif'):
    #print(filnamc)    
    c = gdal.Open(filnamc)
    crs=c.GetGeoTransform()
    proj=c.GetProjection()
    c = c.GetRasterBand(1).ReadAsArray()
   
    for filnamr in glob.glob(rpath + '/*.tif'):
        #print(filnamr)     
        r = gdal.Open(filnamr)
        r= r.GetRasterBand(1).ReadAsArray()
        
        filc=filnamc.split('\\')
        filr=filnamr.split('\\')
    
        if (filr[-1]==filc[-1]):
            print(filnamr)
            print(filnamc)
            
            rc=r*c*ls*k*p 
             
            size=c.shape
            rows = size[0]  # number of pixels in x
            cols = size[1]
            driver = gdal.GetDriverByName('GTiff')
            #dst=dst+'\\'+filr[-1]
            # DataSet = driver.Create(outRaster, cols, rows, no_bands, gdal.GDT_Byte)
            dataset1 = driver.Create(dst+'\\'+filr[-1],cols,rows,1,gdal.GDT_Float32)
            dataset1.SetGeoTransform((crs))  
            dataset1.SetProjection(proj)
            dataset1.GetRasterBand(1).WriteArray(rc)
            dataset1.FlushCache()
            dataset1= None
        
   