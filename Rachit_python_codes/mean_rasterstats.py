# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 19:14:12 2024

@author: gaurishs
"""

import glob
import matplotlib.pyplot as plt
import geopandas as geopandas
from rasterstats import zonal_stats
import pandas as pd
import numpy as np

names=['Asan','Aglar','Gangotri','Hanval','Jhelum','Parbati','Suketi','Tawi','Beas_Manali']#','Aglar','Gangotri','Hanval','Jhelum','Parbati','Suketi','Tawi',Beas_Manali'
for name in names:
    
    #infol=r"D:\Share\NWH\Chirps\\"+name+"_Chirps"
    path = r"E:\Rachit\FLDAS_forcodetest\output\pcp"
    files = sorted(glob.glob(path+'\\*mmAnnual_mm.tif'))
   
    #files.sort(key=lambda x: int(x.split('\\')[-1].split('_')[0]))      # to sort properly according to numbers
    print(files)
    means1 = []
    for i in files:
        stats = zonal_stats(r"\\192.168.9.63\Share\NWH\shapefiles_watershed\\"+name+".shp", i, stats="mean") # we can use: min max min count sum std median majority minority unique range nodata percentile
        means1.append(stats)
        
        months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        list_of_values_Aglar = [item[0]['mean'] for item in means1]
    df1 = pd.DataFrame(list_of_values_Aglar)
    df1.insert(0, "Months", months, True)
    df1.to_csv(r"\\192.168.9.63\Share\NWH\Graphs\\"+name+"_pcp_means.csv")
    df1=None