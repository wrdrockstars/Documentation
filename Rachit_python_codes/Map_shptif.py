# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 14:24:26 2024

@author: gaurishs
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Feb  9 12:46:41 2024

@author: gaurishs
"""

# Install required libraries
# =============================================================================
# import subprocess
# import sys
# 
# def install(package):
#     subprocess.check_call([sys.executable, "-m", "pip", "install", package])
# 
# # Install required libraries
# install("geopandas")
# install("rasterio")
# install("matplotlib")
# =============================================================================

# Main code starts here
import calendar
import geopandas as gpd
import numpy as np
import rasterio
from rasterio.plot import show
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar
#from matplotlib.axis import Axis 
#import matplotlib.ticker as ticker
from sklearn.metrics.pairwise import haversine_distances
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import os
import glob

# Replace these paths with the actual paths to your shapefile and raster folder
shapefile_path = r"\\192.168.9.63\Share\NWH\shapefiles_watershed\Aglar.shp"
raster_folder = r"E:\Rachit\FLDAS_forcodetest\output\Snow_depth\Aglar"
#"\\192.168.9.63\Share\NWH\Montly_average\SM_cliped\Aglar"Asan
#"\\192.168.9.63\Share\NWH\Montly_average\Temp_clipped\Aglar"
#"\\192.168.9.63\Share\NWH\Chirps\Aglar_Chirps"

# Load shapefile
gdf = gpd.read_file(shapefile_path)


global_min = np.inf
global_max = -np.inf

#files = glob.glob(infol + '/**/*.tif',recursive=True)
raster_files =  glob.glob(raster_folder + '/**/*.tif',recursive=True)

#raster_files=[r"\\192.168.9.63\Share\NWH\Chirps\Aglar_Chirps\9_Chirps_precipitation_sum_Mean_4326_mmAnnual_mm.tif"]
# Iterate over each GeoTIFF file
for file in raster_files:
    raster_path =file #os.path.join(raster_folder,file)
    # Open the GeoTIFF file
    with rasterio.open(raster_path) as src:
        # Read the raster data
        data = src.read(1) # Assuming it's a single band raster
        raster_extent = src.bounds
        extent=np.array(raster_extent)

        # Update global minimum and maximum values
        local_min = np.nanmin(data)
        local_max = np.nanmax(data)
        print(local_max,local_min)
        if local_min < global_min:
            global_min = local_min
        if local_max > global_max:
            global_max = local_max


# Get a list of all raster files in the folder


# Loop through each raster file
for raster_file in raster_files:
    # Construct the full path to the current raster file
    current_raster_path =file #os.path.join(raster_folder, raster_file)
    
  
    # Load the current raster file
    raster = rasterio.open(current_raster_path)

    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(10, 10))
    
    #my_cmap_r = reverse_colourmap(my_cmap)

    # Plot raster
    show(raster, ax=ax, cmap='RdYlBu',vmin=global_min,vmax=global_max)#global_max,global_min

    # Plot shapefile
    gdf.plot(ax=ax, facecolor='none', edgecolor='red')

    # Add grid lines with degree-second formatter
    ax.grid(True, linestyle='--', linewidth=0.01, alpha=0.7, which='both')
    

    # Set the grid labels in degree-second format for latitude
    ax.xaxis.set_major_formatter('{x:.2f}°E')

    # Set the grid labels in degree-second format for longitude
    ax.yaxis.set_major_formatter('{x:.2f}°N')
    ax.yaxis.set_major_locator(MultipleLocator(0.05))
    ax.xaxis.set_major_locator(MultipleLocator(0.05))

    #tick_positions = ax.get_yticks()
    ax.tick_params(axis='y', rotation=90,pad=5)
    ax.tick_params(axis='x', pad=10)
    #ax.set_yticks(tick_positions +(-0.0))
  
    # Create custom legend with proxy artists and color bar
# =============================================================================
#     legend_handles = [plt.Line2D([0], [0], color='red', lw=2,label='Parbati Watershed' ),#label='Shapefile'
#                       plt.Line2D([0], [0], color='blue', lw=2, )]#label=f'Raster: {raster_file}')
#     
#     ax.legend(handles=legend_handles, loc='upper right')
# =============================================================================

    rasters_file=raster_file.split('\\')
    mon=rasters_file[-1].split('_')[0]
    var='Snow depth (m)'
    title='Annual'#str(calendar.month_name[int(mon)])#+'\n'+var
    plt.title(label=title,loc='center',fontsize=20)
    
    A=[extent[0]*np.pi/180.,extent[1]*np.pi/180.] #Latitude of interest here 40.7 deg, longitude -74.5
    B=[(extent[0]+1)*np.pi/180.,extent[1]*np.pi/180.] ##Latitude of interest here 40.7 deg, longitude -74.5+1
    dx=(6371000)*haversine_distances([A,B])[0,1]
    scalebar = ScaleBar(dx, 'm',border_pad=0.5,frameon=False,location='lower left')#location='lower left'
    #scalebar.scale_formatter = lambda value, unit: f'{value:f} {unit}'
    ax.add_artist(scalebar)
    
    
    # Add legend
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2.5%", pad=-0.02)#pad=-0.35
    cb = plt.colorbar(ax.get_images()[0], cax=cax, orientation='vertical', extend='both')
    cb.set_label(var)#r+' (% fraction)'
             
    #plt.yticks(rotation=90)
    #Set title (use the current raster file name as part of the title)
    dstname=current_raster_path[:-4]+".png"
    plt.savefig(dstname,dpi=600)
    

    # Show the plot
    plt.show()