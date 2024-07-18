'''
# =============================================================================
# ----LISS 4 Image Processing using XArray and Dask
# (Source: https://spatialthoughts.com/2023/12/25/liss4-processing-xarray/)
# Converting the DN values to TOA Reflectance. We will use Python to process the downloaded scene into a stacked 3-band Cloud-Optimized GeoTIFF files with TOA Reflectance.
# =============================================================================
'''

# =============================================================================
# ----INSTALLATION AND CONFIGURATION
# Preferred method of installing and managing Python packages is using conda. We create a new environment and install the required packages. We are also using the PyEphem package to obtain the Earth-Sun distance that can be installed using pip.
# =============================================================================

# conda create --name liss4
# conda activate liss4
# conda install -c conda-forge rioxarray dask jupyterlab -y
# pip install ephem

import datetime
import ephem
import math
import os
import rioxarray as rxr
import xarray as xr
import zipfile
import rasterio as rio


# =============================================================================
# Dask is a parallell processing library
# =============================================================================
# Next we initiate a local dask cluster. This uses all the available cores on your machine in parallel. You will see a link to the Dask Dashboard. We will use this dashboard later in this tutorial.
from dask.distributed import Client, progress

client = Client()  # set up local cluster on the machine

# =============================================================================
# ----EXTRACT AND READ METADATA
# We first unzip the scene zip file to a folder. 
# But first ensure the directory path where you want to extract the zip file.
# =============================================================================

os.chdir("E:\Rachit\Extras\Liss4_zip")

with zipfile.ZipFile("E:\Rachit\Extras\Liss4_zip\R2F21DEC2023065764009600050SSANSTUC00GTDA.zip") as zf:
  # The LISS4 zip files contain a folder with all the data
  # Get the folder name
  foldername = [
      info.filename for info in zf.infolist()
      if info.is_dir()][0]
  # Extract all the data
  zf.extractall()
 
print(f'Extracted the files to {foldername}.')

# =============================================================================
# Each scene contains a file named BAND_META.txt containing the scene metadata. We parse the file and extract the data as a Python dictionary.
# =============================================================================

metadata_filename = 'BAND_META.txt'
metadata_filepath = os.path.join(foldername, metadata_filename)
 
metadata = {}
with open(metadata_filepath) as f:
  for line in f:
    key, value = line.split('=')
    metadata[key.strip()] = value.strip()
 
scene_id = metadata['OTSProductID']
print(f'Metadata extracted successfully for scene: {scene_id}')

# =============================================================================
# ----READ AND STACK BANDS
# LISS4 images come as 3 individual GeoTIFF rasters for each band. The image files are named BAND2.tif, BAND3.tif and BAND4.tif . We read them using rioxarray. Here we set chunks=True indicating that we want to use Dask and split the dataset into smaller chunks that can be processed in parallel.
# =============================================================================


# =============================================================================
# b2_path = os.path.join(foldername, 'BAND2.tif')
# b3_path = os.path.join(foldername, 'BAND3.tif')
# b4_path = os.path.join(foldername, 'BAND4.tif')
# =============================================================================
# Change the following paths according to the extracted folder location
b2_path = r"E:\Rachit\Extras\Liss4_zip\BH_R2F21DEC2023065764009600050SSANSTUC00GTDA\BAND2.tif"
b3_path = r"E:\Rachit\Extras\Liss4_zip\BH_R2F21DEC2023065764009600050SSANSTUC00GTDA\BAND3.tif"
b4_path = r"E:\Rachit\Extras\Liss4_zip\BH_R2F21DEC2023065764009600050SSANSTUC00GTDA\BAND4.tif"
 
b2_ds = rxr.open_rasterio(b2_path, chunks=True)
b3_ds = rxr.open_rasterio(b3_path, chunks=True)
b4_ds = rxr.open_rasterio(b4_path, chunks=True)

# =============================================================================
# Create a XArray Dataset by stacking individual band images. The scene has a NoData value of 0. So we set the correct NoData value before further processing.
# =============================================================================

scene = xr.concat(
    [b2_ds, b3_ds, b4_ds], dim='band').assign_coords(
    band=['BAND2', 'BAND3', 'BAND4'])



scene = scene.where(scene != 0)
scene.name = scene_id

# =============================================================================
# ----CONVERTING DN TO REFLECTANCES
# The pixel values of the source images are DN values that need to be converted to reflectances before they can be used for analysis.
# 
# The correction formulae and sensor parameters are published in the following paper
# 
# Sharma, Anu & Badarinath, K. & Roy, Parth. (2008). Corrections for atmospheric and adjacency effects on high resolution sensor data a case study using IRS-P6 LISS-IV data. 
# https://www.isprs.org/proceedings/xxxvii/congress/8_pdf/3_wg-viii-3/05.pdf
# =============================================================================

acq_date_str = metadata['DateOfPass']
# Date is in the format 04-MAR-2023
acq_date = datetime.datetime.strptime(
    acq_date_str, '%d-%b-%Y')
 
sun_elevation_angle = metadata['SunElevationAtCenter']
sun_zenith_angle = 90 - float(sun_elevation_angle)

# =============================================================================
# We need to compute the Earth-Sun distance at the date of acquisition. We use the pyephm library to obtain this
# =============================================================================

observer = ephem.Observer()
observer.date = acq_date
sun = ephem.Sun()
sun.compute(observer)
d = sun.earth_distance

# =============================================================================
# Define the Saturation Radiance for each band. These come from the RESOURCESAT-2 Data Users’ Handbook
# (https://www.euromap.de/download/R2_data_user_handbook.pdf') and are available in metadata as well.
# =============================================================================



b2_sr = 47.5614
b3_sr = 45.6990
b4_sr = 31.5000

# =============================================================================
# Define ex-atmospheric solar irradiance values for each band ESUN for each band 
# as published in At-sensor Solar Exo-atmospheric Irradiance, Rayleigh Optical 
# Thickness and Spectral parameters of RS-2 Sensors 
# ('https://www.researchgate.net/profile/Senthil-Kumar-135/post/What-is-the-formula-for-converting-DN-values-to-reflectance-value-for-IRS-R2-LISS-IV-image/attachment/59d61ff4c49f478072e97d4b/AS%3A271753849311232%401441802574228/download/RS2-Spectral+Characteristics.pdf')
# =============================================================================

b2_esun = 181.89
b3_esun = 156.96
b4_esun = 110.48

# =============================================================================
# Define other contants needed for computation.
# =============================================================================
pi = math.pi
sun_zenith_angle_rad = math.radians(sun_zenith_angle)

# =============================================================================
# ----Convert DN to Radiance
# =============================================================================

b2_dn = scene.sel(band='BAND2')
b3_dn = scene.sel(band='BAND3')
b4_dn = scene.sel(band='BAND4')
 
b2_rad = b2_dn*b2_sr/1024
b3_rad = b3_dn*b3_sr/1024
b4_rad = b4_dn*b4_sr/1024

# =============================================================================
# ----Convert Radiance to TOA Reflectance
# =============================================================================

b2_ref = (pi*b2_rad*d*d)/(b2_esun*math.cos(sun_zenith_angle_rad))
b3_ref = (pi*b3_rad*d*d)/(b3_esun*math.cos(sun_zenith_angle_rad))
b4_ref = (pi*b4_rad*d*d)/(b4_esun*math.cos(sun_zenith_angle_rad))

# =============================================================================
# Stack the bands into a single XArray Dataset.
# =============================================================================

reflectance_bands = [b2_ref, b3_ref, b4_ref]
scene_ref = xr.concat(
    reflectance_bands, dim='band').assign_coords(
    band=['BAND2', 'BAND3', 'BAND4']
).chunk('auto')
scene_ref.name = scene_id

# =============================================================================
# ----WRITE RESULTS TO DISK
# Our DataArray is structured to have ‘band’ as a dimension which makes it easy for data manipulation and processing. But for use in standard GIS software – it is better to create an XArray Dataset with each band as a variable.
# =============================================================================

output_ds = scene_ref.to_dataset('band')

# =============================================================================
# Define the options for the output file. We use the COG driver to create a Cloud-Optimized GeoTIFF file.
# =============================================================================

output_file = r"E:\Rachit\Extras\Liss4_zip\R2F21DEC2023065764009600050SSANSTUC00GTDA.tif"
 
output_options = {
    'driver': 'COG',
    'compress': 'deflate',
    'num_threads': 'all_cpus',
    'windowed': False # set True if you run out of RAM
}

# =============================================================================
# Write the raster.
# =============================================================================
output_ds[['BAND2', 'BAND3', 'BAND4']].rio.to_raster(
    output_file, **output_options)
print(f'Output file created {output_file}')
