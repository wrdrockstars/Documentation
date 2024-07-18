# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 11:41:43 2024

@author: rachit
"""

# =============================================================================
# For unzipping all files in a Folder
# =============================================================================
import zipfile
import gzip  # to unzip .gz files
import shutil


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
