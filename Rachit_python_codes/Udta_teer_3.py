import glob
import numpy as np
from datetime import datetime
import pandas as pd
import xarray
import geopandas
from shapely.geometry import mapping
from tqdm import tqdm
import gc

import argparse

def main():
    parser = argparse.ArgumentParser(description="A simple CLI tool for extracting fSCA from NC files")
    parser.add_argument("infol", type=str, help="Input path of the nf files without quotes")
    parser.add_argument("in_sf", type=str, help="Input path of the shapefile without quotes")    
    
    args = parser.parse_args()   


    input_folder = args.infol
    shape_file_path = args.in_sf
  
    # =============================================================================
    # Read shapefile using geopandas
    # =============================================================================
    sf = geopandas.read_file(shape_file_path)
    
    # =============================================================================
    # Pepare a list of paths of the input files
    # =============================================================================
    in_files = sorted(glob.glob(input_folder + "/*.NC"))
    
    # =============================================================================
    # Creating two empty lists to store date and snow covered area
    # =============================================================================
    date_col = []
    snow_area_col = []
    
    # =============================================================================
    # Looping through all the .nc files in the path to extract snow covered area after
    # clipping each one using a shapefile
    # =============================================================================
    for file in tqdm(in_files):
        # continue
    # =============================================================================
    # Generating data from the name of the file
    # =============================================================================
        data_date = file.split('/')[-1][-11:-3]
    # =============================================================================
    # Converting date to proper format
    # =============================================================================
        data_date = datetime.strptime(data_date,"%Y%m%d")    
    # =============================================================================
    # Appending dates to an empty list, to be concatenated to a pandas dataframe later
    # =============================================================================
        date_col.append(data_date.date())    
    # =============================================================================
    # Opening .nc dataset using xarray
    # =============================================================================
        data = xarray.open_dataset(file)
    # =============================================================================
    # Setting the projection of the shape file to .nc file for clipping later;
    # by default, there is no crs in .nc file
    # =============================================================================
        data = data.rio.write_crs(sf.crs, inplace = True)
        
        # data = data.where(np.isfinite(data),0)
    # =============================================================================
    #     This line has just been added to avoid warning message
    # =============================================================================
        data = data.where(~np.isnan(data), drop=False)
    # =============================================================================
    # clipping dataset using shape file; mapping seems to be a utility of shapely
    # which helps to clip nc file; here sf is shape file object opened using geopandas
    # .apply method is bing used to apply the mapping utility of shapely to retrive a 
    # geometry which is to be used to clip the .nc file; a crs is to be applied
    # =============================================================================
        clipped_data = data.rio.clip(sf.geometry.apply(mapping), sf.crs, all_touched = True)
    # =============================================================================
    # retrieving the actual data of fractional snow cover for our study area clipped
    # using the shape file
    # =============================================================================
        fsca = clipped_data['fSCA'][:][:]    
    # =============================================================================
    # Converting data from masked array above to a numpy ndarray
    # =============================================================================
        fsca = fsca.data
    # =============================================================================
    # Generating a binary for snow and non-snow areas
    # =============================================================================
        fsca_binary = np.where((fsca>50)&(fsca<=100),1,0)    
    # =============================================================================
    # Counting snow pixels
    # =============================================================================
        snow_pixels = np.count_nonzero(fsca_binary)
    # =============================================================================
    # Calculating the area of snow by multiplying with the area of each pixel and Converting
    # it to sq.km.
    # =============================================================================
        area = snow_pixels * 250000 / 1000000
    # =============================================================================
    # Appending calculated snow area to an empty list which will be concatenated to the
    # pandas dataframe along with date as generated previously
    # =============================================================================
        snow_area_col.append(area)




        
        data.close()
        del data, clipped_data, fsca, fsca_binary
        # Explicit garbage collection
        gc.collect()





    
    # =============================================================================
    # Generating a dictionary of date and snow area; one of the ways for preparing
    # data to be written in pandas dataframe
    # =============================================================================
    df_data = {'Date': date_col,
                      'Snow Area': snow_area_col}
    # =============================================================================
    # Writing data to pandas dataframe
    # =============================================================================
    df = pd.DataFrame(df_data)
    # =============================================================================
    # Converting pandas dataframe to csv and saving to disk
    # =============================================================================
    df.to_csv(input_folder + '/' + input_folder.split('/')[-1] +"_out_UK.csv")
    
    # =============================================================================
    # Setting all variables to none in order to release memory (RAM)
    # =============================================================================
    input_folder = in_files = date_col = snow_area_col = data_date = df = data = clipped_data = None


if __name__ == "__main__":
    main()
