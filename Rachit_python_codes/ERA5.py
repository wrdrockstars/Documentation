# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 09:45:00 2024

@author: rachit
"""
# =============================================================================
# #  Step-by-step guide for ERA5 precipitation hourly to daily conversion
# (https://confluence.ecmwf.int/display/CKB/ERA5%3A+How+to+calculate+daily+total+precipitation)
# 
# ### Script 1 - Use script below to download daily total precipitation ERA5 data for 1st and 2nd January 2017. This script will download total precipitation, in hourly steps, from CDS (Climate Data Store). Notice to cover total precipitation for 1st January 2017, we need two days of data.
# 
#   a. 1st January 2017 time = 01 - 23  will give you total precipitation data to cover 00 - 23 UTC for 1st January 2017
# 
#   b. 2nd January 2017 time = 00 will give you total precipitation data to cover 23 - 24 UTC for 1st January 2017
# =============================================================================


# Loop for precipitation timeseries daily data download (year by year)
import xarray as xr
import os
os.chdir(r'E:\\Rachit\\IGBP_Project\\Data\\Weather_data\\ERA5\\total_precip\\CDSAPI')
os.getcwd()
#!/usr/bin/env python
"""
Save as get-daily-tp.py and run "python get-daily-tp.py".

Input: None
Output: Daily total precipitation data from 1999 to 2001.
"""
import cdsapi
from datetime import datetime, timedelta
c = cdsapi.Client()

# Define the years you want to download data for
years = list(range(2002, 2024))

# Loop through the years
for year in years:
#     Loop through the months (from January to December)
    for month in range(1, 13):
#         Calculate the number of days in the current month
        if month in [1, 3, 5, 7, 8, 10, 12]:
            days_in_month = 31
        elif month in [4, 6, 9, 11]:
            days_in_month = 30
        else:
    #For leap years
            if year % 4 == 0 and (year % 100 != 0 or year % 400 == 0):
        # Leap year
                days_in_month = 29
            else:
                days_in_month = 28

# Loop through the days
        for day in range(1, days_in_month + 1):
            # Create the date string in YYYYMMDD format
            date_str = f"{year}{month:02d}{day:02d}"
            next_day = day + 1 if day < days_in_month else 1
            # to make the next day according to the days in month
            if day == days_in_month and month != 12:
                # print('Last day of the month',day)
        #             next_month = month + 1 if month < 12 else 1
        #             next_year = year + 1 if month == 12 else year
            # Retrieve data for the current day
                month1 = month
                month2 = month+1
                day1 = day
                day2 = 1
                r1 = c.retrieve(
                    'reanalysis-era5-single-levels', {
                        'variable': 'total_precipitation',
                        'product_type': 'reanalysis',
                        'year' : str(year),
                        'month': f"{month1:02d}",
                        'day'  : [f"{day:02d}"],
                        'time' : [
                            '00:00','01:00','02:00',
                            '03:00','04:00','05:00',
                            '06:00','07:00','08:00',
                            '09:00','10:00','11:00',
                            '12:00','13:00','14:00',
                            '15:00','16:00','17:00',
                            '18:00','19:00','20:00',
                            '21:00','22:00','23:00'
                        ],
                        'format': 'netcdf',
                        'area': [32.125, 72.125, 21.375, 89.875],  # Area bounds/extent
                    })
                r2 = c.retrieve(
                    'reanalysis-era5-single-levels', {
                        'variable': 'total_precipitation',
                        'product_type': 'reanalysis',
                        'year' : str(year),
                        'month': f"{month2:02d}",
                        'day'  : [f"{day2:02d}"],
                        'time' : [
                            '00:00','01:00','02:00',
                            '03:00','04:00','05:00',
                            '06:00','07:00','08:00',
                            '09:00','10:00','11:00',
                            '12:00','13:00','14:00',
                            '15:00','16:00','17:00',
                            '18:00','19:00','20:00',
                            '21:00','22:00','23:00'
                        ],
                        'format': 'netcdf',
                        'area': [32.125, 72.125, 21.125, 90.125],  # Area bounds/extent
                    })

                # Define the output file name
                d = datetime.strptime(str(date_str), '%Y%m%d')
                output_file1 = f'tp_{date_str}.nc'
                output_file2 = f'tp_{(d + timedelta(days=1)).strftime("%Y%m%d")}.nc'

                # Download the data
                r1.download(output_file1)
                r2.download(output_file2)
        #                 print(f'Downloaded data for {date_str} and {(d + timedelta(days=1)).strftime("%Y%m%d")}')


                #Combing the two files
                ds1 = xr.open_dataset(f'E:\\Rachit\\IGBP_Project\\Data\\Weather_data\\ERA5\\total_precip\\CDSAPI\\tp_{date_str}.nc')
                ds2 = xr.open_dataset(f'E:\\Rachit\\IGBP_Project\\Data\\Weather_data\\ERA5\\total_precip\\CDSAPI\\tp_{(d + timedelta(days=1)).strftime("%Y%m%d")}.nc')


                ds_combined = xr.concat([ds1, ds2], dim='time')

                # Save the combined dataset to a new NetCDF file
                ds_combined.to_netcdf(f'E:\\Rachit\\IGBP_Project\\Data\\Weather_data\\ERA5\\total_precip\\CDSAPI\\tp_{date_str}-{(d + timedelta(days=1)).strftime("%Y%m%d")}.nc')

        #                 ds_combined.download(f'E:\\Rachit\\IGBP_Project\\Data\\Weather_data\\ERA5\\total_precip\\ERA5_website_code_4daily_pcp\\correct2\\tp_{date_str}-{(d + timedelta(days=1)).strftime("%Y%m%d")}.nc')

                # Close the datasets
                ds1.close()
                ds2.close()

        #                 print(f'Downloaded data for {date_str}-{(d + timedelta(days=1)).strftime("%Y%m%d")}')

            elif day == days_in_month and month == 12:
                month1 = month
                month2 = 1
                day1 = day
                day2 = 1
                year1 = year
                year2 = year + 1
                r1 = c.retrieve(
                    'reanalysis-era5-single-levels', {
                        'variable': 'total_precipitation',
                        'product_type': 'reanalysis',
                        'year' : str(year1),
                        'month': f"{month1:02d}",
                        'day'  : [f"{day:02d}"],
                        'time' : [
                            '00:00','01:00','02:00',
                            '03:00','04:00','05:00',
                            '06:00','07:00','08:00',
                            '09:00','10:00','11:00',
                            '12:00','13:00','14:00',
                            '15:00','16:00','17:00',
                            '18:00','19:00','20:00',
                            '21:00','22:00','23:00'
                        ],
                        'format': 'netcdf',
                        'area': [32.125, 72.125, 21.125, 90.125],  # Area bounds/extent
                    })
                r2 = c.retrieve(
                    'reanalysis-era5-single-levels', {
                        'variable': 'total_precipitation',
                        'product_type': 'reanalysis',
                        'year' : str(year2),
                        'month': f"{month2:02d}",
                        'day'  : [f"{day2:02d}"],
                        'time' : [
                            '00:00','01:00','02:00',
                            '03:00','04:00','05:00',
                            '06:00','07:00','08:00',
                            '09:00','10:00','11:00',
                            '12:00','13:00','14:00',
                            '15:00','16:00','17:00',
                            '18:00','19:00','20:00',
                            '21:00','22:00','23:00'
                        ],
                        'format': 'netcdf',
                        'area': [32.125, 72.125, 21.125, 90.125],  # Area bounds/extent
                    })

                # Define the output file name
                d = datetime.strptime(str(date_str), '%Y%m%d')
                output_file1 = f'tp_{date_str}.nc'
                output_file2 = f'tp_{(d + timedelta(days=1)).strftime("%Y%m%d")}.nc'

                # Download the data
                r1.download(output_file1)
                r2.download(output_file2)
        #                 print(f'Downloaded data for {date_str} and {(d + timedelta(days=1)).strftime("%Y%m%d")}')


                #Combing the two files
                ds1 = xr.open_dataset(f'E:\\Rachit\\IGBP_Project\\Data\\Weather_data\\ERA5\\total_precip\\CDSAPI\\tp_{date_str}.nc')
                ds2 = xr.open_dataset(f'E:\\Rachit\\IGBP_Project\\Data\\Weather_data\\ERA5\\total_precip\\CDSAPI\\tp_{(d + timedelta(days=1)).strftime("%Y%m%d")}.nc')


                ds_combined = xr.concat([ds1, ds2], dim='time')

                # Save the combined dataset to a new NetCDF file
                ds_combined.to_netcdf(f'E:\\Rachit\\IGBP_Project\\Data\\Weather_data\\ERA5\\total_precip\\CDSAPI\\tp_{date_str}-{(d + timedelta(days=1)).strftime("%Y%m%d")}.nc')

        #                 ds_combined.download(f'E:\\Rachit\\IGBP_Project\\Data\\Weather_data\\ERA5\\total_precip\\ERA5_website_code_4daily_pcp\\correct2\\tp_{date_str}-{(d + timedelta(days=1)).strftime("%Y%m%d")}.nc')

                # Close the datasets
                ds1.close()
                ds2.close()

        #                 print(f'Downloaded data for {date_str}-{(d + timedelta(days=1)).strftime("%Y%m%d")}')
            else:
                r = c.retrieve(
                    'reanalysis-era5-single-levels', {
                        'variable': 'total_precipitation',
                        'product_type': 'reanalysis',
                        'year' : str(year),
                        'month': f"{month:02d}",
                        'day'  : [f"{day:02d}", f"{next_day:02d}"],
                        'time' : [
                            '00:00','01:00','02:00',
                            '03:00','04:00','05:00',
                            '06:00','07:00','08:00',
                            '09:00','10:00','11:00',
                            '12:00','13:00','14:00',
                            '15:00','16:00','17:00',
                            '18:00','19:00','20:00',
                            '21:00','22:00','23:00'
                        ],
                        'format': 'netcdf',
                        'area': [32.125, 72.125, 21.125, 90.125],  # Area bounds/extent
                    })

                # Define the output file name
                d = datetime.strptime(str(date_str), '%Y%m%d')
                output_file = f'tp_{date_str}-{(d + timedelta(days=1)).strftime("%Y%m%d")}.nc'

                # Download the data
                r.download(output_file)
        #                 print(f'Downloaded data for {date_str}-{(d + timedelta(days=1)).strftime("%Y%m%d")}')

        # print('Data download complete.')
        


# =============================================================================
# Time Correction for Indian Scenario UTC +5:30 (not done)
# =============================================================================
# =============================================================================
#### Visualizing the data
# import glob
# import matplotlib.dates as mdates
# import matplotlib.pyplot as plt
# import xarray as xr
# 
# data = xr.open_dataset(r"\\172.16.20.21\wrd_p_igbp\IGBP_PROJECT_2023-24\IGBP_Project_files\Data\Weather_data\ERA5\total_precip\pcp2011_2015.nc")
# print(data)
# subset = data.sel(time=slice('2013-06-13', '2013-06-17'))
# subset
# fig, ax = plt.subplots()
# fig.set_size_inches(50, 5)
# plt.plot(subset.variables['time'][:], subset.variables['tp'][:,5,5] )
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
# ax.xaxis.set_major_locator(mdates.HourLocator(interval=1))
# 
# plt.xticks(rotation=45)  # Rotate x-axis labels for better visibility
# plt.show()
# =============================================================================

# =============================================================================
# # For Time Correction
# import numpy as np
# 
# # Define the timedelta
# delta = np.timedelta64(5, 'h') + np.timedelta64(30, 'm')
# 
# # Add the timedelta to each value in the time variable
# new_time = time + delta
# 
# print(new_time)
# 
# # Update the time coordinate in the dataset
# data = data.assign_coords(time=new_time)
# 
# # Save the updated dataset back to the NetCDF file
# data.to_netcdf(r"E:\Rachit\IGBP_Project\Data\Weather_data\ERA5\total_precip\CDSAPI\tp_19900101-19900102_upd.nc")
# =============================================================================

# =============================================================================
# Script - 2 Run this second script to calculate daily total precipitation. All it does is to add up 24 values for a given day as described in step 1.
# =============================================================================
import os
os.chdir(r'E:\Rachit\IGBP_Project\Data\Weather_data\ERA5\total_precip\CDSAPI')
os.getcwd()
#!/usr/bin/env python 
"""
Save as file calculate-daily-tp.py and run "python calculate-daily-tp.py".

Input files: tp_YYYYMMDD-YYYYMMDD.nc
Output files: daily-tp_YYYYMMDD.nc
"""

import time, sys
from datetime import datetime, timedelta
from netCDF4 import Dataset, date2num, num2date
import numpy as np

# Iterate through the years 
for year in range(2000, 2025):
    # Loop through the months (from January to December)
    for month in range(1, 13):
        # Calculate the number of days in the current month
        if month in [1, 3, 5, 7, 8, 10, 12]:
            days_in_month = 31
        elif month in [4, 6, 9, 11]:
            days_in_month = 30
        else:
            # February
            if year % 4 == 0 and (year % 100 != 0 or year % 400 == 0):
                # Leap year
                days_in_month = 29
            else:
                days_in_month = 28
        
        # Loop through the days
        for day in range(1, days_in_month + 1):
            # Create the date string in YYYYMMDD format
            date_str = f"{year}{month:02d}{day:02d}"
##############################################################################################
            try:
                if day == days_in_month and month != 12:
#                   print('Last day of the month',day)
                    next_day = 1
                    day1 = f"{year}{month:02d}{day:02d}"
                    day2 = f"{year}{month+1:02d}{next_day:02d}"
                    d = datetime.strptime(str(day1), '%Y%m%d')
                    f_in = f'tp_%s-%s.nc' % (day1, day2)  # OR ----> f_in = f'tp_{date_str}-{(d + timedelta(days=1)).strftime("%Y%m%d")}.nc'
                    f_out = f'daily-tp_%s.nc' % day1      # OR ---->  f_out = f'daily-tp_{date_str}.nc'
                elif day == days_in_month and month == 12:
                    next_day = 1
                    nextyearmonth = 1
                    year1 = year
                    year2 = year + 1
                    day1 = f"{year1}{month:02d}{day:02d}"
                    day2 = f"{year2}{nextyearmonth:02d}{next_day:02d}"
                    d = datetime.strptime(str(day1), '%Y%m%d')
                    f_in = f'tp_%s-%s.nc' % (day1, day2)
                    f_out = f'daily-tp_%s.nc' % day1     
                else:
                    # print('Not Last day of the month',day)
                    day = f"{year}{month:02d}{day:02d}"          # Convert year, month, and day to string with zero padding
                    d = datetime.strptime(str(day), '%Y%m%d')
                    f_in = f'tp_%s-%s.nc' % (day, (d + timedelta(days = 1)).strftime('%Y%m%d'))  
                    f_out = f'daily-tp_%s.nc' % day      # OR ---->  f_out = f'daily-tp_{date_str}.nc'
                
                #  % symbol is used as a string formatting operator in Python. 
                #It is used to format strings by inserting values into a string template. 
                #This process is called string interpolation. 
                # %d and %s are placeholders for numerical and string values, respectively
                

                time_needed = []
                for i in range(1, 25):
                    time_needed.append(d + timedelta(hours = i))

                with Dataset(f_in) as ds_src:
                    var_time = ds_src.variables['time']
                    time_avail = num2date(var_time[:], var_time.units,
                            calendar = var_time.calendar)

                    indices = []
                    for tm in time_needed:
                        a = np.where(time_avail == tm)[0]
                        if len(a) == 0:
                            sys.stderr.write('Error: precipitation data is missing/incomplete - %s!\n'
                                    % tm.strftime('%Y%m%d %H:%M:%S'))
                            sys.exit(200)
                        else:
                            print('Found %s' % tm.strftime('%Y%m%d %H:%M:%S'))
                            indices.append(a[0])

                    var_tp = ds_src.variables['tp']
                    tp_values_set = False
                    for idx in indices:
                        if not tp_values_set:
                            data = var_tp[idx, :, :]
                            tp_values_set = True
                        else:
                            data += var_tp[idx, :, :]

                    with Dataset(f_out, mode = 'w', format = 'NETCDF3_64BIT_OFFSET') as ds_dest:
                        # Dimensions
                        for name in ['latitude', 'longitude']:
                            dim_src = ds_src.dimensions[name]
                            ds_dest.createDimension(name, dim_src.size)
                            var_src = ds_src.variables[name]
                            var_dest = ds_dest.createVariable(name, var_src.datatype, (name,))
                            var_dest[:] = var_src[:]
                            var_dest.setncattr('units', var_src.units)
                            var_dest.setncattr('long_name', var_src.long_name)

                        ds_dest.createDimension('time', None)
                        var = ds_dest.createVariable('time', np.int32, ('time',))
                        time_units = 'hours since 1900-01-01 00:00:00'
                        time_cal = 'gregorian'
                        var[:] = date2num([d], units = time_units, calendar = time_cal)
                        var.setncattr('units', time_units)
                        var.setncattr('long_name', 'time')
                        var.setncattr('calendar', time_cal)

                        # Variables
                        var = ds_dest.createVariable(var_tp.name, np.double, var_tp.dimensions)
                        var[0, :, :] = data
                        var.setncattr('units', var_tp.units)
                        var.setncattr('long_name', var_tp.long_name)

                        # Attributes
                        ds_dest.setncattr('Conventions', 'CF-1.6')
                        ds_dest.setncattr('history', '%s %s'
                                % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                                ' '.join(time.tzname)))
                
                # (Replace 'day = 19980101' with 'day = int(date_str)')

                # Existing logic for processing the NetCDF files goes here...

                print(f'Done! Daily total precipitation saved in {f_out}')

            except Exception as e:
                print(f'Error processing {date_str}: {str(e)}')
                
            except SystemExit as e:
                if e.code == 200:
                    print(f'Error: {e}')
                    continue  # Continue with the next iteration of the loop
                else:
                    raise  # Raise other SystemExit exceptions 
                    
# =============================================================================
# Convert to text
# =============================================================================
import numpy as np
import xarray as xr
import pandas as pd
import glob
import xarray as xr

# =============================================================================
# Concatenating nc files
# =============================================================================
# List of filenames to be combined
file_paths = glob.glob(r'E:\Rachit\IGBP_Project\Data\Weather_data\ERA5\total_precip\CDSAPI\daily-tp*.nc')             

# Open all files
datasets = [xr.open_dataset(file) for file in file_paths]

# Concatenate along the desired dimension (e.g., time, lat, lon)
combined_dataset = xr.concat(datasets, dim='time')  # Adjust 'time' with the appropriate dimension in your files

# Save the combined dataset to a new file
combined_dataset.to_netcdf(r'E:\Rachit\IGBP_Project\Data\Weather_data\ERA5\total_precip\CDSAPI\combined_daily_tp_1990_23.nc')

# ds = xr.open_mfdataset(r"E:\Rachit\IGBP_Project\Data\Weather_data\ERA5\total_precip\CDSAPI\daily-tp*.nc", combine='by_coords', concat_dim="time", engine='netcdf4') # open multiple files as a single dataset
# ============================================================================

# =============================================================================
# Converting to SWAT text
# =============================================================================

ncfile = xr.open_dataset(r"E:\Rachit\IGBP_Project\Data\Weather_data\ERA5\total_precip\CDSAPI\combined_daily_tp_1990_23.nc")
# latlon_values = ncfile.variables['latitude'][43]
# latlon_values


df = ncfile.to_dataframe().reset_index()           # convert netcdf to dataframe and to reset the indexes according to netcdf dimensions
df

# Replace very small non-zero values with zero
df['tp'] = df['tp'].apply(lambda x: 0 if abs(x) < 1e-10 else x)
# Replace nan values with zero
df['precip'] = df['tp'].apply(lambda x: 0 if pd.isna(x) else x) #

# non_null_rows = df.dropna(subset=['tp'])         # REMOVED NAN VALUES (values of coordinates within the extent of shapefile)

group1 = df[(df['latitude'] >= 21.375) & (df['latitude'] <= 32.125) & (df['longitude'] >= 72.125) & (df['longitude'] <= 89.875)]
group = group1[['latitude' , 'longitude']].drop_duplicates()  # To store lat lon pair in a separate dataframe and keep only unique combinations

# print(group.iloc[0]['latitude']) # To access the value of a particular entry in a particular column of a DataFrame
# print(group)


for i in range(len(group)):
    filtered_df = df[(df['latitude'] == group.iloc[i]['latitude']) & (df['longitude'] == group.iloc[i]['longitude'])]
    rain = filtered_df["tp"]*1000 # multiplying by 1000 to convert from meters to millimeters
    rain.to_csv(f"E:/Rachit/IGBP_Project/Data/Weather_data/ERA5/total_precip/CDSAPI/rain{group.iloc[i]['latitude']}_{group.iloc[i]['longitude']}.txt".format(i), header=['19900101'], index=None, mode='w')
print('Converted to SWAT text format')

# =============================================================================
# Pcp station input file preparation for SWAT
# =============================================================================
import geopandas as gpd
import pandas as pd

geometry = gpd.points_from_xy(group.longitude, group.latitude)
gdfele = gpd.GeoDataFrame(group, geometry=geometry, crs='EPSG:4326')

# save the GeoDataFrame as a shapefile
gdfele.to_file(r'E:\Rachit\IGBP_Project\Data\Weather_data\ERA5\total_precip\CDSAPI\tp_points.shp', driver='ESRI Shapefile')

# =============================================================================
# Now Go to ArcMap--> Extract values to point using DEM and saved elevattemp.shp
# =============================================================================

# shp to dataframe
eledf = gpd.read_file(r"E:\Rachit\IGBP_Project\Data\Weather_data\ERA5\total_precip\CDSAPI\tp_elevation.shp")
eledf
data = []
for i in range(len(group)):
    data.append({
    'ID': i+1,
    'NAME': f"rain{group.iloc[i]['latitude']}_{group.iloc[i]['longitude']}",
    'LAT': group.iloc[i]['latitude'],
    'LONG': group.iloc[i]['longitude'],
    'ELEVATION': eledf['RASTERVALU'][i]
    })
    datadf = pd.DataFrame(data)
    # (To match with APHRODITE station files) sort first by longitude and then by latitude, ensuring that longitude is sorted in ascending order and latitude is sorted within each longitude group
    datadf.sort_values(by=['LONG', 'LAT'], ascending=[True, True], inplace=True) 
    datadf.to_csv(f"E:\\Rachit\\IGBP_Project\\Data\\Weather_data\\ERA5\\total_precip\\CDSAPI\\tp_stations.txt".format(i), header=True, index=None, mode='w')
    
    
    
    
    
# =============================================================================
# For Temperature    
# =============================================================================
    
import cdsapi

c = cdsapi.Client()

r = c.retrieve(
    'reanalysis-era5-single-levels',
    {
        'product_type': 'reanalysis',
        'variable': [
            'maximum_2m_temperature_since_previous_post_processing', 'minimum_2m_temperature_since_previous_post_processing',
        ],
        'year': ['1993','1994','1995'],
        'month': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
        ],
        'day': [
            '01', '02', '03',
            '04', '05', '06',
            '07', '08', '09',
            '10', '11', '12',
            '13', '14', '15',
            '16', '17', '18',
            '19', '20', '21',
            '22', '23', '24',
            '25', '26', '27',
            '28', '29', '30',
            '31',
        ],
        'time': [
            '00:00', '01:00', '02:00',
            '03:00', '04:00', '05:00',
            '06:00', '07:00', '08:00',
            '09:00', '10:00', '11:00',
            '12:00', '13:00', '14:00',
            '15:00', '16:00', '17:00',
            '18:00', '19:00', '20:00',
            '21:00', '22:00', '23:00',
        ],
        'area': [
            32.125, 72.125, 21.375,
            89.875,
        ],
        'format': 'netcdf',
    })
r.download(r"E:\\Rachit\\IGBP_Project\\Data\\Weather_data\\ERA5\\Temp_2m\\CDSAPI\\temp_1993_1995.nc")
    

h = xr.open_dataset(r"E:\\Rachit\\IGBP_Project\\Data\\Weather_data\\ERA5\\Temp_2m\\CDSAPI\\temp_1990_1992.nc")
h1 = h.variables['mx2t'][0:25,1,1]
    
    
    
    
    
# =============================================================================
# #  Step-by-step guide for ERA5 precipitation hourly to daily conversion
# (https://confluence.ecmwf.int/display/CKB/ERA5%3A+How+to+calculate+daily+total+precipitation)
# 
# ### Script 1 - Use script below to download daily total precipitation ERA5 data for 1st and 2nd January 2017. This script will download total precipitation, in hourly steps, from CDS (Climate Data Store). Notice to cover total precipitation for 1st January 2017, we need two days of data.
# 
#   a. 1st January 2017 time = 01 - 23  will give you total precipitation data to cover 00 - 23 UTC for 1st January 2017
# 
#   b. 2nd January 2017 time = 00 will give you total precipitation data to cover 23 - 24 UTC for 1st January 2017
# =============================================================================


# =============================================================================
# # Loop for temperature timeseries daily data download (year by year)
# =============================================================================


import xarray as xr
import os
os.chdir(r"E:\Rachit\IGBP_Project\Data\Weather_data\ERA5\Temp_2m\CDSAPI")
os.getcwd()
#!/usr/bin/env python
"""
Input: None
Output: Min Max Temperature data from 1990 to 2023.
"""
import cdsapi
from datetime import datetime, timedelta
c = cdsapi.Client()

# Define the years you want to download data for
years = list(range(1999, 2024))

# Loop through the years
for year in years:
#     Loop through the months (from January to December)
    for month in range(1, 13):
#         Calculate the number of days in the current month
        if month in [1, 3, 5, 7, 8, 10, 12]:
            days_in_month = 31
        elif month in [4, 6, 9, 11]:
            days_in_month = 30
        else:
    #For leap years
            if year % 4 == 0 and (year % 100 != 0 or year % 400 == 0):
        # Leap year
                days_in_month = 29
            else:
                days_in_month = 28

# Loop through the days
        for day in range(1, days_in_month + 1):
            # Create the date string in YYYYMMDD format
            date_str = f"{year}{month:02d}{day:02d}"
            
            r1 = c.retrieve(
                'reanalysis-era5-single-levels', {
                    'variable': ['maximum_2m_temperature_since_previous_post_processing', 'minimum_2m_temperature_since_previous_post_processing'],
                    'product_type': 'reanalysis',
                    'year' : str(year),
                    'month': f"{month:02d}",
                    'day'  : [f"{day:02d}"],
                    'time' : [
                        '00:00','01:00','02:00',
                        '03:00','04:00','05:00',
                        '06:00','07:00','08:00',
                        '09:00','10:00','11:00',
                        '12:00','13:00','14:00',
                        '15:00','16:00','17:00',
                        '18:00','19:00','20:00',
                        '21:00','22:00','23:00'
                    ],
                    'format': 'netcdf',
                    'area': [32.125, 72.125, 21.375, 89.875],  # Area bounds/extent
                })

            # Define the output file name
            d = datetime.strptime(str(date_str), '%Y%m%d')
            output_file1 = f'temp_{date_str}.nc'
    

            # Download the data
            r1.download(output_file1)
            
            print(f'Downloaded data for {date_str}')

        print('Data download complete.')
        


# =============================================================================
# Time Correction for Indian Scenario UTC +5:30 (not done)
# =============================================================================
# =============================================================================
#### Visualizing the data
# import glob
# import matplotlib.dates as mdates
# import matplotlib.pyplot as plt
# import xarray as xr
# 
# data = xr.open_dataset(r"\\172.16.20.21\wrd_p_igbp\IGBP_PROJECT_2023-24\IGBP_Project_files\Data\Weather_data\ERA5\total_precip\pcp2011_2015.nc")
# print(data)
# subset = data.sel(time=slice('2013-06-13', '2013-06-17'))
# subset
# fig, ax = plt.subplots()
# fig.set_size_inches(50, 5)
# plt.plot(subset.variables['time'][:], subset.variables['tp'][:,5,5] )
# ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S'))
# ax.xaxis.set_major_locator(mdates.HourLocator(interval=1))
# 
# plt.xticks(rotation=45)  # Rotate x-axis labels for better visibility
# plt.show()
# =============================================================================

# =============================================================================
# # For Time Correction
# import numpy as np
# 
# # Define the timedelta
# delta = np.timedelta64(5, 'h') + np.timedelta64(30, 'm')
# 
# # Add the timedelta to each value in the time variable
# new_time = time + delta
# 
# print(new_time)
# 
# # Update the time coordinate in the dataset
# data = data.assign_coords(time=new_time)
# 
# # Save the updated dataset back to the NetCDF file
# data.to_netcdf(r"E:\Rachit\IGBP_Project\Data\Weather_data\ERA5\total_precip\CDSAPI\tp_19900101-19900102_upd.nc")
# =============================================================================

# =============================================================================
# Script - 2 Run this second script to calculate daily total precipitation. All it does is to add up 24 values for a given day as described in step 1.
# =============================================================================
import os
os.chdir(r'E:\Rachit\IGBP_Project\Data\Weather_data\ERA5\total_precip\CDSAPI')
os.getcwd()
#!/usr/bin/env python 
"""
Save as file calculate-daily-tp.py and run "python calculate-daily-tp.py".

Input files: tp_YYYYMMDD-YYYYMMDD.nc
Output files: daily-tp_YYYYMMDD.nc
"""

import time, sys
from datetime import datetime, timedelta
from netCDF4 import Dataset, date2num, num2date
import numpy as np

# Iterate through the years 
for year in range(2000, 2025):
    # Loop through the months (from January to December)
    for month in range(1, 13):
        # Calculate the number of days in the current month
        if month in [1, 3, 5, 7, 8, 10, 12]:
            days_in_month = 31
        elif month in [4, 6, 9, 11]:
            days_in_month = 30
        else:
            # February
            if year % 4 == 0 and (year % 100 != 0 or year % 400 == 0):
                # Leap year
                days_in_month = 29
            else:
                days_in_month = 28
        
        # Loop through the days
        for day in range(1, days_in_month + 1):
            # Create the date string in YYYYMMDD format
            date_str = f"{year}{month:02d}{day:02d}"
##############################################################################################
            try:
                if day == days_in_month and month != 12:
#                   print('Last day of the month',day)
                    next_day = 1
                    day1 = f"{year}{month:02d}{day:02d}"
                    day2 = f"{year}{month+1:02d}{next_day:02d}"
                    d = datetime.strptime(str(day1), '%Y%m%d')
                    f_in = f'tp_%s-%s.nc' % (day1, day2)  # OR ----> f_in = f'tp_{date_str}-{(d + timedelta(days=1)).strftime("%Y%m%d")}.nc'
                    f_out = f'daily-tp_%s.nc' % day1      # OR ---->  f_out = f'daily-tp_{date_str}.nc'
                elif day == days_in_month and month == 12:
                    next_day = 1
                    nextyearmonth = 1
                    year1 = year
                    year2 = year + 1
                    day1 = f"{year1}{month:02d}{day:02d}"
                    day2 = f"{year2}{nextyearmonth:02d}{next_day:02d}"
                    d = datetime.strptime(str(day1), '%Y%m%d')
                    f_in = f'tp_%s-%s.nc' % (day1, day2)
                    f_out = f'daily-tp_%s.nc' % day1     
                else:
                    # print('Not Last day of the month',day)
                    day = f"{year}{month:02d}{day:02d}"          # Convert year, month, and day to string with zero padding
                    d = datetime.strptime(str(day), '%Y%m%d')
                    f_in = f'tp_%s-%s.nc' % (day, (d + timedelta(days = 1)).strftime('%Y%m%d'))  
                    f_out = f'daily-tp_%s.nc' % day      # OR ---->  f_out = f'daily-tp_{date_str}.nc'
                
                #  % symbol is used as a string formatting operator in Python. 
                #It is used to format strings by inserting values into a string template. 
                #This process is called string interpolation. 
                # %d and %s are placeholders for numerical and string values, respectively
                

                time_needed = []
                for i in range(1, 25):
                    time_needed.append(d + timedelta(hours = i))

                with Dataset(f_in) as ds_src:
                    var_time = ds_src.variables['time']
                    time_avail = num2date(var_time[:], var_time.units,
                            calendar = var_time.calendar)

                    indices = []
                    for tm in time_needed:
                        a = np.where(time_avail == tm)[0]
                        if len(a) == 0:
                            sys.stderr.write('Error: precipitation data is missing/incomplete - %s!\n'
                                    % tm.strftime('%Y%m%d %H:%M:%S'))
                            sys.exit(200)
                        else:
                            print('Found %s' % tm.strftime('%Y%m%d %H:%M:%S'))
                            indices.append(a[0])

                    var_tp = ds_src.variables['tp']
                    tp_values_set = False
                    for idx in indices:
                        if not tp_values_set:
                            data = var_tp[idx, :, :]
                            tp_values_set = True
                        else:
                            data += var_tp[idx, :, :]

                    with Dataset(f_out, mode = 'w', format = 'NETCDF3_64BIT_OFFSET') as ds_dest:
                        # Dimensions
                        for name in ['latitude', 'longitude']:
                            dim_src = ds_src.dimensions[name]
                            ds_dest.createDimension(name, dim_src.size)
                            var_src = ds_src.variables[name]
                            var_dest = ds_dest.createVariable(name, var_src.datatype, (name,))
                            var_dest[:] = var_src[:]
                            var_dest.setncattr('units', var_src.units)
                            var_dest.setncattr('long_name', var_src.long_name)

                        ds_dest.createDimension('time', None)
                        var = ds_dest.createVariable('time', np.int32, ('time',))
                        time_units = 'hours since 1900-01-01 00:00:00'
                        time_cal = 'gregorian'
                        var[:] = date2num([d], units = time_units, calendar = time_cal)
                        var.setncattr('units', time_units)
                        var.setncattr('long_name', 'time')
                        var.setncattr('calendar', time_cal)

                        # Variables
                        var = ds_dest.createVariable(var_tp.name, np.double, var_tp.dimensions)
                        var[0, :, :] = data
                        var.setncattr('units', var_tp.units)
                        var.setncattr('long_name', var_tp.long_name)

                        # Attributes
                        ds_dest.setncattr('Conventions', 'CF-1.6')
                        ds_dest.setncattr('history', '%s %s'
                                % (datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
                                ' '.join(time.tzname)))
                
                # (Replace 'day = 19980101' with 'day = int(date_str)')

                # Existing logic for processing the NetCDF files goes here...

                print(f'Done! Daily total precipitation saved in {f_out}')

            except Exception as e:
                print(f'Error processing {date_str}: {str(e)}')
                
            except SystemExit as e:
                if e.code == 200:
                    print(f'Error: {e}')
                    continue  # Continue with the next iteration of the loop
                else:
                    raise  # Raise other SystemExit exceptions 
                    
# =============================================================================
# Convert to text
# =============================================================================
import numpy as np
import xarray as xr
import pandas as pd
import glob
import xarray as xr

# =============================================================================
# Concatenating nc files
# =============================================================================
# List of filenames to be combined
file_paths = glob.glob(r'E:\Rachit\IGBP_Project\Data\Weather_data\ERA5\total_precip\CDSAPI\daily-tp*.nc')             

# Open all files
datasets = [xr.open_dataset(file) for file in file_paths]

# Concatenate along the desired dimension (e.g., time, lat, lon)
combined_dataset = xr.concat(datasets, dim='time')  # Adjust 'time' with the appropriate dimension in your files

# Save the combined dataset to a new file
combined_dataset.to_netcdf(r'E:\Rachit\IGBP_Project\Data\Weather_data\ERA5\total_precip\CDSAPI\combined_daily_tp_1990_23.nc')

# ds = xr.open_mfdataset(r"E:\Rachit\IGBP_Project\Data\Weather_data\ERA5\total_precip\CDSAPI\daily-tp*.nc", combine='by_coords', concat_dim="time", engine='netcdf4') # open multiple files as a single dataset
# ============================================================================

ncfile = xr.open_dataset(r"E:\Rachit\IGBP_Project\Data\Weather_data\ERA5\total_precip\CDSAPI\combined_daily_tp_1990_23.nc")

df = ncfile.to_dataframe().reset_index()           # convert netcdf to dataframe and to reset the indexes according to netcdf dimensions
df

# Replace very small non-zero values with zero
df['tp'] = df['tp'].apply(lambda x: 0 if abs(x) < 1e-10 else x)

non_null_rows = df.dropna(subset=['tp'])         # REMOVED NAN VALUES (values of coordinates within the extent of shapefile)
non_null_rows

group = non_null_rows[['latitude' , 'longitude']].drop_duplicates() # To store lat lon pair in a separate dataframe and keep only unique combinations

# print(group.iloc[0]['latitude']) # To access the value of a particular entry in a particular column of a DataFrame
# print(group)


for i in range(len(group)):
    filtered_df = non_null_rows[(non_null_rows['latitude'] == group.iloc[i]['latitude']) & (non_null_rows['longitude'] == group.iloc[i]['longitude'])]
    rain = filtered_df["tp"]
    rain.to_csv(f"E:/Rachit/IGBP_Project/Data/Weather_data/ERA5/total_precip/CDSAPI/rain{group.iloc[i]['latitude']}_{group.iloc[i]['longitude']}.txt".format(i), header=['19900101'], index=None, mode='w')
print('Converted to SWAT text format')

# =============================================================================
# Pcp station input file preparation for SWAT
# =============================================================================
import geopandas as gpd
import pandas as pd

geometry = gpd.points_from_xy(group.longitude, group.latitude)
gdfele = gpd.GeoDataFrame(group, geometry=geometry, crs='EPSG:4326')

# save the GeoDataFrame as a shapefile
gdfele.to_file(r'E:\Rachit\IGBP_Project\Data\Weather_data\ERA5\total_precip\CDSAPI\tp_points.shp', driver='ESRI Shapefile')

# =============================================================================
# Now Go to ArcMap--> Extract values to point using DEM and saved elevattemp.shp
# =============================================================================

# shp to dataframe
eledf = gpd.read_file(r"E:\Rachit\IGBP_Project\Data\Weather_data\ERA5\total_precip\CDSAPI\tp_elevation.shp")
eledf
data = []
for i in range(len(group)):
    data.append({
    'ID': i+1,
    'NAME': f"rain{group.iloc[i]['latitude']}_{group.iloc[i]['longitude']}",
    'LAT': group.iloc[i]['latitude'],
    'LONG': group.iloc[i]['longitude'],
    'ELEVATION': eledf['RASTERVALU'][i]
    })
    datadf = pd.DataFrame(data)
    datadf.to_csv(f"E:\\Rachit\\IGBP_Project\\Data\\Weather_data\\ERA5\\total_precip\\CDSAPI\\tp_stations.txt".format(i), header=True, index=None, mode='w')
