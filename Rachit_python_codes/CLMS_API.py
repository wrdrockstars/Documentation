# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 17:05:56 2024

@author: rachit
"""
# https://github.com/eea/clms-api-docs/blob/develop/source/notebooks/clms_portal.ipynb

# =============================================================================
# The goal of this notebook is to help you kickstart your project and to provide an example on how to download a dataset. You can see all the API documentation here:
# https://eea.github.io/clms-api-docs/
# You can view information regarding the Machine2Machine download here:
# https://land.copernicus.eu/en/how-to-guides/download-spatial-data/m2m-download
# =============================================================================


# =============================================================================
# Sign up and token declaration
# =============================================================================
# Go to https://land.copernicus.eu/en
# Register using EU login
# After signing in to the portal navigate to 'My profile' under your username
# Create a new token, press Copy to clipboard and paste it in a .json file called 'my_saved_key.json' in the same folder as this notebook
# =============================================================================


# =============================================================================
# Prepare the environment
# =============================================================================
# To run this notebook you need to install the following dependencies.
# 
# # pip install pyjwt requests cryptography
# After installing them, you are ready to proceed
# =============================================================================

# Copernicus Land Monitoring Service (CLMS) Data Download
import json
import jwt
import time
import requests

base_url = 'https://land.copernicus.eu'

# =============================================================================
# Do the authentication steps
# =============================================================================
# First load the service key downloaded from the website and create a grant request

# Load saved key from filesystem
service_key = json.load(open("E:\Rachit\IGBP_Project\Data\SWE\CLMS_API_Token.json", 'rb'))

private_key = service_key['private_key'].encode('utf-8')

claim_set = {
    "iss": service_key['client_id'],
    "sub": service_key['user_id'],
    "aud": service_key['token_uri'],
    "iat": int(time.time()),
    "exp": int(time.time() + (60 * 60)),
}
grant = jwt.encode(claim_set, private_key, algorithm='RS256')

# With the grant, you can request the access token to use the API

result = requests.post(
        service_key["token_uri"],
        headers={
            "Accept": "application/json",
            "Content-Type": "application/x-www-form-urlencoded",
        },
        data={
            "grant_type": "urn:ietf:params:oauth:grant-type:jwt-bearer",
            "assertion": grant,
        },
)

access_token_info_json = result.json()
access_token = access_token_info_json.get('access_token')
print(access_token)

# =============================================================================
# Find the items that can be downloaded
# =============================================================================
# This step makes a search query and prints the title, UID and link of the first 25 items found.

requests.get('http://nohost/api/@search?portal_type=DataSet&metadata_fields=UID&metadata_fields=dataset_full_format&&metadata_fields=dataset_download_information', headers={'Accept': 'application/json'})


# =============================================================================
# Search prepackaged files that can be downloaded
# =============================================================================
# You can use the @search endpoint to search for datasets

url_downloadable_prepackaged = f"{base_url}/api/@search?portal_type=DataSet&metadata_fields=UID&metadata_fields=downloadable_files"
headers = {'Accept': 'application/json', 'Authorization': f'Bearer {access_token}'}
response_downloadable_prepackaged = requests.get(url_downloadable_prepackaged, headers=headers)

json_downloadable_prepackaged = response_downloadable_prepackaged.json()

for item in json_downloadable_prepackaged.get('items', []):
    print('Product information:')
    print('\n')
    print(3*' '+'Product title: "'+ item['title'].replace(':','')+'"')
    print(3*' '+'UID: "'+item['UID'].replace(' ','')+'"')
    print(3*' '+'Product link: '+ item['@id']+'')
    print('\n')
    
    for downloadable_file in item.get('downloadable_files', {}).get('items', []):
        print('Download option for "'+item['title'].replace(':','')+'":')
        display(downloadable_file)
        print(2*'\n')

    print('\n')
    print(100*'_')

# =============================================================================
# Download example dataset
# =============================================================================
# You are welcome to specify the configuration of your data request by using the data variable. 
If the request is successful, then a task ID will be created under your account and when the data extraction is ready to be downloaded you will receive an email with the download link.

url_download_item = f"{base_url}/api/@datarequest_post"
headers = {'Accept': 'application/json', 'Authorization': f'Bearer {access_token}'}
data = {'Datasets': [{'DatasetID': 'a8d945f0edd143a0a5240c28bafa23da', 'DatasetDownloadInformationID': 'f328dbcc-9069-4240-b8bb-6d5df918671a', 'OutputFormat': 'Geotiff', 'OutputGCS': 'EPSG:4326', 'NUTS': 'ITC11'}, {'DatasetID': 'a8d945f0edd143a0a5240c28bafa23da', 'DatasetDownloadInformationID': 'f328dbcc-9069-4240-b8bb-6d5df918671a', 'OutputFormat': 'Geotiff', 'OutputGCS': 'EPSG:4326', 'NUTS': 'ITC1'}]}
response_download_item = requests.post(url_download_item, headers=headers, json=data)
print('Response:',response_download_item)
print('Message with TaskID:\n'+ response_download_item.text)


# =============================================================================
# To calculate milliseconds for date
# from datetime import datetime, timedelta

# # Epoch start time
# epoch_start = datetime(1970, 1, 1)

# # Define the dates
# start_date = datetime(2006, 1, 1)
# end_date = datetime(2006, 3, 1)

# # Calculate the time difference in milliseconds from the epoch
# start_date_milliseconds = int((start_date - epoch_start).total_seconds() * 1000)
# end_date_milliseconds = int((end_date - epoch_start).total_seconds() * 1000)

# print("Start Date (Jan 1, 2006) in milliseconds since epoch:", start_date_milliseconds)
# print("End Date (March 1, 2006) in milliseconds since epoch:", end_date_milliseconds)
# =============================================================================

requests.post('http://nohost/api/@datarequest_post', 
              headers={'Accept': 'application/json', 'Content-Type': 'application/json', 'Authorization': 'Bearer <REDACTED>'}, 
              json={'Datasets': [{'DatasetID': 'a8d945f0edd143a0a5240c28bafa23da', 'DatasetDownloadInformationID': 'f328dbcc-9069-4240-b8bb-6d5df918671a', 'OutputFormat': 'Geotiff', 'OutputGCS': 'EPSG:4326', 'BoundingBox': [21.125, 90.125, 32.125, 72.125], 'TemporalFilter': {'StartDate': 1136073600000, 'EndDate': 1141171200000}}]})
