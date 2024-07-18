# -*- coding: utf-8 -*-
"""
Created on Wed May 22 19:16:11 2024

@author: rachit
"""
# =============================================================================
# Source: ChatGPT
# =============================================================================
import os
import glob

# Specify the directory containing the text files
directory = r'\\172.16.20.21\wrd_p_igbp\IGBP_PROJECT_2023-24\IGBP_Project_files\Data\Weather_data\IMD\temp\Yamuna_test'  # Replace with the path to your directory

# Use glob to find all text files in the directory
file_paths = glob.glob(os.path.join(directory, '*.txt'))
a = file_paths
for file_path in file_paths:
    # Open the file for reading
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Remove the trailing comma from the first line if it exists
    if lines:
        lines[0] = lines[0].rstrip(',\n') + '\n'

    # Write the modified content back to the file
    with open(file_path, 'w') as file:
        file.writelines(lines)

print("Comma removed from the first line of all files.")
