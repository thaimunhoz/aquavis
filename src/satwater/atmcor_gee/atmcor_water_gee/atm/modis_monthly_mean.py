import os
import json
import numpy as np
import rasterio

def find_missing_dois(base_dir, expected_dois):

    missing_dois_by_year = {}

    # Iterate through each year folder
    for year in sorted(os.listdir(base_dir)):

        year_path = os.path.join(base_dir, year)

        if os.path.isdir(year_path):

            # Get the list of DOIs in this year folder
            existing_dois = set(os.listdir(year_path))

            int_list = sorted([int(s) for s in existing_dois])

            # Find missing DOIs
            missing_dois = sorted(set(expected_dois) - set(int_list))

            if missing_dois:

                missing_dois_by_year[year] = sorted(missing_dois)

    return missing_dois_by_year



import os
import shutil
from pathlib import Path
from datetime import datetime

def doy_to_month(doy, year):

    """Convert Day of Year (DOY) to a month name."""

    date = datetime.strptime(f"{year}-{doy:03}", "%Y-%j")  # Convert DOY to date

    return date.strftime("%B")  # Get the full month name


def group_files_by_month(base_dir):

    years_dict = {}

    # Iterate over each year folder
    for year in sorted(os.listdir(base_dir)):

        year_path = os.path.join(base_dir, year)

        if os.path.isdir(year_path):

            monthly_path = {}

            # Iterate over each file in the year folder
            for file_name in os.listdir(year_path):

                # Extract the DOY from the file name (assuming format like 'DOY_001', modify if necessary)
                doy_str = file_name  # Extract DOY part
                doy = int(doy_str)  # Convert DOY to integer

                # Get the month for this DOY
                month_name = doy_to_month(doy, int(year))

                # Create target folder for the month if it doesn't exist
                #target_month_folder = os.path.join(target_dir, month_name)
                #Path(target_month_folder).mkdir(parents=True, exist_ok=True)

                # Move the file to the correct month folder
                src_file = os.path.join(year_path, file_name)
                #dst_file = os.path.join(target_month_folder, file_name)

                image_files = os.listdir(src_file + '/' + '061')

                for image in image_files:

                    if month_name in monthly_path:
                        monthly_path[month_name].append(src_file + '/' + '061' + '/' + image)
                    else:
                        monthly_path[month_name] = [src_file + '/' + '061' + '/' + image]

                    #with rasterio.open(src_file + '/' + '061' + '/' + image) as src:
                    #    out_image = src.read(1)
                    #    out_meta = src.meta.copy()

                    #    out_meta.update({
                    #        "driver": "GTiff",
                    #        "height": out_image.shape[0],
                    #        "width": out_image.shape[1],
                    #        "transform": src.transform
                    #    })

                        # Save the clipped image
                    #    with rasterio.open(target_month_folder + '/' + image, "w", **out_meta) as dest:
                    #        dest.write(out_image, 1)

            years_dict[year] = monthly_path

            print(f"Year {year} processed.")

    return years_dict

# Example usage
#base_directory = r"Z:\dbcenter\products\atm\modis\C61\MOD08_D3"

#years_dict = group_files_by_month(base_directory)

#with open(r'Z:\guser\tml\mypapers\hls_synthetic\mod08_monthly_path.json', 'w') as json_file:
#    json.dump(years_dict, json_file)

# *********************************************************************************************************
import rioxarray
import xarray

def monthly_mean_modis(json_path, atm_parameter):

    # 1. Organizing the paths according to the month
    with open(json_path, 'r') as file:
        data = json.load(file)

    merged_data = {}

    for year, months in data.items():

        for month, values in months.items():

            if month in merged_data:
                merged_data[month].extend(values)
            else:
                merged_data[month] = values.copy()

    # 2. Create a stack for each month and calculte the mean
    for month, paths in merged_data.items():

        images_per_month = []

        for path in paths:

            if atm_parameter in path and 'xml' not in path:

                xda = rioxarray.open_rasterio(path, chunks={"x": 2080, "y": 2080})

                images_per_month.append(xda)

        images_stack = xarray.concat(images_per_month, dim="time")

        monthly_mean = images_stack.mean(dim="time")

        monthly_mean.rio.to_raster(r'Z:\guser\tml\mypapers\hls_synthetic\modis_monthly_mean' + '/' + 'MOD08_' + atm_parameter + '_' + month + '_mean.tiff')

        #print(f"Image saved in {r'Z:\guser\tml\mypapers\hls_synthetic\modis_monthly_mean' + '/' + 'MOD08_' + atm_parameter + '_' + month + '_mean.tiff'}")