import os
import ast
import datetime
import pandas as pd
from typing import List

select_sat = 'sentinel'
period = ('20230101', '20231231')
main_folder = '/ddnlus/scratch/r3693/hls_water/sentinel_water'
output_dir = '/ddnlus/scratch/r3693/hls_water/HLS_DATASET/sentinel'
sentinel_tile_mapping = pd.read_csv('/ddnlus/r3693/hls_water/scripts/hls_water/src/satwater/auxfiles/tiles/sentinel_landsat_intersections.csv')

tiles = os.listdir(main_folder)

paths_table = pd.DataFrame(columns=['tile', 'path', 'output_path'])


def check_date_range(all_imgs: List[str], dates_period: List[str], select_sat: str = 'landsat') -> List[str]:
    """
    Filter images based on whether their dates fall within a given date range.

    Args:
        all_imgs (list): List of image file paths.
        dates_period (list): Start and end dates in the format ['YYYYMMDD', 'YYYYMMDD'].
        select_sat (str): Satellite type ('landsat' or 'sentinel').

    Returns:
        list: Filtered list of image file paths within the date range.
    """

    filtered_imgs = []
    start_date = datetime.datetime.strptime(dates_period[0], '%Y%m%d')
    end_date = datetime.datetime.strptime(dates_period[1], '%Y%m%d')

    for img in all_imgs:
        # Extract the date string from the image filename
        if select_sat == 'sentinel':
            date_str = os.path.basename(img).split('_')[2].split('T')[0]
        else:
            date_str = os.path.basename(img).split('_')[3]

        # Convert date string to datetime object
        img_date = datetime.datetime.strptime(date_str, '%Y%m%d')

        # Check if the date falls within the range
        if start_date <= img_date <= end_date:
            filtered_imgs.append(img)

    return filtered_imgs


for tile in tiles:

    if select_sat == 'landsat':

        # Match Sentinel tiles to Landsat tiles
        sent_tile = sentinel_tile_mapping[sentinel_tile_mapping['sentinel_tile'] == tile]
        landsat_tiles = ast.literal_eval(sent_tile["landsat_tiles"].values[0])
        path_rows = [f"{int(lt.split('_')[0]):03d}_{int(lt.split('_')[1]):03d}" for lt in landsat_tiles]

        # Collect image paths for Landsat tiles
        all_imgs_aux = []
        for path_row in path_rows:
            src_dir = os.path.join(main_folder, path_row, "L89")
            if not os.path.exists(src_dir):
                print(f"No images found for tile {path_row}. Skipping...")
                continue

            all_imgs_aux.extend([os.path.join(src_dir, img) for img in os.listdir(src_dir)])
    else:

        # Collect image paths for Sentinel tiles
        src_dir = os.path.join(main_folder, tile)
        all_imgs_aux = [
            os.path.join(src_dir, sub_dir, os.listdir(os.path.join(src_dir, sub_dir))[0])
            for sub_dir in os.listdir(src_dir)
        ]

    # Filter images by date range
    all_imgs = check_date_range(all_imgs_aux, period, select_sat)

    # Prepare output paths
    output_paths = [
        os.path.join(output_dir, "atmcor", select_sat, tile, os.path.splitext(os.path.basename(img))[0])
        for img in all_imgs
    ]

    # Append rows to paths_table
    for i in range(len(all_imgs)):
        paths_table = paths_table.append({'tile': tile, 'path': all_imgs[i], 'output_path': output_paths[i]},
                                         ignore_index=True)

    print(f"Found {len(all_imgs)} images for tile {tile}")

paths_table.to_csv('/ddnlus/r3693/hls_water/scripts/hls_water/src/satwater/auxfiles/tiles/paths_sentinel_toa.txt', index=False)