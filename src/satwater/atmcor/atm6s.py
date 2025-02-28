import os
import ast
import datetime
import logging
import pandas as pd
from multiprocessing import Pool
from typing import List, Tuple
from src.satwater.atmcor.atmcor_water.run_atmcor_water import run_gceratmos

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

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
        try:
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
        except (IndexError, ValueError) as e:
            logging.warning(f"Failed to process image {img}: {e}")

    return filtered_imgs

def run(select_sat: str, params: dict) -> None:

    """
    Execute atmospheric correction processing for satellite data.

    Args:
        select_sat (str): Selected satellite type ('landsat' or 'sentinel').
        params (dict): Dictionary of parameters including input/output directories, tile info, and processing settings.
    """

    tiles = [params[select_sat].get('tiles', [])]
    sentinel_tile_mapping = pd.read_csv(r'C:\Users\tml411\Documents\Python Scripts\hls_water\src\satwater\auxfiles\tiles\sentinel_landsat_intersections.csv')

    for tile in tiles:
        try:
            if select_sat == 'landsat':

                # Match Sentinel tiles to Landsat tiles
                sent_tile = sentinel_tile_mapping[sentinel_tile_mapping['sentinel_tile'] == tile]
                landsat_tiles = ast.literal_eval(sent_tile["landsat_tiles"].values[0])
                path_rows = [f"{int(lt.split('_')[0]):03d}_{int(lt.split('_')[1]):03d}" for lt in landsat_tiles]

                # Collect image paths for Landsat tiles
                all_imgs_aux = []
                for path_row in path_rows:
                    src_dir = os.path.join(params[select_sat]["input_dir"], path_row, params[select_sat]["generation"])
                    if not os.path.exists(src_dir):
                        logging.info(f"No images found for tile {path_row}. Skipping...")
                        continue

                    all_imgs_aux.extend([os.path.join(src_dir, img) for img in os.listdir(src_dir)])
            else:

                # Collect image paths for Sentinel tiles
                src_dir = os.path.join(params[select_sat]["input_dir"], tile)
                all_imgs_aux = [
                    os.path.join(src_dir, sub_dir, os.listdir(os.path.join(src_dir, sub_dir))[0])
                    for sub_dir in os.listdir(src_dir)
                ]

            # Filter images by date range
            all_imgs = check_date_range(all_imgs_aux, params['aux_info']['period'], select_sat)
            if not all_imgs:
                logging.info(f"No images found for tile {tile} in the specified date range.")
                continue

            # Prepare output paths
            output_paths = [
                os.path.join(params["output_dir"], "atmcor", select_sat, tile, os.path.splitext(os.path.basename(img))[0])
                for img in all_imgs
            ]

            # Run atmospheric correction in a for loop:
            # for img, output_path in zip(all_imgs, output_paths):
            #     try:
            #         run_gceratmos(img, output_path, select_sat)
            #     except Exception as e:
            #         print(f"An error occurred while processing image {img}: {e}")

                #logging.info(f"Atmospheric correction completed for tile {tile}.")

            # Run in a for loop
            for img, output_path in zip(all_imgs, output_paths):
                run_gceratmos(img, output_path, select_sat)

            # Run atmospheric correction in parallel
            # with Pool(processes=params['aux_info']['n_cores']) as pool:
            #
            #     args = zip(all_imgs, output_paths, [select_sat] * len(all_imgs))
            #     results = pool.starmap(run_gceratmos, args)
            #     logging.info(f"Atmospheric correction completed for tile {tile}. Results: {results}")

        except Exception as e:

            logging.error(f"An error occurred while processing tile {tile}: {e}")