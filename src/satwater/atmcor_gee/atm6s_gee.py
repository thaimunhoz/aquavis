import os
import ast
import datetime
import logging
import pandas as pd
from itertools import chain
from datetime import datetime
from multiprocessing import Pool
from typing import List, Tuple

from src.satwater.input_gee import toa_gee
from src.satwater.atmcor_gee.gceratmos_hls.run_gceratmos_gee import run_gceratmos_gee

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run(select_sat: str, params: dict) -> None:

    """
    Execute atmospheric correction processing for satellite data.

    Args:
        select_sat (str): Selected satellite type ('landsat' or 'sentinel').
        params (dict): Dictionary of parameters including input/output directories, tile info, and processing settings.
    """

    tiles = [params[select_sat].get('tiles', [])]
    sentinel_tile_mapping = pd.read_csv(r'C:\Users\tml411\Documents\Python Scripts\hls_water\src\satwater\auxfiles\tiles\sentinel_landsat_intersections.csv')

    # Collect image paths for Landsat tiles
    all_imgs = []

    for tile in tiles:

        try:
            if select_sat == 'landsat':

                # Match Sentinel tiles to Landsat tiles
                sent_tile = sentinel_tile_mapping[sentinel_tile_mapping['sentinel_tile'] == tile]
                landsat_tiles = ast.literal_eval(sent_tile["landsat_tiles"].values[0])
                path_rows = [f"{int(lt.split('_')[0]):03d}_{int(lt.split('_')[1]):03d}" for lt in landsat_tiles]

                for path_row in path_rows:

                    all_images_aux = []

                    start_date = datetime.strptime(params['aux_info']['period'][0], '%Y%m%d').strftime('%Y-%m-%d')
                    end_date = datetime.strptime(params['aux_info']['period'][1], '%Y%m%d').strftime('%Y-%m-%d')

                    all_images_aux.append(toa_gee.get_available_dates_landsat(path_row, start_date, end_date))

                    if not all_images_aux:
                        logging.info(f"No images found for tile {path_row} in the specified date range.")
                        continue
                    else:
                        # Dictonary with available path_row as keys the dates as list values
                        all_imgs.append(all_images_aux[0])

                all_imgs = list(chain.from_iterable(all_imgs))
                output_paths = [os.path.join(params["output_dir"], "atmcor", img) for img in all_imgs]

            else:

                all_images_aux = []

                start_date = datetime.strptime(params['aux_info']['period'][0], '%Y%m%d').strftime('%Y-%m-%d')
                end_date = datetime.strptime(params['aux_info']['period'][1], '%Y%m%d').strftime('%Y-%m-%d')

                all_images_aux.append(toa_gee.get_available_dates_sentinel(tile, start_date, end_date))

                if not all_images_aux:
                    logging.info(f"No images found for tile {path_row} in the specified date range.")
                    continue
                else:
                    # Dictonary with available path_row as keys the dates as list values
                    all_imgs.append(all_images_aux[0])

                all_imgs = list(chain.from_iterable(all_imgs))
                output_paths = [os.path.join(params["output_dir"], "atmcor", img) for img in all_imgs]

            # Run atmospheric correction in a for loop:
            for img, output_path in zip(all_imgs, output_paths):
                run_gceratmos_gee(img, output_path, select_sat)
                logging.info(f"Atmospheric correction completed for tile {tile}.")

            # Run atmospheric correction in parallel
            # with Pool(processes=params['aux_info']['n_cores']) as pool:
            #
            #     args = zip(all_imgs, output_paths, [select_sat] * len(all_imgs))
            #     results = pool.starmap(run_gceratmos_gee, args)
            #     logging.info(f"Atmospheric correction completed for tile {tile}. Results: {results}")

        except Exception as e:

            logging.error(f"An error occurred while processing tile {tile}: {e}")