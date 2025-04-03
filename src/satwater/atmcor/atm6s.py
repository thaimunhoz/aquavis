import os
import ast
import logging
import datetime
import pandas as pd
from multiprocessing import Pool
from typing import List, Dict, Optional

from src.satwater.atmcor.atmcor_water.run_atmcor_water import AtmosphericCorrector

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class AtmcorProcessor:

    """Handles satellite image processing including date filtering and atmospheric correction."""

    def __init__(self, select_sat: str, params: Dict):

        """
        Args:
            select_sat (str): Satellite type ('landsat' or 'sentinel')
            params (Dict): Dictionary of processing parameters
        """
        self.select_sat = select_sat
        self.params = params
        self.sentinel_tile_mapping = self._load_tile_mapping()
        self.corrector = AtmosphericCorrector()

    def _load_tile_mapping(self) -> pd.DataFrame:
        """Load the Sentinel-Landsat tile mapping CSV."""
        mapping_path = r'C:\Users\tml411\Documents\Python Scripts\hls_water\src\satwater\auxfiles\tiles\sentinel_landsat_intersections.csv'
        return pd.read_csv(mapping_path)

    @staticmethod
    def _parse_date_from_filename(filename: str, satellite_type: str) -> datetime.datetime:
        """
        Extract and parse date from satellite image filename.

        Args:
            filename (str): Image filename
            satellite_type (str): 'landsat' or 'sentinel'

        Returns:
            datetime.datetime: Parsed date from filename
        """
        basename = os.path.basename(filename)

        if satellite_type == 'sentinel':
            date_str = basename.split('_')[2].split('T')[0]
        else:
            date_str = basename.split('_')[3]

        return datetime.datetime.strptime(date_str, '%Y%m%d')

    def filter_images_by_date(self, image_paths: List[str]) -> List[str]:
        """
        Filter images based on date range specified in parameters.

        Args:
            image_paths (List[str]): List of image file paths

        Returns:
            List[str]: Filtered list of image paths within date range
        """
        if not image_paths:
            return []

        start_date, end_date = [
            datetime.datetime.strptime(d, '%Y%m%d')
            for d in self.params['aux_info']['period']
        ]

        filtered_images = []

        for img_path in image_paths:
            try:
                img_date = self._parse_date_from_filename(img_path, self.select_sat)
                if start_date <= img_date <= end_date:
                    filtered_images.append(img_path)
            except (IndexError, ValueError) as e:
                logger.warning(f"Could not parse date from {img_path}: {e}")

        return filtered_images

    def _get_landsat_path_rows(self, tile: str) -> List[str]:
        """Convert Sentinel tile to corresponding Landsat path/rows."""
        sent_tile = self.sentinel_tile_mapping[self.sentinel_tile_mapping['sentinel_tile'] == tile]
        landsat_tiles = ast.literal_eval(sent_tile["landsat_tiles"].values[0])
        return [
            f"{int(lt.split('_')[0]):03d}_{int(lt.split('_')[1]):03d}"
            for lt in landsat_tiles
        ]

    def _collect_landsat_images(self, path_rows: List[str]) -> List[str]:
        """Collect Landsat image paths for given path/rows."""
        all_images = []
        generation = self.params['landsat']['generation']
        input_dir = self.params['landsat']['input_dir']

        for path_row in path_rows:
            src_dir = os.path.join(input_dir, path_row, generation)

            if not os.path.exists(src_dir):
                logger.warning(f"No images found for tile {path_row}. Skipping...")
                continue

            all_images.extend([
                os.path.join(src_dir, img)
                for img in os.listdir(src_dir)
            ])

        return all_images

    def _collect_sentinel_images(self, tile: str) -> List[str]:
        """Collect Sentinel image paths for given tile."""
        input_dir = os.path.join(self.params['sentinel']['input_dir'], tile)

        if not os.path.exists(input_dir):
            logger.warning(f"No images found for tile {tile}")
            return []

        return [
            os.path.join(input_dir, sub_dir, os.listdir(os.path.join(input_dir, sub_dir))[0])
            for sub_dir in os.listdir(input_dir)
        ]

    def collect_images(self, tile: str) -> List[str]:
        """
        Collect all images for a given tile based on satellite type.

        Args:
            tile (str): Tile identifier

        Returns:
            List[str]: List of image paths
        """
        if self.select_sat == 'landsat':
            path_rows = self._get_landsat_path_rows(tile)
            return self._collect_landsat_images(path_rows)
        else:
            return self._collect_sentinel_images(tile)

    def process_tile(self, tile: str) -> None:
        """
        Process all images for a given tile through atmospheric correction.

        Args:
            tile (str): Tile identifier
        """
        # Collect and filter images
        all_images = self.collect_images(tile)
        filtered_images = self.filter_images_by_date(all_images)

        if not filtered_images:
            logger.info(f"No images found for tile {tile} in the specified date range.")
            return

        # Prepare output paths
        output_dir = os.path.join(self.params["output_dir"], "atmcor")
        os.makedirs(output_dir, exist_ok=True)

        output_paths = [
            os.path.join(output_dir, os.path.splitext(os.path.basename(img))[0])
            for img in filtered_images
        ]

        # Process images
        self._process_images(filtered_images, output_paths, tile)

    def _process_images(self, input_paths: List[str], output_paths: List[str], tile: str) -> None:
        """
        Process images through atmospheric correction (sequential or parallel).

        Args:
            input_paths (List[str]): Input image paths
            output_paths (List[str]): Corresponding output paths
            tile (str): Tile identifier
        """
        # Sequential processing
        for img_path, out_path in zip(input_paths, output_paths):
            if os.path.exists(out_path):
                logger.info(f"Skipping {out_path}, already exists.")
                continue

            try:

                self.corrector.run_correction(img_path, out_path, self.select_sat, tile)
                logger.info(f"Successfully processed {img_path}")
            except Exception as e:
                logger.error(f"Failed to process {img_path}: {e}")

        # Alternative parallel processing (commented out)
        # if self.params['aux_info'].get('n_cores', 1) > 1:
        #     self._process_images_parallel(input_paths, output_paths, tile)

    def _process_images_parallel(self, input_paths: List[str], output_paths: List[str], tile: str) -> None:
        """Process images in parallel using multiprocessing."""
        n_cores = self.params['aux_info'].get('n_cores', 1)

        with Pool(processes=n_cores) as pool:
            args = zip(input_paths, output_paths, [self.select_sat] * len(input_paths))
            results = pool.starmap(self.corrector.run_correction, args)
            logger.info(f"Atmospheric correction completed for tile {tile}. Results: {results}")


def run(select_sat: str, params: Dict) -> None:
    """
    Main function to execute atmospheric correction processing.

    Args:
        select_sat (str): Satellite type ('landsat' or 'sentinel')
        params (Dict): Dictionary of processing parameters
    """

    processor = AtmcorProcessor(select_sat, params)

    # Get tiles from parameters (default to empty list if not found)
    tiles = params[select_sat].get('tiles', [])

    for tile in tiles:
        processor.process_tile(tile)