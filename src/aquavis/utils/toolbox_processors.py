import os
import ast
import datetime
import pandas as pd
from typing import List, Dict, Tuple
from src.aquavis.config.config import SatWaterConfig

'''Function used in the processor module'''

def _process_images(correction, input_paths: List, output_paths: List) -> None:
    """Process images (sequential or parallel)."""

    # Sequential processing
    for img_path, out_path in zip(input_paths, output_paths):

        # check if the image_path exists
        if os.path.exists(img_path):
            continue

        correction.run(img_path, out_path)

    # Alternative parallel processing (commented out)
    # with Pool(processes=12) as pool:
    #     args = zip(input_paths, output_paths, [params["select_sat"]] * len(input_paths), params)
    #     pool.starmap(corrector, args)

def atmcor_input_files(select_sat: str, input_dir: str, s2_tile: str, period: Tuple):
    """Search for input files in the database based on MSI tile and date range."""
    all_images = collect_images(select_sat, input_dir, s2_tile)
    return filter_images_by_date(select_sat, period, all_images)

def atmcor_output_files(select_sat: str, input_dir: str, s2_tile: str, period: Tuple, output_dir: str):
    """Prepare output paths for the processed images."""
    filtered_images = atmcor_input_files(select_sat, input_dir, s2_tile, period)

    return [
        os.path.join(output_dir, os.path.splitext(os.path.basename(img))[0])
        for img in filtered_images
    ]

def input_files(base_path: str) -> List:
    """Search for input files from the atmospheric correction step."""
    return [os.path.join(base_path, scene) for scene in os.listdir(base_path)]

def output_files(base_path: str, scene_images: List[str]) -> List:
    """Create output directory for adjacency correction."""
    output_paths = [
        os.path.join(base_path, os.path.basename(img))
        for img in scene_images
    ]
    return output_paths

def collect_images(select_sat: str, input_dir: str, tile: str) -> List[str]:
    """Collect all images for a given tile based on satellite type."""
    if select_sat == 'landsat':
        path_rows = _get_landsat_path_rows(tile)
        return _collect_landsat_images(input_dir, path_rows)
    else:
        return _collect_sentinel_images(input_dir, tile)

def filter_images_by_date(select_sat: str, period: Tuple[str, str], image_paths: List[str]) -> List[str]:
    """Filter images based on date range specified in parameters."""
    if not image_paths:
        return []

    start_date, end_date = [
        datetime.datetime.strptime(d, '%Y%m%d')
        for d in period
    ]

    filtered_images = []

    for img_path in image_paths:
        img_date = _parse_date_from_filename(img_path, select_sat)
        if start_date <= img_date <= end_date: filtered_images.append(img_path)
    return filtered_images

def _get_landsat_path_rows(tile: str) -> List[str]:
    """Convert Sentinel tile to corresponding Landsat path/rows."""
    sentinel_tile_mapping = _load_tile_mapping()
    sent_tile = sentinel_tile_mapping[sentinel_tile_mapping['sentinel_tile'] == tile]
    landsat_tiles = ast.literal_eval(sent_tile["landsat_tiles"].values[0])
    return [
        f"{int(lt.split('_')[0]):03d}_{int(lt.split('_')[1]):03d}"
        for lt in landsat_tiles
    ]

def _load_tile_mapping() -> pd.DataFrame:
    """Load the Sentinel-Landsat tile mapping CSV."""
    config = SatWaterConfig()._load_paths()
    mapping_path = config["sentinel_landsat_intersection"]
    return pd.read_csv(mapping_path)

def _parse_date_from_filename(filename: str, satellite_type: str) -> datetime.datetime:
    """Extract and parse date from satellite image filename."""
    basename = os.path.basename(filename)
    if satellite_type == 'sentinel':
        date_str = basename.split('_')[2].split('T')[0]
    else:
        date_str = basename.split('_')[3]
    return datetime.datetime.strptime(date_str, '%Y%m%d')

def _collect_landsat_images(input_dir: str, path_rows: List[str], landsat_generation: str = "L89") -> List[str]:
    """Collect Landsat image paths for given path/rows."""
    all_images = []
    for path_row in path_rows:
        src_dir = os.path.join(input_dir, path_row, landsat_generation)

        all_images.extend([
            os.path.join(src_dir, img)
            for img in os.listdir(src_dir)
        ])

    return all_images

def _collect_sentinel_images(input_dir, tile: str) -> List[str]:
    """Collect Sentinel image paths for given tile."""
    input_dir = os.path.join(input_dir, tile)

    if not os.path.exists(input_dir):
        print(f"No images found for tile {tile}")
        return []

    return [
        os.path.join(input_dir, sub_dir, os.listdir(os.path.join(input_dir, sub_dir))[0])
        for sub_dir in os.listdir(input_dir)
    ]