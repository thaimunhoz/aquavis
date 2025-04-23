import os
import glob
import pandas as pd
import geopandas as gpd
import xml.etree.ElementTree as ET
from src.aquavis.config.config import SatWaterConfig

def get_sentinel_tile_info(tile):

    config = SatWaterConfig()
    paths = config._load_paths()

    sentinel_landsat_int_path = paths["sentinel_landsat_intersection"]
    sentinel_landsat_int_df = pd.read_csv(sentinel_landsat_int_path)

    return sentinel_landsat_int_df[sentinel_landsat_int_df['sentinel_tile'] == tile].iloc[0]

def get_tile_shp(sen_tile_target, epsg_code):

    config = SatWaterConfig()
    paths = config._load_paths()

    sentinel_tiles_path = paths["tiles_sentinel"]

    df_sentinel = gpd.read_file(sentinel_tiles_path, encoding="UTF-8")
    df_select = df_sentinel[df_sentinel['Name'] == sen_tile_target]
    return df_select.to_crs(f'EPSG:{epsg_code}')

def tile_epsg(tile):

    sentinel_tile_info = get_sentinel_tile_info(tile)
    tile_epsg = sentinel_tile_info['sentinel_epsg']
    return tile_epsg

def read_metadata(self):
    tree = ET.parse(os.path.join(self.path_main, 'MTD.xml'))
    root = tree.getroot()
    return xml_to_dict(root)

def xml_to_dict(element):

    if len(element) == 0:
        return element.text

    result = {}

    for child in element:

        child_dict = xml_to_dict(child)

        if child.tag in result:
            if isinstance(result[child.tag], list):
                result[child.tag].append(child_dict)
            else:
                result[child.tag] = [result[child.tag], child_dict]
        else:
            result[child.tag] = child_dict

    return result

def create_dir(path):

    if not os.path.exists(path):

        os.makedirs(path)

    return path

def _get_scene_name(select_sat, scene_path: str) -> str:
    """Get scene name based on satellite type."""
    if select_sat == "sentinel":
        return glob.glob(os.path.join(scene_path, '*.SAFE*'))[0]
    return glob.glob(os.path.join(scene_path, 'LC*'))[0]





