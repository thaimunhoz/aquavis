import os
import glob
import shutil
import rasterio
import numpy as np
import pandas as pd
from typing import List
import geopandas as gpd
from pathlib import Path
from rasterio.mask import mask
from datetime import datetime, timedelta

class MCDExtractWindow:

    def __init__(
            self,
            dir_mod08: str,
            dir_mde: str,
            dir_temp: str,
            ini_date: str,
            end_date: str,
            bounding_shp: str):

        self.dir_mod08 = Path(dir_mod08)
        self.dir_mde = Path(dir_mde)
        self.dir_temp = os.path.join(dir_temp, self._generate_random_string(10))
        self.ini_date = datetime.strptime(ini_date, "%Y-%m-%d")
        self.end_date = datetime.strptime(end_date, "%Y-%m-%d")
        self.bounding_shp = bounding_shp

    def run_extraction_mod08d3(self):

        all_scenes = self._get_available_images(self.dir_mod08)  # dobson units
        all_scenes = self._check_for_dates(all_scenes)

        dataset_info = []

        for scene in all_scenes:

            dict_values = {}
            path_oz = glob.glob(os.path.join(scene, '*Total_Ozone_Mean.tif'))
            path_aod = glob.glob(os.path.join(scene, '*AOD_550_Dark_Target_Deep_Blue_Combined_Mean.tif'))
            path_wv = glob.glob(os.path.join(scene, '*Atmospheric_Water_Vapor_Mean.tif'))

            if path_oz and path_aod and path_wv:

                dict_values["filename"] = os.path.basename(path_oz[0])
                date_string = os.path.basename(path_oz[0]).split('.')[1][1:]
                date_object = datetime.strptime(date_string, '%Y%j')
                formatted_date = date_object.strftime('%Y-%m-%d')
                dict_values["date"] = formatted_date
                value_mean, value_std, valid_pixels = self._get_avg_std_raster(path_oz[0])

                if value_mean == -9999:  # no valid data on my intersection
                    continue

                dict_values["OZ_mean"] = value_mean
                dict_values["OZ_std"] = value_std
                dict_values["Valid_pixels"] = valid_pixels

                path_cwv = path_aod[0]
                value_mean, value_std, _ = self._get_avg_std_raster(path_cwv)
                dict_values["AOD_mean"] = value_mean
                dict_values["AOD_std"] = value_std

                path_cwv = path_wv[0]
                value_mean, value_std, _ = self._get_avg_std_raster(path_cwv)
                dict_values["WV_mean"] = value_mean
                dict_values["WV_std"] = value_std

                image_dataset = pd.DataFrame([dict_values])
                dataset_info.append(image_dataset)

        if dataset_info:

            dataset_info = pd.concat(dataset_info, axis=0).reset_index(drop=True)

        return dataset_info

    def run_extract_mde(self):

        os.makedirs(self.dir_temp, exist_ok=True)
        path_vrt_temp = os.path.join(self.dir_mde, 'global_mde.vrt')

        value_mean, value_std, valid_pixels = self._get_avg_std_raster(path_vrt_temp)

        dict_values = {}
        dict_values["MDE_mean"] = value_mean
        dict_values["MDE_std"] = value_std
        dataset_info = pd.DataFrame([dict_values])

        shutil.rmtree(self.dir_temp)

        return dataset_info

    def _get_date_time(self, path_raster: str):

        date_string = os.path.basename(path_raster).split('.')[-4]
        year = int(date_string[:4])
        day_of_year = int(date_string[4:7])
        hour = int(date_string[7:9])
        minutes = int(date_string[9:11])
        dt = datetime(year=year, month=1, day=1) + timedelta(days=day_of_year - 1, hours=hour, minutes=minutes)

        return dt

    def _get_avg_std_raster(self, path_raster: str):

        # Open the GeoTIFF file
        with rasterio.open(path_raster) as src:

            gdf = self.bounding_shp
            gdf = gdf.to_crs(src.crs)

            raster_bounds = src.bounds

            # Iterate over each polygon in the shapefile
            for index, row in gdf.iterrows():
                geom = gpd.GeoDataFrame({'geometry': [row['geometry']]}).geometry

                # Get the bounding box of the polygon
                xmin, ymin, xmax, ymax = geom[0].bounds

                # check if the bounding box of the polygon intersects with the raster bounds
                if (xmin > raster_bounds[2] or
                        xmax < raster_bounds[0] or
                        ymin > raster_bounds[3] or
                        ymax < raster_bounds[1]):
                    return -9999, -9999, -9999

                else:  # yes, there is an intersection
                    from rasterio.mask import mask
                    # Read the subset of the raster data within the polygon area
                    subset, _ = mask(src, geom, nodata=-9999, crop=True, indexes=1)
                    array_1d = subset.flatten()
                    array_1d = array_1d[array_1d != -9999]

                    if array_1d.size == 0:
                        mean_i = np.nan
                        std_i = np.nan
                        valid_pixels = np.nan
                    else:
                        mean_i = np.mean(array_1d)
                        std_i = np.std(array_1d)
                        valid_pixels = array_1d.size / len(subset.flatten()) * 100

        return mean_i, std_i, valid_pixels

    def _get_available_images(self, dir_product: Path) -> List[str]:

        """Returns a list of all available images in the folder."""

        all_scenes = []
        for year_i in os.listdir(dir_product):
            for date in os.listdir(dir_product / year_i):
                for tile in os.listdir(dir_product / year_i / date):
                    path_tile = dir_product / year_i / date / tile
                    all_scenes.append(str(path_tile))
        return all_scenes

    def _check_for_dates(self, all_scenes: List[str]) -> List[str]:

        if len(all_scenes) == 0:
            return list()
        dates = []
        for scene in all_scenes:
            split_path = scene.split('\\')
            year_day_str = split_path[-3] + split_path[-2]
            dates.append(datetime.strptime(year_day_str, '%Y%j'))
        # Filter by date
        for i in range(len(dates)):
            if dates[i] < self.ini_date or dates[i] > self.end_date:
                all_scenes[i] = None
        all_scenes = [p for p in all_scenes if p is not None]
        return all_scenes

    def _generate_random_string(self, length_i):

        import random, string
        letters = string.ascii_letters + string.digits
        random_string = ''.join(random.choice(letters) for _ in range(length_i))
        return random_string

    def get_modis_monthly_mean(self, month, atm_parameter, geometry_gdf):

        median_modis_path = r'Z:\guser\tml\mypapers\hls_synthetic\modis_monthly_mean'
        files = os.listdir(median_modis_path)

        filtered_files = [f for f in files if month in f and atm_parameter in f]

        with rasterio.open(median_modis_path + '/' + filtered_files[0]) as src:

            out_meta = src.meta.copy()

            gdf_reprojected = gpd.GeoDataFrame(geometry=geometry_gdf.geometry.to_crs(out_meta['crs']))
            geometries = [feature['geometry'] for feature in gdf_reprojected.iterfeatures()]

            out_image, out_transform = mask(src, geometries, crop=True)

            parameter_value = out_image.mean()

        return parameter_value