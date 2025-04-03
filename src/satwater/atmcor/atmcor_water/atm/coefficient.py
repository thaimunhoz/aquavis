import string
import random
import rasterio
import numpy as np
import pandas as pd
import geopandas as gpd
from pathlib import Path
from rasterio.mask import mask
from rasterio.io import DatasetReader
from datetime import datetime, timedelta
from typing import List, Dict, Tuple, Optional, Union

class MCDExtractWindow:

    """Handles extraction of MODIS atmospheric data for specified time windows and regions."""

    # Constants for file patterns
    OZONE_FILE_PATTERN = '*Total_Ozone_Mean.tif'
    AOD_FILE_PATTERN = '*AOD_550_Dark_Target_Deep_Blue_Combined_Mean.tif'
    WV_FILE_PATTERN = '*Atmospheric_Water_Vapor_Mean.tif'
    MDE_FILE = 'global_mde.vrt'
    TEMP_DIR_LENGTH = 10

    def __init__(
            self,
            dir_mod08: Union[str, Path],
            dir_mde: Union[str, Path],
            dir_temp: Union[str, Path],
            ini_date: str,
            end_date: str,
            bounding_shp: gpd.GeoDataFrame
    ):
        """
        Initialize MODIS data extractor.

        Args:
            dir_mod08: Path to MOD08_D3 data directory
            dir_mde: Path to DEM data directory
            dir_temp: Base path for temporary directory
            ini_date: Start date in 'YYYY-MM-DD' format
            end_date: End date in 'YYYY-MM-DD' format
            bounding_shp: GeoDataFrame defining region of interest
        """
        self.dir_mod08 = Path(dir_mod08)
        self.dir_mde = Path(dir_mde)
        self.dir_temp = Path(dir_temp) / self._generate_random_string(self.TEMP_DIR_LENGTH)
        self.ini_date = datetime.strptime(ini_date, "%Y-%m-%d")
        self.end_date = datetime.strptime(end_date, "%Y-%m-%d")
        self.bounding_shp = bounding_shp

    def run_extraction_mod08d3(self) -> pd.DataFrame:

        """
                Extract MOD08_D3 atmospheric parameters for the specified date range and region.

                Returns:
                    DataFrame containing ozone, AOD, and water vapor statistics
                """

        all_scenes = self._get_available_modis_scenes()
        valid_scenes = self._filter_scenes_by_date(all_scenes)

        dataset_info = []

        for scene in valid_scenes:
            scene_data = self._process_modis_scene(scene)
            if scene_data is not None:
                dataset_info.append(scene_data)

        return pd.concat(dataset_info, ignore_index=True) if dataset_info else pd.DataFrame()

    def run_extract_mde(self) -> pd.DataFrame:
        """
                Extract DEM statistics for the specified region.

                Returns:
                    DataFrame containing DEM statistics
                """
        self.dir_temp.mkdir(parents=True, exist_ok=True)
        dem_path = self.dir_mde / self.MDE_FILE

        try:
            mean_val, std_val, _ = self._calculate_raster_stats(dem_path)
            result = pd.DataFrame({
                "MDE_mean": [mean_val],
                "MDE_std": [std_val]
            })
        finally:
            self._cleanup_temp_dir()

        return result

    def _generate_random_string(self, length_i):

        import random, string
        letters = string.ascii_letters + string.digits
        random_string = ''.join(random.choice(letters) for _ in range(length_i))
        return random_string

    def get_modis_monthly_mean(self, month: str, parameter: str, geometry: gpd.GeoDataFrame) -> float:

        """
                Get monthly mean value for a specific atmospheric parameter.

                Args:
                    month: Month name (e.g., 'January')
                    parameter: Atmospheric parameter ('AOD', 'Water_Vapor', 'Total_Ozone')
                    geometry: GeoDataFrame defining region of interest

                Returns:
                    Mean parameter value for the specified month and region
        """
        base_path = Path(r'Z:\guser\tml\mypapers\hls_synthetic\modis_monthly_mean')
        matching_files = [
            f for f in base_path.iterdir()
            if month in f.name and parameter in f.name
        ]

        if not matching_files:
            raise FileNotFoundError(f"No files found for {month} and {parameter}")

        with rasterio.open(matching_files[0]) as src:
            stats = self._calculate_raster_stats(src, geometry)
            return stats[0]  # Return mean value

    def _process_modis_scene(self, scene_path: Path) -> Optional[pd.DataFrame]:
        """Process a single MODIS scene and extract atmospheric parameters."""
        ozone_file = next(scene_path.glob(self.OZONE_FILE_PATTERN), None)
        aod_file = next(scene_path.glob(self.AOD_FILE_PATTERN), None)
        wv_file = next(scene_path.glob(self.WV_FILE_PATTERN), None)

        if not all([ozone_file, aod_file, wv_file]):
            return None

        scene_data = self._extract_scene_metadata(ozone_file)
        if scene_data is None:
            return None

        # Process ozone data
        oz_mean, oz_std, valid_pixels = self._calculate_raster_stats(ozone_file)
        if oz_mean == -9999:
            return None

        scene_data.update({
            "OZ_mean": oz_mean,
            "OZ_std": oz_std,
            "Valid_pixels": valid_pixels
        })

        # Process AOD data
        aod_mean, aod_std, _ = self._calculate_raster_stats(aod_file)
        scene_data.update({
            "AOD_mean": aod_mean,
            "AOD_std": aod_std
        })

        # Process water vapor data
        wv_mean, wv_std, _ = self._calculate_raster_stats(wv_file)
        scene_data.update({
            "WV_mean": wv_mean,
            "WV_std": wv_std
        })

        return pd.DataFrame([scene_data])

    def _extract_scene_metadata(self, ozone_file: Path) -> Optional[Dict]:
        """Extract metadata from ozone filename."""
        try:
            filename = ozone_file.name
            date_string = filename.split('.')[1][1:]
            date_obj = datetime.strptime(date_string, '%Y%j')

            return {
                "filename": filename,
                "date": date_obj.strftime('%Y-%m-%d')
            }
        except (IndexError, ValueError):
            return None

    def _calculate_raster_stats(
            self,
            raster_path: Union[Path, DatasetReader],
            geometry: Optional[gpd.GeoDataFrame] = None
    ) -> Tuple[float, float, float]:
        """
        Calculate statistics for raster data within specified geometry.

        Args:
            raster_path: Path to raster file or opened raster DatasetReader
            geometry: Optional GeoDataFrame (uses class bounding_shp if None)

        Returns:
            Tuple of (mean, std_dev, valid_pixels_percentage)
        """
        gdf = geometry if geometry is not None else self.bounding_shp

        with rasterio.open(raster_path) if isinstance(raster_path, Path) else nullcontext(raster_path) as src:
            gdf = gdf.to_crs(src.crs)
            raster_bounds = src.bounds

            for _, row in gdf.iterrows():
                geom = gpd.GeoDataFrame({'geometry': [row.geometry]}, crs=gdf.crs).geometry
                xmin, ymin, xmax, ymax = geom.bounds[0]

                if not (xmin > raster_bounds[2] or xmax < raster_bounds[0] or
                        ymin > raster_bounds[3] or ymax < raster_bounds[1]):
                    try:
                        subset, _ = mask(
                            src,
                            geom,
                            nodata=-9999,
                            crop=True,
                            indexes=1,
                            all_touched=True
                        )
                        flat_data = subset.flatten()
                        valid_data = flat_data[flat_data != -9999]

                        if valid_data.size == 0:
                            return np.nan, np.nan, np.nan

                        return (
                            np.mean(valid_data),
                            np.std(valid_data),
                            valid_data.size / flat_data.size * 100
                        )
                    except Exception as e:
                        print(f"Error processing raster: {e}")

        return -9999, -9999, -9999

    def _get_available_modis_scenes(self) -> List[Path]:
        """Get all available MODIS scenes in the data directory."""
        scenes = []
        for year_dir in self.dir_mod08.iterdir():
            if not year_dir.is_dir():
                continue
            for date_dir in year_dir.iterdir():
                if not date_dir.is_dir():
                    continue
                for tile_dir in date_dir.iterdir():
                    if tile_dir.is_dir():
                        scenes.append(tile_dir)
        return scenes

    def _filter_scenes_by_date(self, scenes: List[Path]) -> List[Path]:
        """Filter scenes by date range."""
        valid_scenes = []
        for scene in scenes:
            try:
                parts = scene.parts
                year_day = parts[-3] + parts[-2]
                scene_date = datetime.strptime(year_day, '%Y%j')
                if self.ini_date <= scene_date <= self.end_date:
                    valid_scenes.append(scene)
            except (IndexError, ValueError):
                continue
        return valid_scenes

    def _cleanup_temp_dir(self) -> None:
        """Clean up temporary directory if it exists."""
        try:
            if self.dir_temp.exists():
                import shutil
                shutil.rmtree(self.dir_temp)
        except Exception as e:
            print(f"Error cleaning up temp directory: {e}")

    @staticmethod
    def _generate_random_string(length: int) -> str:
        """Generate a random alphanumeric string of specified length."""
        chars = string.ascii_letters + string.digits
        return ''.join(random.choice(chars) for _ in range(length))


class nullcontext:
    """Context manager that does nothing, used for already opened raster files."""

    def __init__(self, enter_result=None):
        self.enter_result = enter_result

    def __enter__(self):
        return self.enter_result

    def __exit__(self, *excinfo):
        pass