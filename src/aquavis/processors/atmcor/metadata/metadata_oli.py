import os
import glob
import numpy as np
import rioxarray as rxr
from datetime import datetime
from dataclasses import dataclass

from src.aquavis.processors.atmcor import toolbox as tool
from src.aquavis.processors.data_class import AquaVisDataLoader

@dataclass
class DateTime:
    """Stores date and time values with type hints."""
    day: int = np.nan
    month: int = np.nan
    year: int = np.nan
    time_hh: float = np.nan

    def __str__(self) -> str:
        return f'day: {self.day}, month: {self.month}, year: {self.year}, time_hh: {self.time_hh}'

@dataclass
class Constants:
    MTL_FILE_PATTERN = '*.xml'
    BAND_ID = '_B'
    MTL_ID = 'MTL'
    ANGLE_FILE_PATTERNS = {
        'solar_azimuth': '*SAA.tif',
        'solar_zenith': '*SZA.tif',
        'view_azimuth': '*VAA.tif',
        'view_zenith': '*VZA.tif'
    }
    OPTICAL_BANDS_ORDER = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7']
    EXCLUDED_BANDS = ['B9', 'B10', 'B11', 'B12']

class Metadata_OLI_L89:
    """Processes Landsat 8/9 OLI metadata and atmospheric parameters."""

    def __init__(self, satellite: str, path_main: str):

        self.loader = AquaVisDataLoader()
        self.params = self.loader.load_aquavis_data()

        self.params.sensor = satellite
        self.path_main = path_main

        self.constants = Constants()

    def run(self) -> None:
        """Main method to process all metadata."""
        self._process_band_names()
        self._load_metadata_file()
        self._rescale_factors()
        self._determine_roi()
        self._process_datetime()
        self._process_geometry()

        # Save the processed parameters
        self.loader.save_aquavis_data(self.params)

    def _process_band_names(self) -> None:
        """Process and order optical band names, excluding thermal bands."""
        all_bands = [
            f for f in os.listdir(self.path_main)
            if self.constants.BAND_ID in f and
               not any(excluded in f for excluded in self.constants.EXCLUDED_BANDS)
        ]

        self.params.bandname = sorted(
            all_bands,
            key=lambda x: next((i for i, band in enumerate(self.constants.OPTICAL_BANDS_ORDER)
                                if band in x), float('inf'))
        )

    def _load_metadata_file(self) -> None:
        """Load and parse the MTL metadata file."""
        mtl_files = glob.glob(os.path.join(self.path_main, self.constants.MTL_FILE_PATTERN))
        mtl_file = next((f for f in mtl_files if self.constants.MTL_ID in f), None)

        if not mtl_file:
            raise FileNotFoundError(f"No MTL file found in {self.path_main}")

        self.params.dict_metadata = tool.xml_to_json(mtl_file)
        self.params.type = str(self.params.dict_metadata['LANDSAT_METADATA_FILE'][
                            'PRODUCT_CONTENTS']['LANDSAT_PRODUCT_ID'][0:4])

    def _rescale_factors(self) -> None:
        """Extract rescaling factors for all bands."""

        rescale_aux = {}
        for i in range(1, 9):  # Bands 1-8
            add_key = f'REFLECTANCE_ADD_BAND_{i}'
            mult_key = f'REFLECTANCE_MULT_BAND_{i}'

            rescale_aux[i - 1] = {
                'add': float(self.params.dict_metadata['LANDSAT_METADATA_FILE'][
                                 'LEVEL1_RADIOMETRIC_RESCALING'][add_key]),
                'mult': float(self.params.dict_metadata['LANDSAT_METADATA_FILE'][
                                  'LEVEL1_RADIOMETRIC_RESCALING'][mult_key])
            }

        self.params.rescale = rescale_aux

    def _determine_roi(self) -> None:
        """Determine region of interest, with fallback to full image."""
        self.params.roi = tool.return_water(self.path_main, self.params.rescale, self.params.tile)

        if len(self.params.roi) == 0:
            sample_band = glob.glob(os.path.join(self.path_main, self.params.bandname[0]))[0]
            self.params.roi = tool.return_bbox(sample_band)

    def _process_datetime(self) -> None:
        """Extract and process acquisition datetime from metadata."""
        img_attrs = self.params.dict_metadata['LANDSAT_METADATA_FILE']['IMAGE_ATTRIBUTES']

        date_acquired = img_attrs['DATE_ACQUIRED']
        time_parts = img_attrs['SCENE_CENTER_TIME'][0:16].split(':')

        time_hh = (int(time_parts[0]) +
                   (float(time_parts[1]) / 60) +
                   (float(time_parts[2]) / 3600))

        date = datetime.strptime(date_acquired, '%Y-%m-%d').timetuple()
        self.params.datetime = DateTime(
            day=date.tm_mday,
            month=date.tm_mon,
            year=date.tm_year,
            time_hh=time_hh
        )

    def _process_geometry(self) -> None:
        """Process sun and view angles for all bands."""
        angle_data = {}

        # Load all angle files
        for angle_type, pattern in self.constants.ANGLE_FILE_PATTERNS.items():
            file_path = glob.glob(os.path.join(self.path_main, pattern))[0]
            angle_data[angle_type] = rxr.open_rasterio(file_path).values.astype(float)
            angle_data[angle_type][angle_data[angle_type] == 0] = np.nan
            angle_data[angle_type][angle_data[angle_type] == -9999] = np.nan

        # Calculate mean angles (divided by 100 as per original code)
        mean_angles = {
            'solar_az': np.nanmean(angle_data['solar_azimuth']) / 100,
            'solar_zn': np.nanmean(angle_data['solar_zenith']) / 100,
            'view_az': np.nanmean(angle_data['view_azimuth']) / 100,
            'view_zn': np.nanmean(angle_data['view_zenith']) / 100
        }

        # Apply same angles to all bands (1-8)
        metadata_aux = {}
        for i in range(8):
            metadata_aux[i] = mean_angles.copy()

        self.params.geometry = metadata_aux