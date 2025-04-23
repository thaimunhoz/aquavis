import os
import glob
import numpy as np
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
    MTD_TL = '/MTD_TL.xml'
    BAND_ID = '_B'
    GRANULE_DIR = '/GRANULE'
    IMG_DATA_DIR = '/IMG_DATA'
    MTD_MSIL1C = '/MTD_MSIL1C.xml'
    BANDS_ORDER = ['B01', 'B02', 'B03', 'B04', 'B05', 'B06',
                   'B07', 'B08', 'B8A', 'B09', 'B10', 'B11', 'B12']


class Metadata_MSI_S2:
    """Processes Sentinel-2 MSI metadata and atmospheric parameters."""

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
        self.rescale_factor()
        self._determine_roi()
        self._process_datetime()
        self._process_geometry()

        # Save the processed parameters
        self.loader.save_aquavis_data(self.params)

    def _process_band_names(self) -> None:
        """Process and order band names."""
        self.granule_path = glob.glob(os.path.join(self.path_main + self.constants.GRANULE_DIR, '*L1C_*'))[0]
        self.s2path = f"{self.granule_path}{self.constants.IMG_DATA_DIR}"

        band_files = [f for f in os.listdir(self.s2path) if self.constants.BAND_ID in f]
        band_files.insert(8, band_files.pop(-1))  # Adjust B8A position

        self.params.bandname = sorted(
            band_files,
            key=lambda x: next((i for i, band in enumerate(self.constants.BANDS_ORDER)
                                if band in x), float('inf'))
        )

    def _load_metadata_file(self) -> None:
        """Load and parse the MTL metadata file."""
        self.params.dict_metadata = tool.xml_to_json(f"{self.granule_path}{self.constants.MTD_TL}")
        self.params.type = str(self.params.dict_metadata['n1:Level-1C_Tile_ID'][
                            'n1:General_Info']['TILE_ID']['#text'][0:3])

    def rescale_factor(self) -> None:
        """Calculate rescale factors for all bands."""
        msi_metadata = tool.xml_to_json(self.path_main + self.constants.MTD_MSIL1C)
        quant_value = float(msi_metadata['n1:Level-1C_User_Product'][
                                'n1:General_Info']['Product_Image_Characteristics'][
                                'QUANTIFICATION_VALUE']['#text'])

        offset_list = msi_metadata['n1:Level-1C_User_Product'][
            'n1:General_Info']['Product_Image_Characteristics'][
            'Radiometric_Offset_List']['RADIO_ADD_OFFSET']

        rescale_aux = {}
        for i in range(13):
            try:
                offset = float(offset_list[i]['#text'])
            except (IndexError, KeyError, TypeError):
                offset = 0.0

            rescale_aux[i] = {'qvalue': quant_value, 'offset': offset}

        self.params.rescale = rescale_aux

    def _determine_roi(self) -> None:
        """Determine region of interest, with fallback to full image."""
        self.params.roi = tool.return_water(self.s2path, self.params.rescale, self.params.tile)

        if len(self.params.roi) == 0:
            sample_band = glob.glob(os.path.join(self.path_main, self.params.bandname[0]))[0]
            self.params.roi = tool.return_bbox(sample_band)

    def _process_datetime(self) -> None:
        """Extract and process date/time from metadata."""
        time_data = self.params.dict_metadata["n1:Level-1C_Tile_ID"][
            'n1:General_Info']['SENSING_TIME']['#text']

        date_acquired = time_data[0:10]
        date = datetime.strptime(date_acquired, '%Y-%m-%d').timetuple()
        scene_time = time_data[11:-1].split(':')

        time_hh = (int(scene_time[0]) +
                   (float(scene_time[1]) / 60) +
                   (float(scene_time[2]) / 3600))

        self.params.datetime = DateTime(
            day=date.tm_mday,
            month=date.tm_mon,
            year=date.tm_year,
            time_hh=time_hh
        )

    def _process_geometry(self) -> None:
        """Process sun and view angles for all bands."""
        tile_angles = self.params.dict_metadata['n1:Level-1C_Tile_ID'][
            'n1:Geometric_Info']['Tile_Angles']

        sun_az = float(tile_angles['Mean_Sun_Angle']['AZIMUTH_ANGLE']['#text'])
        sun_zn = float(tile_angles['Mean_Sun_Angle']['ZENITH_ANGLE']['#text'])

        view_angles = tile_angles["Mean_Viewing_Incidence_Angle_List"][
            'Mean_Viewing_Incidence_Angle']

        geometry_aux = {}
        for i in range(13):
            geometry_aux[i] = {
                'solar_az': sun_az,
                'solar_zn': sun_zn,
                'view_az': float(view_angles[i]['AZIMUTH_ANGLE']['#text']),
                'view_zn': float(view_angles[i]['ZENITH_ANGLE']['#text'])
            }

        self.params.geometry = geometry_aux