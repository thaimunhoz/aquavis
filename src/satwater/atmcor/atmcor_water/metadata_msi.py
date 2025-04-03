import os
import glob
import calendar
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Any

from src.satwater.atmcor.atmcor_water import toolbox as tool
from src.satwater.atmcor.atmcor_water.atm.coefficient import MCDExtractWindow

@dataclass
class DateTime:
    """Stores date and time values with type hints."""
    day: int = np.nan
    month: int = np.nan
    year: int = np.nan
    time_hh: float = np.nan

    def __str__(self) -> str:
        return f'day: {self.day}, month: {self.month}, year: {self.year}, time_hh: {self.time_hh}'


class Metadata_MSI_S2:
    """Processes Sentinel-2 MSI metadata and atmospheric parameters."""

    # Constants
    MTD_TL = '/MTD_TL.xml'
    BAND_ID = '_B'
    GRANULE_DIR = '/GRANULE'
    IMG_DATA_DIR = '/IMG_DATA'
    MTD_MSIL1C = '/MTD_MSIL1C.xml'
    MOD08_D3_PATH = ':/dbcenter/products/atm/modis/C61/MOD08_D3'
    MDE_PATH = ':/dbcenter/products/land/dem30m'
    TEMP_COEF_PATH = ':/public/temp_dir'

    # Band processing order
    BANDS_ORDER = ['B01', 'B02', 'B03', 'B04', 'B05', 'B06',
                   'B07', 'B08', 'B8A', 'B09', 'B10', 'B11', 'B12']

    def __init__(
            self,
            path_main: str,
            path_dest: str,
            networkdrive_letter: str,
            satellite: str,
            aero_type: str,
            mode: Optional[str] = None,
            msi_tile: Optional[str] = None
    ):
        """Initialize metadata processor."""
        self.path_main = path_main
        self.path_dest = path_dest
        self.networkdrive_letter = networkdrive_letter
        self.satellite = satellite
        self.aero_type = aero_type
        self.mode = mode
        self.msi_tile = msi_tile

        # Initialize attributes with proper types
        self.type: str = "nan"
        self.bandname: List[str] = ["nan"]
        self.s2path: str = "nan"
        self.aod: float = np.nan
        self.water_vapour: float = np.nan
        self.ozone: float = np.nan
        self.altitude: float = np.nan
        self.geometry: Dict[int, Dict[str, float]] = {}
        self.datetime: DateTime = DateTime()
        self.rescale: Dict[int, Dict[str, float]] = {}
        self.roi: Any = None
        self.dict_metadata: Dict[str, Any] = {}

    def run(self) -> None:
        """Main method to process all metadata."""
        self._process_image_files()
        self._process_metadata()
        self._save_atmospheric_parameters()

    def _process_image_files(self) -> None:
        """Process image files and bands."""
        granule_path = self._find_granule_path()
        self.s2path = f"{granule_path}{self.IMG_DATA_DIR}"
        self.rescale_factor()

        self.roi = tool.return_water(self.s2path, self.rescale, self.msi_tile)
        if not self.roi:
            sample_file = os.path.join(self.s2path, os.listdir(self.s2path)[3])
            self.roi = tool.return_bbox(sample_file)

        self._process_band_names()

    def _find_granule_path(self) -> str:
        """Find and return the granule path."""
        granule_paths = glob.glob(os.path.join(self.path_main + self.GRANULE_DIR, '*L1C_*'))
        if not granule_paths:
            raise FileNotFoundError("No granule directory found with 'L1C_' pattern")
        return granule_paths[0]

    def _process_band_names(self) -> None:
        """Process and order band names."""
        band_files = [f for f in os.listdir(self.s2path) if self.BAND_ID in f]
        band_files.insert(8, band_files.pop(-1))  # Adjust B8A position

        # Sort bands according to predefined order
        self.bandname = sorted(
            band_files,
            key=lambda x: next((i for i, band in enumerate(self.BANDS_ORDER)
                                if band in x), float('inf'))
        )

    def _process_metadata(self) -> None:
        """Process all metadata components."""
        granule_path = self._find_granule_path()
        self.dict_metadata = tool.xml_to_json(f"{granule_path}{self.MTD_TL}")

        self.type = str(self.dict_metadata['n1:Level-1C_Tile_ID'][
                            'n1:General_Info']['TILE_ID']['#text'][0:3])

        self.date_and_time()
        self._process_geometry()
        self._read_atmospheric_coefficients()

    def date_and_time(self) -> None:
        """Extract and process date/time from metadata."""
        time_data = self.dict_metadata["n1:Level-1C_Tile_ID"][
            'n1:General_Info']['SENSING_TIME']['#text']

        date_acquired = time_data[0:10]
        date = datetime.strptime(date_acquired, '%Y-%m-%d').timetuple()
        scene_time = time_data[11:-1].split(':')

        time_hh = (int(scene_time[0]) +
                   (float(scene_time[1]) / 60) +
                   (float(scene_time[2]) / 3600))

        self.datetime = DateTime(
            day=date.tm_mday,
            month=date.tm_mon,
            year=date.tm_year,
            time_hh=time_hh
        )

    def _process_geometry(self) -> None:
        """Process sun and view angles for all bands."""
        tile_angles = self.dict_metadata['n1:Level-1C_Tile_ID'][
            'n1:Geometric_Info']['Tile_Angles']

        sun_az = float(tile_angles['Mean_Sun_Angle']['AZIMUTH_ANGLE']['#text'])
        sun_zn = float(tile_angles['Mean_Sun_Angle']['ZENITH_ANGLE']['#text'])

        view_angles = tile_angles["Mean_Viewing_Incidence_Angle_List"][
            'Mean_Viewing_Incidence_Angle']

        for i in range(13):
            self.geometry[i] = {
                'solar_az': sun_az,
                'solar_zn': sun_zn,
                'view_az': float(view_angles[i]['AZIMUTH_ANGLE']['#text']),
                'view_zn': float(view_angles[i]['ZENITH_ANGLE']['#text'])
            }

    def rescale_factor(self) -> None:
        """Calculate rescale factors for all bands."""
        msi_metadata = tool.xml_to_json(self.path_main + self.MTD_MSIL1C)
        quant_value = float(msi_metadata['n1:Level-1C_User_Product'][
                                'n1:General_Info']['Product_Image_Characteristics'][
                                'QUANTIFICATION_VALUE']['#text'])

        offset_list = msi_metadata['n1:Level-1C_User_Product'][
            'n1:General_Info']['Product_Image_Characteristics'][
            'Radiometric_Offset_List']['RADIO_ADD_OFFSET']

        for i in range(13):
            try:
                offset = float(offset_list[i]['#text'])
            except (IndexError, KeyError, TypeError):
                offset = 0.0

            self.rescale[i] = {'qvalue': quant_value, 'offset': offset}

    def _read_atmospheric_coefficients(self) -> None:
        """Read atmospheric coefficients with fallback strategies."""
        date_str = f"{self.datetime.year}-{self.datetime.month}-{self.datetime.day}"

        # Try daily values first
        self._try_daily_coefficients(date_str)

        # Fallback to weekly averages if needed
        if self._needs_fallback():
            self._try_weekly_coefficients(date_str)

        # Final fallback to monthly averages
        if self._needs_fallback():
            self._try_monthly_coefficients()

    def _try_daily_coefficients(self, date_str: str) -> None:
        """Attempt to get daily atmospheric coefficients."""
        mcd_scanner = self._create_mcd_scanner(date_str, date_str)
        self._extract_coefficients(mcd_scanner)

    def _try_weekly_coefficients(self, date_str: str) -> None:
        """Attempt to get weekly atmospheric coefficients."""
        date = datetime.strptime(date_str, '%Y-%m-%d')
        start_of_week = date - timedelta(days=date.weekday())
        end_of_week = start_of_week + timedelta(days=6)

        mcd_scanner = self._create_mcd_scanner(
            start_of_week.strftime('%Y-%m-%d'),
            end_of_week.strftime('%Y-%m-%d')
        )
        self._extract_coefficients(mcd_scanner, include_altitude=False)

    def _try_monthly_coefficients(self) -> None:
        """Attempt to get monthly atmospheric coefficients."""
        month_name = calendar.month_name[self.datetime.month]
        mcd_scanner = self._create_mcd_scanner("", "")

        self.aod = mcd_scanner.get_modis_monthly_mean(month_name, 'AOD', self.roi)
        self.water_vapour = mcd_scanner.get_modis_monthly_mean(month_name, 'Water_Vapor', self.roi)
        self.ozone = mcd_scanner.get_modis_monthly_mean(month_name, 'Total_Ozone', self.roi) / 1000

    def _create_mcd_scanner(self, start_date: str, end_date: str) -> MCDExtractWindow:
        """Create MCDExtractWindow instance with current configuration."""
        return MCDExtractWindow(
            dir_mod08=f"{self.networkdrive_letter}{self.MOD08_D3_PATH}",
            dir_mde=f"{self.networkdrive_letter}{self.MDE_PATH}",
            dir_temp=f"{self.networkdrive_letter}{self.TEMP_COEF_PATH}",
            ini_date=start_date,
            end_date=end_date,
            bounding_shp=self.roi
        )

    def _extract_coefficients(
            self,
            mcd_scanner: MCDExtractWindow,
            include_altitude: bool = True
    ) -> None:
        """Extract coefficients from MCD scanner."""
        mod08_data = mcd_scanner.run_extraction_mod08d3()

        self.aod = mod08_data['AOD_mean'].mean()
        self.water_vapour = mod08_data['WV_mean'].mean()
        self.ozone = mod08_data['OZ_mean'].mean() / 1000  # Convert to cm_atm

        if include_altitude:
            mde_data = mcd_scanner.run_extract_mde()
            self.altitude = mde_data['MDE_mean'].mean() / 1000  # Convert to km

    def _needs_fallback(self) -> bool:
        """Check if we need to try fallback coefficient sources."""
        return (np.isnan(self.aod) or self.aod == 0.0 or
                np.isnan(self.water_vapour) or self.water_vapour == 0.0 or
                np.isnan(self.ozone) or self.ozone == 0.0)

    def _save_atmospheric_parameters(self) -> None:
        """Save atmospheric parameters to CSV file."""
        params = {
            'img': [self.path_main],
            'aod': [self.aod],
            'wv': [self.water_vapour],
            'oz': [self.ozone],
            'alt': [self.altitude]
        }

        output_path = os.path.join(self.path_dest, 'atm_parameters.csv')
        pd.DataFrame(params).to_csv(output_path, index=False)