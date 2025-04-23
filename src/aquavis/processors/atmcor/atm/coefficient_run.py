import os
import calendar
import numpy as np
import pandas as pd
from datetime import datetime, timedelta

from src.aquavis.config.config import SatWaterConfig
from src.aquavis.processors.data_class import AquaVisDataLoader
from src.aquavis.processors.atmcor.atm.coefficient import MCDExtractWindow

class RunAtmCoefficients:

    def __init__(self, path_main: str, output_dir: str):

        self.loader = AquaVisDataLoader()
        self.params = self.loader.load_aquavis_data()
        self.paths = SatWaterConfig()._load_paths()

        self.path_main = path_main
        self.output_dir = output_dir

    def _read_atmospheric_coefficients(self) -> None:
        """Read atmospheric coefficients with fallback strategies."""
        date_str = f"{self.params.datetime.year}-{self.params.datetime.month}-{self.params.datetime.day}"

        # Try daily values first
        self._try_daily_coefficients(date_str)

        # Fallback to weekly averages if needed
        if self._needs_fallback():
            self._try_weekly_coefficients(date_str)

        # Final fallback to monthly averages
        if self._needs_fallback():
            self._try_monthly_coefficients()

        # Save the atmospheric parameters to a CSV file
        self._save_atmospheric_parameters()

        self.loader.save_aquavis_data(self.params)

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
        month_name = calendar.month_name[self.params.datetime.month]
        mcd_scanner = self._create_mcd_scanner("", "")

        self.params.aod = mcd_scanner.get_modis_monthly_mean(month_name, 'AOD', self.params.roi)
        self.params.water_vapour = mcd_scanner.get_modis_monthly_mean(month_name, 'Water_Vapor', self.params.roi)
        self.params.ozone = mcd_scanner.get_modis_monthly_mean(month_name, 'Total_Ozone', self.params.roi) / 1000

    def _create_mcd_scanner(self, start_date: str, end_date: str) -> MCDExtractWindow:
        """Create MCDExtractWindow instance with current configuration."""
        return MCDExtractWindow(
            dir_mod08=f"{self.params.networkdrive_letter}{self.paths['MOD08_D3_PATH']}",
            dir_mde=f"{self.params.networkdrive_letter}{self.paths['MDE_PATH']}",
            dir_temp=f"{self.params.networkdrive_letter}{self.paths['TEMP_COEF_PATH']}",
            ini_date=start_date,
            end_date=end_date,
            bounding_shp=self.params.roi
        )

    def _extract_coefficients(
            self,
            mcd_scanner: MCDExtractWindow,
            include_altitude: bool = True
    ) -> None:
        """Extract coefficients from MCD scanner."""
        mod08_data = mcd_scanner.run_extraction_mod08d3()

        self.params.aod = mod08_data['AOD_mean'].mean()
        self.params.water_vapour = mod08_data['WV_mean'].mean()
        self.params.ozone = mod08_data['OZ_mean'].mean() / 1000  # Convert to cm_atm

        if include_altitude:
            mde_data = mcd_scanner.run_extract_mde()
            self.params.altitude = mde_data['MDE_mean'].mean() / 1000  # Convert to km

    def _needs_fallback(self) -> bool:
        """Check if we need to try fallback coefficient sources."""
        return (np.isnan(self.params.aod) or self.params.aod == 0.0 or
                np.isnan(self.params.water_vapour) or self.params.water_vapour == 0.0 or
                np.isnan(self.params.ozone) or self.params.ozone == 0.0)

    def _save_atmospheric_parameters(self) -> None:
        """Save atmospheric parameters to CSV file."""
        params = {
            'img': [self.path_main],
            'aod': [self.params.aod],
            'wv': [self.params.water_vapour],
            'oz': [self.params.ozone],
            'alt': [self.params.altitude]
        }

        output_path = os.path.join(self.output_dir, 'atm_parameters.csv')
        pd.DataFrame(params).to_csv(output_path, index=False)