import warnings
import numpy as np
from typing import Dict, Any
from dataclasses import dataclass
from math import pi, cos, log, exp

class Correction:
    """
      Performs atmospheric correction on satellite imagery using radiative transfer calculations.

      Implements corrections based on Gordon and Wang (1994) and Gordon (1997).
      """

    def __init__(self, parameters, atmosphere, arr: np.ndarray, index: int):
        """
        Initialize the atmospheric corrector.

        Args:
            parameters: Atmospheric parameters container
            atmosphere: Atmospheric transmission values container
            arr: Input array of TOA (Top of Atmosphere) values
            index: Band index for processing
        """
        self.params = parameters
        self.atmos = atmosphere
        self.arr = arr
        self.index = index
        self.arr_sr = np.full_like(arr, np.inf)  # Surface reflectance result

    def run(self) -> None:
        """Execute the atmospheric correction process."""
        # Calculate total gas transmissions
        Tg_OG = self._calculate_total_gas_transmission()
        Tg_H20 = self._get_atmospheric_value('Tg_H20')
        Tg_O3 = self._get_atmospheric_value('Tg_O3')

        # Calculate optical depths and diffuse transmittance
        optical_depth_O3 = -log(Tg_O3)
        optical_depth_Rayleigh = self._get_atmospheric_value('optical_depth_rayleigh')
        tdif_upward = self._calculate_diffuse_transmittance(optical_depth_Rayleigh, optical_depth_O3)

        # Get atmospheric intrinsic reflectance
        r_atm = self._get_atmospheric_value('p_atm')

        # Convert to TOA reflectance
        arr_toa = self._convert_to_toa_reflectance()

        # Apply atmospheric correction
        self.arr_sr = ((arr_toa / (Tg_OG * Tg_O3 * Tg_H20)) - r_atm) / tdif_upward

    def _calculate_total_gas_transmission(self) -> float:
        """Calculate combined transmission of various atmospheric gases."""
        gases = ['tg_OG_co', 'tg_OG_c02', 'tg_OG_o2', 'tg_OG_no2', 'tg_OG_ch4']
        transmission = 1.0

        for gas in gases:
            transmission *= self._get_atmospheric_value(gas)

        return transmission

    def _get_atmospheric_value(self, key: str) -> float:
        """Safely get atmospheric value with type conversion."""
        return float(self.atmos.values[self.index][key])

    def _calculate_diffuse_transmittance(self, rayleigh_depth: float, ozone_depth: float) -> float:
        """Calculate upward diffuse transmittance."""
        view_zenith_rad = np.radians(self._get_atmospheric_value('view_zn'))
        return exp(-((rayleigh_depth / 2) + ozone_depth) / cos(view_zenith_rad))

    def _convert_to_toa_reflectance(self) -> np.ndarray:
        """
        Convert input array to Top of Atmosphere reflectance based on satellite type.

        Returns:
            Array of TOA reflectance values
        """
        satellite_converters = {
            'OLCI_S3': self._convert_olci,
            'MSI_S2': self._convert_msi,
            'OLI_L8/9': self._convert_oli,
            'OCI_PACE': self._convert_oci
        }

        converter = satellite_converters.get(self.params.satellite)
        if converter is None:
            warnings.warn(f"Unsupported sensor type: {self.params.satellite}", UserWarning)
            return self.arr

        return converter()

    def _convert_olci(self) -> np.ndarray:
        """Convert OLCI_S3 data to TOA reflectance."""
        solar_zenith_rad = np.radians(self.params.geometry[self.index]['solar_zn'])
        return (self.arr * pi) / (float(self.params.rescale[self.index]) * cos(solar_zenith_rad))

    def _convert_msi(self) -> np.ndarray:
        """Convert MSI_S2 data to TOA reflectance."""
        nodata_mask = np.where(self.arr == 0, 0, 1)  # 0 indicates NaN values
        arr_toa = (self.arr + float(self.params.rescale[self.index]['offset'])) / float(
            self.params.rescale[self.index]['qvalue'])
        return np.where((arr_toa * nodata_mask) == 0, -9999, (arr_toa * nodata_mask))

    def _convert_oli(self) -> np.ndarray:
        """Convert OLI_L8/9 data to TOA reflectance."""
        nodata_mask = np.where(self.arr == 0, 0, 1)  # 0 indicates NaN values
        solar_zenith_rad = np.radians(self.params.geometry[self.index]['solar_zn'])
        arr_toa = ((self.arr * float(self.params.rescale[self.index]['mult'])) + float(
            self.params.rescale[self.index]['add'])) / cos(solar_zenith_rad)
        return np.where((arr_toa * nodata_mask) == 0, -9999, (arr_toa * nodata_mask))

    def _convert_oci(self) -> np.ndarray:
        """Convert OCI_PACE data (already in TOA reflectance)."""
        return self.arr
