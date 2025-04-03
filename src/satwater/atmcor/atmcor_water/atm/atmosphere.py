import os
import warnings
import numpy as np
import pandas as pd
import Py6S as py6s
from typing import Dict, List, Tuple

class Atmosphere:
    """Handles atmospheric correction using 6S radiative transfer model."""

    # Constants
    RSR_FILE_PATHS = {
        'MSI_S2': {
            'S2A': 'rsr_S2A_MSI.txt',
            'S2B': 'rsr_S2B_MSI.txt',
            'wavelength_range': [0.412, 2.321],
            'band_count': 13
        },
        'OLI_L8/9': {
            'LC08': 'rsr_L8_OLI.txt',
            'OLI2': 'rsr_L9_OLI2.txt',
            'wavelength_range': [0.427, 2.356],
            'band_count': 8
        }
    }

    BASE_RSR_PATH = os.path.join('src', 'satwater', 'atmcor', 'atmcor_water', 'atm')

    def __init__(self, parameters):
        """Initialize with atmospheric parameters."""
        self.parameters = parameters
        self.values: Dict[int, Dict[str, float]] = {}
        self.values_adjcorr: Dict[int, Dict[str, float]] = {}

    def run(self) -> None:
        """Run atmospheric correction for all bands."""
        # Initialize 6S model
        s = self._initialize_sixs_model()

        # Load and process RSR data
        bands_data = self._load_rsr_data()

        # Run simulation for each band
        for band_idx in range(len(bands_data)):
            band_data = bands_data[band_idx + 1]  # Bands are 1-indexed
            self._run_band_simulation(s, band_idx, band_data)

    def _initialize_sixs_model(self) -> py6s.SixS:
        """Initialize and configure the 6S model."""
        s = py6s.SixS()

        # Set atmospheric profile
        s.atmos_profile = py6s.AtmosProfile.UserWaterAndOzone(
            self.parameters.water_vapour,
            self.parameters.ozone
        )

        # Set aerosol properties
        s.aot550 = self.parameters.aod
        aero_profile = (
            py6s.AeroProfile.Continental if self.parameters.aero_type == 'Continental'
            else py6s.AeroProfile.Maritime
        )
        s.aero_profile = py6s.AeroProfile.PredefinedType(aero_profile)

        # Set altitude configuration
        s.altitudes = py6s.Altitudes()
        s.altitudes.set_sensor_satellite_level()
        s.altitudes.set_target_custom_altitude(self.parameters.altitude)

        return s

    def _load_rsr_data(self) -> Dict[int, Tuple[List[float], List[float]]]:
        """Load and process relative spectral response data."""
        satellite_config = self.RSR_FILE_PATHS.get(self.parameters.satellite)
        if not satellite_config:
            warnings.warn(f"Sensor type {self.parameters.satellite} not recognized", UserWarning)
            return {}

        rsr_file = satellite_config.get(self.parameters.type)
        if not rsr_file:
            warnings.warn(f"Sensor subtype {self.parameters.type} not recognized", UserWarning)
            return {}

        rsr_path = os.path.join(self.BASE_RSR_PATH, rsr_file)
        try:
            rsr_data = pd.read_csv(
                rsr_path,
                sep='\t',
                skiprows=5,
                names=['wavelength', 'rsr', 'band', 'id']
            )
            return self._interpolate_rsr(
                rsr_data,
                satellite_config['band_count'],
                satellite_config['wavelength_range']
            )
        except FileNotFoundError:
            warnings.warn(f"RSR file not found at {rsr_path}", UserWarning)
            return {}

    def _run_band_simulation(
            self,
            s: py6s.SixS,
            band_idx: int,
            band_data: Tuple[List[float], List[float]]
    ) -> None:
        """Run 6S simulation for a single band."""
        # Set geometry parameters
        s.geometry = py6s.Geometry.User()
        s.geometry.day = self.parameters.datetime.day
        s.geometry.month = self.parameters.datetime.month
        s.geometry.solar_a = self.parameters.geometry[band_idx]['solar_az']
        s.geometry.solar_z = self.parameters.geometry[band_idx]['solar_zn']
        s.geometry.view_a = self.parameters.geometry[band_idx]['view_az']
        s.geometry.view_z = self.parameters.geometry[band_idx]['view_zn']

        # Set wavelength parameters
        min_wavelength = np.min(band_data[1])
        max_wavelength = np.max(band_data[1])
        s.wavelength = py6s.Wavelength(min_wavelength, max_wavelength, band_data[0])

        # Run simulation
        s.run()

        # Store results
        self._store_simulation_results(s, band_idx)

    def _store_simulation_results(self, s: py6s.SixS, band_idx: int) -> None:
        """Store simulation results for a band."""
        geo = self.parameters.geometry[band_idx]

        self.values[band_idx] = {
            'view_zn': geo['view_zn'],
            'view_az': geo['view_az'],
            'solar_zn': geo['solar_zn'],
            'solar_az': geo['solar_az'],
            'tg_OG_co': float(s.outputs.transmittance_co.total),
            'tg_OG_c02': float(s.outputs.transmittance_co2.total),
            'tg_OG_o2': float(s.outputs.transmittance_oxygen.total),
            'tg_OG_no2': float(s.outputs.transmittance_no2.total),
            'tg_OG_ch4': float(s.outputs.transmittance_ch4.total),
            'Tg_H20': float(s.outputs.transmittance_water.total),
            'Tg_O3': float(s.outputs.transmittance_ozone.total),
            'p_atm': float(s.outputs.atmospheric_intrinsic_reflectance),
            'optical_depth_rayleigh': float(s.outputs.optical_depth_total.rayleigh)
        }

        self.values_adjcorr[band_idx] = {
            'view_z': geo['view_zn'],
            'view_az': geo['view_az'],
            'solar_zn': geo['solar_zn'],
            'solar_az': geo['solar_az'],
            'optical_depth__total_Ray': s.outputs.optical_depth_total.rayleigh,
            'rayleigh_scatransmi_upward': s.outputs.transmittance_rayleigh_scattering.upward,
            'optical_depth__total_Aero': s.outputs.optical_depth_total.aerosol,
            'aerosol_scatransmi_upward': s.outputs.transmittance_aerosol_scattering.upward,
            'optical_depth__total_AeroRay': float(s.outputs.optical_depth_total.aerosol + s.outputs.optical_depth_total.rayleigh),
            'total_scattering_transmittance_upward': s.outputs.transmittance_total_scattering.upward,
            'optical_depth_total': float(s.outputs.optical_depth_total.total)
        }

    @staticmethod
    def _interpolate_rsr(
            rsr: pd.DataFrame,
            band_number: int,
            wavelength_range: List[float]
    ) -> Dict[int, Tuple[List[float], List[float]]]:
        """
        Interpolate relative spectral response data to 2.5 nm resolution.

        Args:
            rsr: DataFrame containing RSR data
            band_number: Number of bands to process
            wavelength_range: Min and max wavelength values

        Returns:
            Dictionary mapping band numbers to (response, wavelength) pairs
        """
        output = {}
        w_min, w_max = wavelength_range

        for band_id in range(1, band_number + 1):
            band_data = rsr.loc[rsr['id'] == band_id]

            # Original wavelength and response values
            orig_wavelengths = band_data['wavelength'].values
            orig_response = band_data['rsr'].values

            # Create complete wavelength range with 0 response outside band range
            complete_wavelengths = np.arange(w_min, w_max, 0.001)
            complete_response = np.zeros_like(complete_wavelengths)

            # Find overlapping indices
            mask = (complete_wavelengths >= orig_wavelengths.min()) & \
                   (complete_wavelengths <= orig_wavelengths.max())
            complete_response[mask] = np.interp(
                complete_wavelengths[mask],
                orig_wavelengths,
                orig_response
            )

            # Interpolate to 2.5 nm resolution
            interp_wavelengths = np.around(np.arange(w_min, w_max, 0.0025), 4)
            interp_response = np.interp(
                interp_wavelengths,
                complete_wavelengths,
                complete_response
            )

            output[band_id] = (interp_response.tolist(), interp_wavelengths.tolist())

        return output