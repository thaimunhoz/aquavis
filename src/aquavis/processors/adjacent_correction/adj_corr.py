import os
import glob
import shutil
import logging
import rasterio
import numpy as np
import xml.etree.ElementTree as ET
from typing import Dict, Any, Optional

from src.aquavis.utils import io
from src.aquavis.processors.ProcessorABC import Processors
from src.aquavis.processors.adjacent_correction import adjcorr_toolbox as adjc_fun
from src.aquavis.processors.data_class import AquaVisDataLoader

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class AdjCorrClass(Processors):

    """Handles adjacent correction for satellite imagery."""

    # Constants for band selection
    AQUAVIS_BANDS = {
        "MSI_S2": ["B2", "B02", "B3", "B03", "B4", "B04", "B5", "B05", "B8A", "B8A", "B11"],
        "default": ["B2", "B02", "B3", "B03", "B4", "B04", "B5", "B05", "B6", "B06"]
    }

    def __init__(self):

        """Initialize the adjacent correction processor."""

        self.metadata = AquaVisDataLoader().load_aquavis_data()

    def Adjacency_correction(self, band, array_band, adjc_array):

        """
        Removes the adjacency effect of the image, using the equation described in the Vermote et al. (1997).
        """

        # Zenith view angle -> degree:
        view_z = float(self.atmosphere_parameters["InputData"]['sixSV_params'][band]['view_z'])

        # Atmopheric optical depth (Rayleigh + Aerosol) -> Atmospheric_OpticalDepth:
        Atmospheric_OpticalDepth = float(self.atmosphere_parameters["InputData"]['sixSV_params'][band]['optical_depth__total_AeroRay'])

        # Total transmittance UPWARD (Rayleigh + Aerosol) -> T_upward:
        T_upward = float(self.atmosphere_parameters["InputData"]['sixSV_params'][band]['total_scattering_transmittance_upward'])

        # Total transmittance UPWARD direct (Rayleigh + Aerosol) -> T_upward_dirAeroRay:
        T_upward_dirAeroRay = np.exp(-Atmospheric_OpticalDepth / np.cos(view_z * (np.pi / 180)))

        # Total transmittance UPWARD diffuse (Rayleigh + Aerosol) -> T_upward_difAeroRay
        T_upward_difAeroRay = T_upward - T_upward_dirAeroRay

        # Surface reflectance without adjacency effect - Vermote et al. (1997):
        sr_corr = array_band - (adjc_array * T_upward_difAeroRay)

        return sr_corr

    def atmospheric_point_scattering_function(self, band):

        """Calculates the Fr weight per distance."""

        # Converts the grid distance to km:
        if self.atmosphere_parameters["General_Info"]["satellite"] == "MSI_S2":
            if band == "B11":  # 20m
                grid_matrix = adjc_fun.create_grid(30, 20)
            else:  # 10m
                grid_matrix = adjc_fun.create_grid(60, 10)
        else:  # 30m
            grid_matrix = adjc_fun.create_grid(20, 30)

        radius_km = grid_matrix / 1000

        # Zenith view angle -> degree:
        view_z = float(self.atmosphere_parameters["InputData"]['sixSV_params'][band]['view_z'])

        # Rayleigh UPWARD diffuse transmittance -> T_upward_difRayleigh:
        Rayleigh_OpticalDepth = float(self.atmosphere_parameters["InputData"]['sixSV_params'][band]['optical_depth__total_Ray'])
        T_upward_Rayleigh = float(self.atmosphere_parameters["InputData"]['sixSV_params'][band]['rayleigh_scatransmi_upward'])
        T_upward_dirRayleigh = np.exp(-Rayleigh_OpticalDepth / np.cos(view_z * (np.pi / 180)))
        T_upward_difRayleigh = T_upward_Rayleigh - T_upward_dirRayleigh

        # Aerosol UPWARD diffuse transmittance -> T_upward_difAerosol:
        Aerosol_OpticalDepth = float(self.atmosphere_parameters["InputData"]['sixSV_params'][band]['optical_depth__total_Aero'])
        T_upward_Aerosol = float(self.atmosphere_parameters["InputData"]['sixSV_params'][band]['aerosol_scatransmi_upward'])
        T_upward_dirAerosol = np.exp(-Aerosol_OpticalDepth / np.cos(view_z * (np.pi / 180)))
        T_upward_difAerosol = T_upward_Aerosol - T_upward_dirAerosol

        # Calculates the Aerosol's Fr and Rayleigh's Fr functions using the equation described by Vermote et al.(2006):
        FrRayleigh = ((0.930 * np.exp(-0.08 * radius_km)) + (0.070 * np.exp(-1.10 * radius_km)))
        FrAerosol = ((0.448 * np.exp(-0.27 * radius_km)) + (0.552 * np.exp(-2.83 * radius_km)))

        # Calculates the APSF (Fr) -> Atmospheric Point Scattering Function:
        Fr = (T_upward_difRayleigh * FrRayleigh + T_upward_difAerosol * FrAerosol) / (T_upward_difRayleigh + T_upward_difAerosol)

        return Fr

    def _get_band_info(self) -> Dict[str, str]:
        """Get dictionary of band keys and filenames for processing."""
        satellite = self.atmosphere_parameters["General_Info"]["satellite"]
        bands = self.AQUAVIS_BANDS.get(satellite, self.AQUAVIS_BANDS["default"])

        # Filter and update band names
        band_dict = {
            key: value for key, value in self.atmosphere_parameters["General_Info"]["bandname"].items()
            if any(band in value for band in bands)
        }

        # Standardize file extensions
        ext = ".jp2" if satellite == "MSI_S2" else ".TIF"
        return {
            key: value.replace(ext, ".tif")
            for key, value in band_dict.items()
        }

    def _process_single_band(self, input_path: str, output_path: str, band_key: str, band_file: str) -> None:

        """Process adjacent correction for a single band."""

        # Read and prepare input image
        with rasterio.open(os.path.join(input_path, band_file)) as src:
            image_data = adjc_fun._prepare_image_data(src)
            crs, transform, nodata = src.crs, src.transform, src.nodata

        # Calculate environmental function and apply convolution
        Fr = self.atmospheric_point_scattering_function(band_key)
        convolved_image = adjc_fun.fast_convolution(image_data, Fr)

        # Calculate rho_env
        ones_array = np.full_like(image_data, 1)
        convolved_Fr = adjc_fun.fast_convolution(ones_array, Fr)
        rho_env = convolved_image / convolved_Fr

        # Apply adjacency correction
        sr_corr = self.Adjacency_correction(band_key, image_data, rho_env)

        # Save corrected image
        output_file = os.path.join(output_path, band_file)
        adjc_fun._save_corrected_image(sr_corr, crs, transform, nodata, output_file)

    def run(self, input_path: str, output_path: str) -> None:

        """Apply adjacent correction to satellite imagery."""
        os.makedirs(output_path, exist_ok=True)

        if self.metadata.select_sat == "sentinel":
            scene_name = glob.glob(os.path.join(input_path, '*.SAFE*'))[0]
        else:
            scene_name = glob.glob(os.path.join(input_path, 'LC*'))[0]

        io.validate_file(self.metadata, scene_name)

        tree = ET.parse(os.path.join(scene_name, "MTD.xml"))
        root = tree.getroot()
        shutil.copy2(os.path.join(scene_name, "MTD.xml"), os.path.join(output_path, "MTD.xml"))
        self.atmosphere_parameters = adjc_fun.xml_to_dict(root)

        band_info = self._get_band_info()

        for band_key, band_file in band_info.items():
            self._process_single_band(scene_name, output_path, band_key, band_file)

        logger.info("Adjacency correction completed successfully")