import os
import glob
import shutil
import rasterio
import numpy as np
from typing import Dict, List, Any
import xml.etree.ElementTree as ET
from scipy.signal import fftconvolve

from src.satwater.utils import satwutils
from src.satwater.adjacent_correction import adjcorr_functions as adjc_fun

class AdjCorrClass:
    """Handles adjacent correction for satellite imagery."""

    # Constants for band selection
    AQUAVIS_BANDS = {
        "MSI_S2": ["B2", "B02", "B3", "B03", "B4", "B04", "B5", "B05", "B8A", "B8A", "B11"],
        "default": ["B2", "B02", "B3", "B03", "B4", "B04", "B5", "B05", "B6", "B06"]
    }

    def __init__(self):
        """Initialize the adjacent correction processor."""
        self.tree = None
        self.output_path = None
        self.image_path = None

    def fast_convolution(self, image: np.ndarray, kernel: np.ndarray) -> np.ndarray:
        """
        Perform fast convolution using FFT.

        Args:
            image: Input image array
            kernel: Convolution kernel

        Returns:
            Convolved image with same dimensions as input
        """
        # Handle NaNs by converting to zeros
        image = np.nan_to_num(image, nan=0)
        kernel = np.nan_to_num(kernel, nan=0)

        # Calculate padding needed
        image_h, image_w = image.shape
        kernel_h, kernel_w = kernel.shape
        pad_h, pad_w = kernel_h // 2, kernel_w // 2

        # Pad image and perform FFT convolution
        padded_image = np.pad(image, ((pad_h, pad_h), (pad_w, pad_w)),mode='constant', constant_values=0)
        convolved = fftconvolve(padded_image, kernel, mode='same')

        # Crop to original dimensions
        return convolved[pad_h:pad_h + image_h, pad_w:pad_w + image_w]

    def apply_adjcorr(self, input_path: str, output_path: str, params: Dict[str, Any]) -> None:
        """
        Apply adjacent correction to satellite imagery.

        Args:
            input_path: Path to input scene directory
            output_path: Path for output corrected files
            params: Processing parameters dictionary
        """
        self._setup_paths(input_path, output_path)
        self._parse_metadata()

        # Update parameters with metadata info
        params['aux_info']['date_time_info'] = self.metadata["General_Info"]["datetime_image"]
        self._copy_metadata_file()

        # Process each band
        band_info = self._get_band_info()
        for band_key, band_file in band_info.items():
            self._process_single_band(band_key, band_file)

    def _setup_paths(self, input_path: str, output_path: str) -> None:
        """Initialize paths for processing."""
        self.image_path = input_path
        self.output_path = output_path
        os.makedirs(output_path, exist_ok=True)

    def _parse_metadata(self) -> None:
        """Parse and store XML metadata."""
        metadata_file = os.path.join(self.image_path, "MTD.xml")
        self.tree = ET.parse(metadata_file)
        self.metadata = adjc_fun.xml_to_dict(self.tree.getroot())

    def _copy_metadata_file(self) -> None:
        """Copy metadata file to output directory."""
        shutil.copy(os.path.join(self.image_path, "MTD.xml"),
                    os.path.join(self.output_path, "MTD.xml"))

    def _get_band_info(self) -> Dict[str, str]:
        """Get dictionary of band keys and filenames for processing."""
        satellite = self.metadata["General_Info"]["satellite"]
        bands = self.AQUAVIS_BANDS.get(satellite, self.AQUAVIS_BANDS["default"])

        # Filter and update band names
        band_dict = {
            key: value for key, value in self.metadata["General_Info"]["bandname"].items()
            if any(band in value for band in bands)
        }

        # Standardize file extensions
        ext = ".tif" if satellite == "MSI_S2" else ".TIF"
        return {
            key: value.replace(ext, ".tif")
            for key, value in band_dict.items()
        }

    def _process_single_band(self, band_key: str, band_file: str) -> None:
        """Process adjacent correction for a single band."""
        # Read and prepare input image
        with rasterio.open(os.path.join(self.image_path, band_file)) as src:
            image_data = self._prepare_image_data(src)
            crs, transform, nodata = src.crs, src.transform, src.nodata

        # Calculate environmental function and apply convolution
        Fr = adjc_fun.atmospheric_point_scattering_function(self.metadata, band_key)
        convolved_image = self.fast_convolution(image_data, Fr)

        # Calculate rho_env
        ones_array = np.full_like(image_data, 1)
        convolved_Fr = self.fast_convolution(ones_array, Fr)
        rho_env = convolved_image / convolved_Fr

        # Apply adjacency correction
        sr_corr = adjc_fun.Adjacency_correction(
            self.metadata, band_key, image_data, rho_env
        )

        # Save corrected image
        output_file = os.path.join(self.output_path, band_file)
        self._save_corrected_image(sr_corr, crs, transform, nodata, output_file)

    def _prepare_image_data(self, src: rasterio.DatasetReader) -> np.ndarray:
        """Prepare image data by handling nodata values."""
        image_data = src.read(1)
        return np.where(np.isin(image_data, [-9999, 0]), np.nan, image_data)

    def _save_corrected_image(self, data: np.ndarray, crs: Any,
                              transform: Any, nodata: Any, path: str) -> None:
        """Save corrected image to GeoTIFF."""
        with rasterio.open(
                path,
                'w',
                driver='GTiff',
                count=1,
                dtype=data.dtype,
                width=data.shape[1],
                height=data.shape[0],
                crs=crs,
                transform=transform,
                compress='lzw',
                nodata=nodata
        ) as dst:
            dst.write(data, 1)

    def run(self, params: Dict[str, Any]) -> None:
        """
        Run adjacent correction on all scenes in the output directory.

        Args:
            params: Dictionary containing processing parameters
        """
        # Set up input and output paths
        base_path = os.path.join(params['output_dir'], 'atmcor')
        params['output_dir_adjcorr'] = os.path.join(params['output_dir'], 'adjcorr')
        satwutils.create_dir(params['output_dir_adjcorr'])

        # Get input scenes based on satellite type
        input_paths = self._get_input_paths(params, base_path)

        for scene_path in input_paths:
            scene_name = self._get_scene_name(params, scene_path)
            output_path = os.path.join(params['output_dir_adjcorr'], os.path.basename(scene_name))

            if os.path.exists(output_path):
                print(f"Skipping {output_path}, already exists.")
                continue

            satwutils.create_dir(output_path)
            self.apply_adjcorr(scene_name, output_path, params)

    def _get_input_paths(self, params: Dict[str, Any], base_path: str) -> List[str]:
        """Get list of input scene paths based on satellite type."""
        if params['aux_info']['sat_name'] == "sentinel":
            return [os.path.join(base_path, scene) for scene in os.listdir(base_path)]
        return [os.path.join(base_path, scene) for scene in os.listdir(base_path)]

    def _get_scene_name(self, params: Dict[str, Any], scene_path: str) -> str:
        """Get scene name based on satellite type."""
        if params['aux_info']['sat_name'] == "sentinel":
            return glob.glob(os.path.join(scene_path, '*.SAFE*'))[0]
        return glob.glob(os.path.join(scene_path, 'LC*'))[0]