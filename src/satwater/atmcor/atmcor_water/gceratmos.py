import os
import glob
import warnings
import rasterio
import numpy as np
import rioxarray as rxr
from typing import Optional
from datetime import datetime

from src.satwater.atmcor.atmcor_water import toolbox as tool
from src.satwater.atmcor.atmcor_water.atm.atmosphere import Atmosphere
from src.satwater.atmcor.atmcor_water.atm.correction import Correction
from src.satwater.atmcor.atmcor_water.metadata_msi import Metadata_MSI_S2
from src.satwater.atmcor.atmcor_water.metadata_oli import Metadata_OLI_L89

class GCERAtmos:

    """Handles atmospheric correction for satellite imagery."""

    GRANULE_DIR = '/GRANULE'
    IMG_DATA_DIR = '/IMG_DATA'

    def __init__(
            self,
            path_main: str,
            path_dest: str,
            networkdrive_letter: str,
            satellite: str,
            aero_type: str,
            mode: Optional[str] = None,
            msi_tile: Optional[str] = None,
            path_buffer: Optional[str] = None
    ):
        """
        Initialize the atmospheric corrector.

        Args:
            path_main: Path to main input directory
            path_dest: Path to output directory
            networkdrive_letter: Network drive letter
            satellite: Satellite type ('MSI_S2' or 'OLI_L8/9')
            aero_type: Aerosol type
            mode: Processing mode (optional)
            msi_tile: MSI tile identifier (optional)
            path_buffer: Buffer path (optional)
        """
        self.path_main = path_main
        self.path_dest = path_dest
        self.networkdrive_letter = networkdrive_letter
        self.satellite = satellite.upper()
        self.aero_type = aero_type
        self.mode = mode
        self.msi_tile = msi_tile
        self.path_buffer = path_buffer

    def run(self) -> None:
        """Execute the atmospheric correction process."""
        processor = {
            'MSI_S2': self._process_sentinel,
            'OLI_L8/9': self._process_landsat
        }.get(self.satellite, self._handle_unknown_satellite)

        processor()

    def _process_sentinel(self) -> None:
        """Process Sentinel-2 MSI data."""
        print(self.path_main[-65:])

        # Setup directories
        dest_dir = tool.newdirectory(self.path_dest, self.path_main[-65:])
        temp_dir = tool.newdirectory(dest_dir, 'tempdir')

        # Process JP2 files
        self._convert_jp2_to_tiff(temp_dir)

        # Process metadata and atmospheric correction
        meta = Metadata_MSI_S2(
            self.path_main, dest_dir, self.networkdrive_letter,
            self.satellite, self.aero_type, self.mode, self.msi_tile
        )
        meta.run()

        self._perform_atmospheric_correction_msi(meta, temp_dir, dest_dir)
        # shutil.rmtree(temp_dir)  # Uncomment when ready to clean up

    def _process_landsat(self) -> None:
        """Process Landsat 8/9 OLI data."""
        print(self.path_main[-40:])

        output_dir = tool.newdirectory(self.path_dest, self.path_main[-40:])

        # Process metadata
        meta = Metadata_OLI_L89(
            self.path_main, output_dir, self.networkdrive_letter,
            self.satellite, self.aero_type, self.mode, self.msi_tile
        )
        meta.run()

        self._perform_atmospheric_correction_oli(meta, output_dir)

    def _handle_unknown_satellite(self) -> None:
        """Handle unknown satellite types."""
        warnings.warn(
            f"The sensor type '{self.satellite}' was not identified in GCERATMOS.",
            UserWarning
        )

    def _convert_jp2_to_tiff(self, temp_dir: str) -> None:
        """Convert JP2 files to TIFF format."""
        granule_path = glob.glob(os.path.join(self.path_main + self.GRANULE_DIR, '*L1C_*'))[0]
        band_jp2_files = glob.glob(os.path.join(granule_path + self.IMG_DATA_DIR, '*_B*.jp2'))

        for jp2_file in band_jp2_files:
            output_name = os.path.join(temp_dir, f"{os.path.basename(jp2_file)[:-4]}.tif")
            tool.jp2_to_tiff_xarray(jp2_file, output_name)

    def _perform_atmospheric_correction_msi(self, meta: object, temp_dir: str, dest_dir: str) -> None:
        """Perform atmospheric correction on all bands."""
        atmos_param = Atmosphere(meta)
        atmos_param.run()

        for index, band in enumerate(meta.bandname):
            input_path = os.path.join(temp_dir, f"{band[:-4]}.tif")
            output_path = os.path.join(dest_dir, f"{band[:-4]}.tif")

            print(input_path)  # Debug output

            arr = rxr.open_rasterio(input_path).squeeze().values.astype(float)
            corr = Correction(meta, atmos_param, arr, index)
            corr.run()

            arr_corrected = np.where(corr.arr_sr < 0, -9999, corr.arr_sr)
            tool.export(arr_corrected, band, input_path, output_path)

        tool.export_meta(meta, atmos_param, dest_dir)

    def _perform_atmospheric_correction_oli(self, meta: object, output_dir: str) -> None:
        """Perform atmospheric correction on all bands."""
        atmos_param = Atmosphere(meta)
        atmos_param.run()

        for index, band in enumerate(meta.bandname):
            input_path = os.path.join(meta.path_main, f"{band[:-4]}.tif")
            output_path = os.path.join(output_dir, f"{band[:-4]}.tif")

            with rasterio.open(input_path) as src:
                arr = src.read(1)
                out_meta = src.meta.copy()

            corr = Correction(meta, atmos_param, arr, index)
            corr.run()

            arr_corrected = np.where(corr.arr_sr < 0, -9999, corr.arr_sr)

            self._save_corrected_image(arr_corrected, out_meta, output_path)

        tool.export_meta(meta, atmos_param, output_dir)

    def _save_corrected_image(
            self,
            array: np.ndarray,
            metadata: dict,
            output_path: str
    ) -> None:
        """Save corrected image with updated metadata."""
        metadata.update({
            "driver": "GTiff",
            "height": array.shape[0],
            "width": array.shape[1],
            "dtype": 'float64'
        })

        with rasterio.open(output_path, "w", **metadata) as dest:
            dest.write(array, 1)

    @staticmethod
    def _parse_date_from_filename(filename: str, satellite_type: str) -> datetime:
        """Extract and parse date from satellite image filename."""
        basename = os.path.basename(filename)

        if satellite_type.upper() == 'SENTINEL':
            date_str = basename.split('_')[2].split('T')[0]
        else:  # Landsat
            date_str = basename.split('_')[3]

        return datetime.strptime(date_str, '%Y%m%d')