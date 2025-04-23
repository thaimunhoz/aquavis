import os
import glob
import shutil
import logging
import os.path
import rasterio
import rioxarray
import numpy as np
from pathlib import Path
import xml.etree.ElementTree as ET
from rasterio.warp import Resampling
from typing import Dict, List, Optional
from src.aquavis.utils import satwutils

from src.aquavis.utils import io
from src.aquavis.processors.data_class import AquaVisDataLoader
from src.aquavis.processors.qa_flag_band import qa_flag_band
from src.aquavis.processors.ProcessorABC import Processors

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class AQUAVisProductGenerator(Processors):

    def __init__(self):
        self.params = AquaVisDataLoader().load_aquavis_data()

    def convert_float_to_int16(self, input_raster, output_raster, scale_factor=10000, nodata=-9999):
        '''Convert a float raster to int16 by multiplying by a scale factor and setting nodata values.'''

        with rasterio.open(input_raster) as src:
            arr = src.read(1).astype('float32')

            # Clean and prepare data
            arr = np.where((arr > -0.1) & (arr < 0), 0.0001, arr)
            arr[arr <= 0] = np.nan

            # Apply output type scaling
            if self.params.output_type != 'rho':
                arr = arr / np.pi  # Convert to Rrs

            arr *= scale_factor
            arr = np.where(np.isnan(arr), nodata, arr)

            # Write output with updated metadata
            kwargs = src.meta.copy()
            kwargs.update({
                'dtype': rasterio.int16,
                'nodata': nodata,
                'compress': 'lzw'
            })

            with rasterio.open(output_raster, 'w', **kwargs) as dst:
                dst.write(arr.astype(rasterio.int16), 1)

    def _parse_datetime(self) -> str:
        """Parse and format datetime from parameters."""
        dt_info = self.img_params['General_Info']['datetime_image']
        day_part, time_part = dt_info.split("T")
        year, month, day = day_part.split("-")
        hour, minute, second = time_part.split(":")
        return f"{year}{month}{day}T{hour}{minute}{second}"

    def _generate_product_name(self, band: Optional[str] = None) -> str:
        """Generate standardized product filename."""
        date_time = self._parse_datetime()
        parts = [
            "AQUAVis",
            f"T{self.params.tile}",
            date_time,
            self.params.ncode,
            band if band else "QAflag",
            "v1.0"
        ]
        return "_".join(parts) + ".tif"

    def _process_cloud_data(self, scene_path: str) -> None:
        """Process and save cloud information for the scene."""
        flags_dir = Path(self.params.output_dir_qaflag) / Path(scene_path).name
        flags_dir.mkdir(parents=True, exist_ok=True)
        output_path = flags_dir / f"{Path(scene_path).name}_cloud_shadow.tif"

        try:
            if self.params.select_sat == 'landsat':
                cloud_path = os.path.join(self.img_params["General_Info"]["atmcor_folder"], "SatClouds", "temp", "cloud.tif")
                shutil.copy(cloud_path, output_path)
            else:
                #scene_name = Path(scene_path).name.replace('.SAFE', '')
                cloud_path = os.path.join(self.img_params["General_Info"]["atmcor_folder"], "SatClouds", "temp", "cloud10.tif")
                data = rioxarray.open_rasterio(cloud_path)
                resampled_data = data.rio.reproject(
                    data.rio.crs,
                    resolution=(30, 30),
                    resampling=Resampling.bilinear
                )
                resampled_data.rio.to_raster(output_path)
        except Exception as e:
            print(f"Cloud processing failed for {scene_path}: {str(e)}")

    def _process_bands(self, sentinel_bands: List, scene_path: str, out_dir_hls: str) -> None:
        """Process all bands for a scene."""
        scene_all_bands = sorted(
            glob.glob(f"{scene_path}/**.tif"),
            key=lambda x: next((i for i, band in enumerate(sentinel_bands)
                                if band in x), float('inf'))
        )

        for file_band in scene_all_bands:
            band_nmi = Path(file_band).stem.split("_")[-1]
            band_nmi = 'B05' if band_nmi == 'B8A' else band_nmi
            band_nmi = satwutils.fix_band_name(band_nmi, self.params.select_sat)

            out_filename = self._generate_product_name(band_nmi)
            out_path = Path(out_dir_hls) / out_filename
            self.convert_float_to_int16(file_band, str(out_path))

    def _clean_working_directories(self) -> None:
        """Clean up working directories after processing."""
        folders_to_clean = [
            self.params.output_dir_atmcor,
            self.params.output_dir_qaflag,
            self.params.output_dir_adjcorr,
            self.params.output_dir_glintcorr,
            self.params.output_dir_tiling,
            self.params.output_dir_water_mask
        ]

        for folder in folders_to_clean:
            try:
                os.chmod(folder, 0o777)
                shutil.rmtree(folder)
            except (OSError, shutil.Error):
                continue

    def _get_scenes(self, params: Dict[str, str]) -> List[str]:

        sat = self.params.select_sat

        if sat == 'landsat':
            ncode = 'L30'
            path_pr = fr'{params["output_dir_wm"]}'
            scenes = [fr'{path_pr}\{i}' for i in os.listdir(path_pr)]
        else:
            ncode = 'S30'
            path_pr = fr'{params["output_dir_wm"]}'
            scenes = [fr'{path_pr}\{i}' for i in os.listdir(path_pr)]

        params['ncode'] = ncode
        params["sat"] = sat

        return scenes

    def run(self, scene_path: str, output_path: str) -> None:

        io.validate_file(self.params, scene_path)

        mtd_file_path = os.path.join(self.params.output_dir_adjcorr, os.path.basename(scene_path), "MTD.xml")
        tree = ET.parse(mtd_file_path)
        root = tree.getroot()
        self.img_params = satwutils.xml_to_dict(root)

        date_time = self._parse_datetime()
        out_dir_hls = os.path.join(self.params.output_dir_aquavis, f"AQUAVis_T{self.params.tile}_{date_time}_{self.params.ncode}_v1.0")
        os.makedirs(out_dir_hls, exist_ok=True)

        sentinel_bands = ['B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B8A', 'B11', 'B12']
        self._process_bands(sentinel_bands, scene_path, out_dir_hls)

        self._process_cloud_data(scene_path)

        # Generate QA flag band
        # qa_band = qa_flag_band.run_qa_flag(self.params, scene_path)
        # qa_path = out_dir_hls / self._generate_product_name()
        # qa_band.rio.to_raster(str(qa_path))

        logger.info("AQUAVis completed successfully!")