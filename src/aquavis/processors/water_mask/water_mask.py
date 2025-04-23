import os
import logging
import rasterio
import numpy as np
import geopandas as gpd
import rioxarray as rxr
from rasterio.crs import CRS
from shapely.geometry import shape
from rasterio.features import shapes
from rasterio.transform import Affine
from rasterio.features import rasterize

from src.aquavis.utils import io
from src.aquavis.config.config import SatWaterConfig
from src.aquavis.processors.ProcessorABC import Processors
from src.aquavis.processors.data_class import AquaVisDataLoader

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class RunWaterMask(Processors):

    def __init__(self):
        self.params = AquaVisDataLoader().load_aquavis_data()
        self.paths = SatWaterConfig()._load_paths()
        pass

    def raster2shp(self, raster_numpy: np.ndarray, raster_transform: Affine, raster_crs: CRS) -> gpd.GeoDataFrame:

        """Convert a raster array to a shapefile."""

        results = list(
            {"properties": {"raster_val": v}, "geometry": s}
            for s, v in shapes(np.asarray(raster_numpy, dtype=np.int16), transform=raster_transform)
            if v  # Only take shapes with raster_val = True (i.e., v=1)
        )

        geometries = [shape(feature["geometry"]) for feature in results]

        return gpd.GeoDataFrame(geometry=geometries, crs=raster_crs)

    def get_mask_water(self, img_path: str) -> gpd.GeoDataFrame:
        """Generate a water mask using an image-based approach - combination of SWIR threshold and MNDWI index."""

        band_path = [i for i in os.listdir(img_path) if i.endswith(".tif")]
        green_band = next((band for band in band_path if "B3" in band or "B03" in band), None)
        swir_band = next((band for band in band_path if "B11" in band or "B6" in band), None)

        xda_green = rxr.open_rasterio(os.path.join(img_path, green_band))
        xda_swir = rxr.open_rasterio(os.path.join(img_path, swir_band))

        xda_green_matched = xda_green.rio.reproject_match(xda_swir)

        xda_swir = xda_swir.where(xda_swir > 0, 0.0001)
        xda_green_matched = xda_green_matched.where(xda_green_matched > 0, 0.0001)

        #swir_mask = xda_swir < 0.03
        mndwi = (xda_green_matched - xda_swir) / (xda_green_matched + xda_swir)
        mndwi_mask = mndwi > 0
        #water_mask = swir_mask & mndwi_mask

        with rasterio.open(os.path.join(img_path, green_band)) as src:
            _transform = src.transform
            _crs = src.crs

        # Include adjacency buffer
        water_shp = self.raster2shp(mndwi_mask, _transform, _crs)

        return water_shp

    def clip_water(self, img_path: str) -> None:
        """Return only the water extend as valid observation. Land pixels are set no value."""

        water_aux = self.get_mask_water(img_path)

        #water_shp_buffered = gpd.GeoDataFrame(geometry=[water_aux], crs=water_aux.crs)

        band_paths = [i for i in os.listdir(img_path) if i.endswith(".tif")]

        for band in band_paths:

            with rasterio.open(os.path.join(img_path, band)) as src:

                band_data = src.read(1)
                out_meta = src.meta.copy()
                out_meta.update(nodata=-9999, dtype=band_data.dtype)

                band_data[band_data <= 0] = 0.0001
                band_data = np.nan_to_num(band_data, nan=0.0001)

                # Rasterize the water mask to match the raster size
                water_raster = rasterize(
                    [(geom, 1) for geom in water_aux.geometry],  # Assign value 1 to water pixels
                    out_shape=(src.height, src.width),
                    transform=src.transform,
                    fill=0,  # Background (land) will be 0
                    dtype="uint8"
                )

                masked_band = np.where(water_raster == 1, band_data, -9999)
                final_image = np.where(masked_band < 0, -9999, masked_band)

                profile = src.profile
                profile.update(dtype=band_data.dtype, nodata=src.nodata)

                # Save clipped image
                if not os.path.exists(os.path.join(self.params.output_dir_water_mask, os.path.basename(img_path))):
                    os.makedirs(os.path.join(self.params.output_dir_water_mask, os.path.basename(img_path)))

                if not os.path.exists(os.path.join(self.params.output_dir_qaflag, os.path.basename(img_path))):
                    os.makedirs(os.path.join(self.params.output_dir_qaflag, os.path.basename(img_path)))

                output_path = os.path.join(self.params.output_dir_water_mask, os.path.basename(img_path), os.path.basename(band))

                with rasterio.open(output_path, "w", **profile) as dst:
                    dst.write(final_image, 1)

                with rasterio.open(os.path.join(os.path.join(self.params.output_dir_qaflag, os.path.basename(img_path)), f"{os.path.basename(img_path)}_water.tif"), "w", **profile) as dst:
                    dst.write(water_raster, 1)

    def run(self, path_main: str, path_dest: str) -> None:

        io.validate_file(self.params, path_main)

        self.clip_water(path_main)

        logger.info("Water mask completed successfully")