import os
import ee
import geemap
import logging
import rasterio
import numpy as np
import geopandas as gpd
import rioxarray as rxr
from typing import Optional
from rasterio.crs import CRS
from rasterio.mask import mask
from rasterio.merge import merge
from rasterio.features import shapes
from rasterio.transform import Affine
from shapely.geometry import box, shape
from rasterio.features import rasterize

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class WaterMaskClass:

    """
    A class to extract water masks from satellite imagery an image-based approach.
    """

    def __init__(self):
        pass

    def raster2shp(
        self, raster_numpy: np.ndarray, raster_transform: Affine, raster_crs: CRS
    ) -> gpd.GeoDataFrame:
        """
        Convert a raster array to a shapefile.

        Args:
            raster_numpy (np.ndarray): Raster array.
            raster_transform (Affine): Raster transform.
            raster_crs (CRS): Raster coordinate reference system.

        Returns:
            gpd.GeoDataFrame: GeoDataFrame containing the vectorized raster.
        """
        results = list(
            {"properties": {"raster_val": v}, "geometry": s}
            for s, v in shapes(np.asarray(raster_numpy, dtype=np.int16), transform=raster_transform)
            if v  # Only take shapes with raster_val = True (i.e., v=1)
        )
        geometries = [shape(feature["geometry"]) for feature in results]

        return gpd.GeoDataFrame(geometry=geometries, crs=raster_crs)

    def buffer_into_polygon(self, gdf, buffer_size):

        '''
        Buffer the water mask and keep only the biggest polygon
        Input:
            gdf (GeoDataFrame): input shapefile
            buffer_size (int): buffer size
        Return:
            str: path to the modified water mask shapefile
        '''

        gdf['geometry'] = gdf.geometry.buffer(-buffer_size)
        gdf_cleaned = gdf[~gdf['geometry'].is_empty & gdf['geometry'].is_valid]

        # Keep only the polygons with the biggest area
        gdf_cleaned['area'] = gdf_cleaned['geometry'].area
        gdf_largest = gdf_cleaned.sort_values(by='area', ascending=False).head(1)

        return gdf_largest

    def get_mask_water(self, img_path: str) -> None:

        """
        Generate a water mask using an image-based approach - combination of SWIR threshold and MNDWI index.

        Args:
            img_path (str): Path to the input image folder.

        Returns:
            np.ndarray: Water mask array.
        """

        band_path = [i for i in os.listdir(img_path) if i.endswith(".TIF")]
        green_band = next((band for band in band_path if "B3" in band or "B03" in band), None)
        swir_band = next((band for band in band_path if "B11" in band), (band for band in band_path if "B6" in band or "B06" in band))

        xda_green = rxr.open_rasterio(os.path.join(img_path, green_band))
        xda_swir = rxr.open_rasterio(os.path.join(img_path, swir_band))

        xda_green_matched = xda_green.rio.reproject_match(xda_swir)

        swir_mask = xda_swir < 0.03
        mndwi = (xda_green_matched - xda_swir) / (xda_green_matched + xda_swir)
        mndwi_mask = mndwi > 0.3
        water_mask = swir_mask & mndwi_mask

        with rasterio.open(os.path.join(img_path, green_band)) as src:
            _transform = src.transform
            _crs = src.crs

        # Include adjacency buffer
        water_shp = self.raster2shp(water_mask, _transform, _crs)

        water_shp_buffered = self.buffer_into_polygon(water_shp, 300)

        return water_shp_buffered

    def clip_water(self, img_path: str, dest_path: str) -> None:

        """
        Return only the water extend as valid observation. Land pixels are set no value.

        Args:
            img_path (str): Path to the input image folder.
        Returns:
            clipped water images
        """

        files = [folder for folder in os.listdir(img_path) if os.path.isdir(os.path.join(img_path, folder)) and "SatCloud" not in folder][0]

        img_path = os.path.join(img_path, files)
        water_shp_buffered = self.get_mask_water(img_path)
        band_paths = [i for i in os.listdir(img_path) if i.endswith(".TIF")]

        swir_band = next((band for band in band_paths if "B11" in band), (band for band in band_paths if "B06" in band or "B6" in band))
        xda_swir = rxr.open_rasterio(os.path.join(img_path, swir_band))

        for band in band_paths:

            with rasterio.open(band) as src:

                band_data = src.read(1)
                out_meta = src.meta.copy()
                out_meta.update(nodata=-9999, dtype=band_data.dtype)

                # Rasterize the water mask to match the raster size
                water_raster = rasterize(
                    [(geom, 1) for geom in water_shp_buffered.geometry],  # Assign value 1 to water pixels
                    out_shape=(src.height, src.width),
                    transform=src.transform,
                    fill=0,  # Background (land) will be 0
                    dtype="uint8"
                )

                # Apply glint correction based on SWIR subtraction
                glint_corr = np.where((band_data - xda_swir) < 0, band_data, (band_data - xda_swir))

                masked_band = np.where(water_raster == 1, glint_corr, -9999)

                profile = src.profile
                profile.update(dtype=band_data.dtype, nodata=src.nodata)

                # Save clipped image
                output_path = os.path.join(dest_path, os.path.basename(band))

                with rasterio.open(output_path, "w", **profile) as dst:
                    dst.write(masked_band, 1)
