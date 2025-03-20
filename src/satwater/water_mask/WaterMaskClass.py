import os
import logging
import rasterio
import numpy as np
import multiprocessing
import geopandas as gpd
import rioxarray as rxr
from rasterio.crs import CRS
from shapely.geometry import shape
from rasterio.features import shapes
from rasterio.transform import Affine
from rasterio.features import rasterize

from src.satwater.water_mask import generate_water_mask as wm

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

        band_path = [i for i in os.listdir(img_path) if i.endswith(".tif")]
        green_band = next((band for band in band_path if "B3" in band or "B03" in band), None)
        swir_band = next((band for band in band_path if "B11" in band or "B6" in band), None)

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

        #water_shp_buffered = self.buffer_into_polygon(water_shp, 300)

        return water_shp

    def clip_water(self, img_path: str, params) -> None:

        """
        Return only the water extend as valid observation. Land pixels are set no value.

        Args:
            img_path (str): Path to the input image folder.
        Returns:
            clipped water images
        """
        select_sat = params['aux_info']['sat_name']

        tiles = [params[select_sat].get('tiles', [])][0]
        sentinel_tiles_path = r"C:\Users\tml411\Documents\Python Scripts\hls_water\src\satwater\auxfiles\tiles\MGRS_tiles.shp"
        msi_tiles_gdf = gpd.read_file(sentinel_tiles_path)
        tile_gdf = msi_tiles_gdf[msi_tiles_gdf['Name'] == tiles]

        water_aux = self.get_mask_water(img_path)
        water_1 = water_aux.geometry.unary_union
        water_2 = wm.create_water_mask(img_path, tile_gdf).geometry.unary_union

        final_polygon = water_1.union(water_2)
        water_shp_buffered = gpd.GeoDataFrame(geometry=[final_polygon], crs=water_aux.crs)

        band_paths = [i for i in os.listdir(img_path) if i.endswith(".tif")]

        for band in band_paths:

            with rasterio.open(os.path.join(img_path, band)) as src:

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

                masked_band = np.where(water_raster == 1, band_data, -9999)
                final_image = np.where(masked_band < 0, -9999, masked_band)

                profile = src.profile
                profile.update(dtype=band_data.dtype, nodata=src.nodata)

                # Save clipped image
                if not os.path.exists(os.path.join(params['output_dir_wm'], os.path.basename(img_path))): os.makedirs(os.path.join(params['output_dir_wm'], os.path.basename(img_path)))

                output_path = os.path.join(params['output_dir_wm'], os.path.basename(img_path), os.path.basename(band))

                with rasterio.open(output_path, "w", **profile) as dst:
                    dst.write(final_image, 1)

    def run(self, params):

        """
        Runs a given set of parameters to initiate a selection process for satellite data.
        """

        sat = params['aux_info']['sat_name']

        if sat == 'landsat':
            sen_tiles = [params[sat]['tiles']]
        else:
            sen_tiles = [params['sentinel']['tiles']]

        for sen_tile_target in sen_tiles:

            params['sen_tile_target'] = sen_tile_target

            # Create the output directory if it does not exist
            params['output_dir_wm'] = fr'{params["output_dir"]}\water_mask'
            os.makedirs(fr'{params["output_dir_wm"]}', exist_ok=True)

            if sat == 'landsat':
                ncode = 'L30'
                path_pr = fr'{params["output_dir"]}\tiling\{params["sen_tile_target"]}\landsat'
                scenes = [fr'{path_pr}\{i}' for i in os.listdir(path_pr)]
            else:
                ncode = 'S30'
                path_pr = fr'{params["output_dir"]}\tiling\{params["sen_tile_target"]}\sentinel'
                scenes = [fr'{path_pr}\{i}' for i in os.listdir(path_pr)]

            params['ncode'] = ncode
            params["sat"] = sat

            n_params = [params] * len(scenes)

            # Process water mask extraction
            for scene in scenes:
                self.clip_water(scene, params)

            # with multiprocessing.Pool(processes=params['aux_info']['n_cores']) as pool:
            #     results = pool.starmap_async(self.clip_water, zip(scenes, n_params)).get()
            #     print(results)
            #     pool.close()