import os
import rasterio
import numpy as np
import rioxarray as rxr
import geopandas as gpd
from shapely.geometry import shape
import xml.etree.ElementTree as ET
from rasterio.features import shapes
from scipy.ndimage import binary_dilation as bn

def raster2shp(raster_numpy, raster_transform, raster_crs):

    results = list(
        {"properties": {"raster_val": v}, "geometry": s}
        for s, v in shapes(np.asarray(raster_numpy, dtype=np.int16), transform=raster_transform)
        if v  # Only take shapes with raster_val = True (i.e., v=1)
    )
    geometries = [shape(feature["geometry"]) for feature in results]

    gdf = gpd.GeoDataFrame(geometry=geometries, crs=raster_crs)

    # Keep only the polygons with the biggest area
    # gdf['area'] = gdf['geometry'].area
    # gdf_largest = gdf.sort_values(by='area', ascending=False).head(1)

    return gdf

def buffer_into_polygon(gdf, buffer_size):

    gdf['geometry'] = gdf.geometry.buffer(buffer_size)
    gdf_cleaned = gdf[~gdf['geometry'].is_empty & gdf['geometry'].is_valid]

    # Keep only the polygons with the biggest area
    # gdf_cleaned['area'] = gdf_cleaned['geometry'].area
    # gdf_largest = gdf_cleaned.sort_values(by='area', ascending=False).head(1)

    return gdf_cleaned


def create_inward_buffer(input_raster_path, buffer_distance_meters=3000, pixel_resolution=30):
    band_path = [i for i in os.listdir(input_raster_path)]
    green_band = next((band for band in band_path if "B3" in band or "B03" in band), None)
    swir_band = next((band for band in band_path if "B11" in band or "B6" in band or "B06" in band), None)

    xda_green = rxr.open_rasterio(os.path.join(input_raster_path, green_band))
    xda_swir = rxr.open_rasterio(os.path.join(input_raster_path, swir_band))

    xda_green_matched = xda_green.rio.reproject_match(xda_swir)

    mndwi = (xda_green_matched - xda_swir) / (xda_green_matched + xda_swir)
    mndwi_mask = mndwi > 0.2
    water_mask = mndwi_mask.astype(int).values[0, :, :]

    # Calculate buffer in pixels (buffer distance converted to pixels)
    buffer_pixels = int(abs(buffer_distance_meters) / pixel_resolution)

    # Invert the water mask (0 -> water, 1 -> non-water)
    inverted_mask = np.where(water_mask == 1, 0, 1)

    # Apply dilation to the inverted mask to expand the non-water area (expand the outer non-water body)
    dilated_mask = bn(inverted_mask, structure=np.ones((buffer_pixels, buffer_pixels)))

    # Invert back the dilated mask to create the inward buffer (1 for inside the buffer, 0 for outside)
    buffer_mask = np.where(np.logical_and(dilated_mask == 1, water_mask == 1), 1, 0)

    return

def get_mask_water(img_path):

    band_path = [i for i in os.listdir(img_path)]
    green_band = next((band for band in band_path if "B3" in band or "B03" in band), None)
    swir_band = next((band for band in band_path if "B11" in band or "B6" in band or "B06" in band), None)

    xda_green = rxr.open_rasterio(os.path.join(img_path, green_band))
    xda_swir = rxr.open_rasterio(os.path.join(img_path, swir_band))

    xda_green_matched = xda_green.rio.reproject_match(xda_swir)

    mndwi = (xda_green_matched - xda_swir) / (xda_green_matched + xda_swir)
    mndwi_mask = mndwi > 0.2
    water_mask = mndwi_mask

    with rasterio.open(os.path.join(img_path, green_band)) as src:
        _transform = src.transform
        _crs = src.crs

    # Include adjacency buffer
    water_shp = raster2shp(water_mask, _transform, _crs)
    water_copy = water_shp.copy()

    water_shp_buffered = buffer_into_polygon(water_copy, -3000)

    # Get only the difference between the original and the buffered
    water_shp_final = gpd.GeoDataFrame(geometry=water_shp.difference(water_shp_buffered), crs=water_shp.crs)

    # water_shp_final_copy = water_shp_final.copy()
    # water_adjc = buffer_into_polygon(water_shp_final_copy, 500)

    return water_shp_final

def create_grid(size, pixel_size):

    """
    Create a grid with a given size and pixel resolution.
    """
    center = size // 2
    x, y = np.meshgrid(np.arange(size), np.arange(size))
    distance_map = np.maximum(np.abs(x - center), np.abs(y - center)) * pixel_size
    return distance_map

def xml_to_dict(element):

    if len(element) == 0:
        return element.text

    result = {}

    for child in element:

        child_dict = xml_to_dict(child)

        if child.tag in result:
            if isinstance(result[child.tag], list):
                result[child.tag].append(child_dict)
            else:
                result[child.tag] = [result[child.tag], child_dict]
        else:
            result[child.tag] = child_dict

    return result




