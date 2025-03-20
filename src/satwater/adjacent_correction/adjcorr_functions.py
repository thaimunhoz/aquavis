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


def atmospheric_point_scattering_function(atmosphere_parameters, band):

    """
    Calculates the Fr weight per distance.
    """

    # Converts the grid distance to km:
    if atmosphere_parameters["General_Info"]["satellite"] == "MSI_S2":
        if band == "B11":  # 20m
            grid_matrix = create_grid(30, 20)
        else:  # 10m
            grid_matrix = create_grid(60, 10)
    else:  # 30m
        grid_matrix = create_grid(20, 30)

    radius_km = grid_matrix / 1000

    # Zenith view angle -> degree:
    view_z = float(atmosphere_parameters["InputData"]['sixSV_params'][band]['view_z'])

    # Rayleigh UPWARD diffuse transmittance -> T_upward_difRayleigh:
    Rayleigh_OpticalDepth = float(atmosphere_parameters["InputData"]['sixSV_params'][band]['optical_depth__total_Ray'])
    T_upward_Rayleigh = float(atmosphere_parameters["InputData"]['sixSV_params'][band]['rayleigh_scatransmi_upward'])
    T_upward_dirRayleigh = np.exp(-Rayleigh_OpticalDepth / np.cos(view_z * (np.pi / 180)))
    T_upward_difRayleigh = T_upward_Rayleigh - T_upward_dirRayleigh

    # Aerosol UPWARD diffuse transmittance -> T_upward_difAerosol:
    Aerosol_OpticalDepth = float(atmosphere_parameters["InputData"]['sixSV_params'][band]['optical_depth__total_Aero'])
    T_upward_Aerosol = float(atmosphere_parameters["InputData"]['sixSV_params'][band]['aerosol_scatransmi_upward'])
    T_upward_dirAerosol = np.exp(-Aerosol_OpticalDepth / np.cos(view_z * (np.pi / 180)))
    T_upward_difAerosol = T_upward_Aerosol - T_upward_dirAerosol

    # Calculates the Aerosol's Fr and Rayleigh's Fr functions using the equation described by Vermote et al.(2006):
    FrRayleigh = ((0.930 * np.exp(-0.08 * radius_km)) + (0.070 * np.exp(-1.10 * radius_km)))
    FrAerosol = ((0.448 * np.exp(-0.27 * radius_km)) + (0.552 * np.exp(-2.83 * radius_km)))

    # Calculates the APSF (Fr) -> Atmospheric Point Scattering Function:
    Fr = (T_upward_difRayleigh * FrRayleigh + T_upward_difAerosol * FrAerosol) / (
                T_upward_difRayleigh + T_upward_difAerosol)

    return Fr


def Adjacency_correction(atmosphere_parameters, band, array_band, adjc_array):

    """
    Removes the adjacency effect of the image, using the equation described in the Vermote et al. (1997).
    """

    # Zenith view angle -> degree:
    view_z = float(atmosphere_parameters["InputData"]['sixSV_params'][band]['view_z'])

    # Atmopheric optical depth (Rayleigh + Aerosol) -> Atmospheric_OpticalDepth:
    Atmospheric_OpticalDepth = float(
        atmosphere_parameters["InputData"]['sixSV_params'][band]['optical_depth__total_AeroRay'])

    # Total transmittance UPWARD (Rayleigh + Aerosol) -> T_upward:
    T_upward = float(atmosphere_parameters["InputData"]['sixSV_params'][band]['total_scattering_transmittance_upward'])

    # Total transmittance UPWARD direct (Rayleigh + Aerosol) -> T_upward_dirAeroRay:
    T_upward_dirAeroRay = np.exp(-Atmospheric_OpticalDepth / np.cos(view_z * (np.pi / 180)))

    # Total transmittance UPWARD diffuse (Rayleigh + Aerosol) -> T_upward_difAeroRay
    T_upward_difAeroRay = T_upward - T_upward_dirAeroRay

    # Surface reflectance without adjacency effect - Vermote et al. (1997):
    sr_corr = array_band - (adjc_array * T_upward_difAeroRay)

    return sr_corr