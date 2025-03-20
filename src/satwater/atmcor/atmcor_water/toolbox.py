import os
import json
import pathlib
import datetime
import rasterio
import xmltodict
import numpy as np
import xarray as xr
import geopandas as gpd
import rioxarray as rxr
from shapely.geometry import shape
import xml.etree.ElementTree as ET
from rasterio.features import shapes

def xml_to_json(path_metadata: str):
    """
    Converts a file from .xml to .json.
    """
    with open(path_metadata) as xml_file:
        # Converts the .xml to dict(obj):
        file_dictionary = xmltodict.parse(xml_file.read())
        xml_file.close()
        # Converts the dict to .json:
        return json.loads(json.dumps(file_dictionary))

def export(array: float, index: str, reference: str, dest: str) -> None:

    """
    Exports a single-band raster to the specified destination using rioxarray.
    """

    #filename_out_factor = f"{dest}/{index[:-4]}.tif"
    filename_out_factor = dest

    # Open the reference raster to copy geospatial information
    reference_raster = rxr.open_rasterio(reference, masked=True)

    # Convert the array to an xarray DataArray and assign spatial metadata
    output_raster = xr.DataArray(
        array,
        dims=("y", "x"),
        coords={"y": reference_raster.y, "x": reference_raster.x},
        attrs=reference_raster.attrs  # Copy metadata including CRS
    )

    # Set CRS and transform info
    output_raster.rio.write_crs(reference_raster.rio.crs, inplace=True)
    output_raster.rio.write_transform(reference_raster.rio.transform(), inplace=True)

    # Export to GeoTIFF
    output_raster.rio.to_raster(filename_out_factor, dtype="float32", compress="lzw")

    print(f"Exported raster saved to {filename_out_factor}")

def newdirectory(path: str, name: str) -> str:
    """
    Creates a new directory in the specified path.
    """
    saved_path = path + '/' + name
    pathlib.Path(saved_path).mkdir(parents=True, exist_ok=True)
    return saved_path

def dict_to_xml(dictionary, parent=None):

    """
    Converts the dict into .xml.
    """

    if parent is None:
        parent = ET.Element('root')

    for key, value in dictionary.items():
        if isinstance(value, dict):
            dict_to_xml(value, ET.SubElement(parent, key))
        else:
            child = ET.SubElement(parent, key)
            child.text = str(value)
    return parent


def export_meta(meta, atmos_param, dest: str) -> None:

    """
    It writes and exports the metadata as a file .xml.
    """

    current_datetime = datetime.datetime.now()
    image_datetime = datetime.datetime(meta.datetime.year, meta.datetime.month, meta.datetime.day, int(meta.datetime.time_hh), int((meta.datetime.time_hh - int(meta.datetime.time_hh)) * 60))

    geometry = {'B' + str(key): value for key, value in meta.geometry.items()}
    rescale = {'B' + str(key): value for key, value in meta.rescale.items()}
    bands = {'B' + str(i): meta.bandname[i] for i in range(len(meta.bandname))}
    sixsv_params = {'B' + str(key): value for key, value in atmos_param.values_adjcorr.items()}

    metadata_dict = {
        'General_Info': {
            'software': 'GCERatmos',
            'version': '1',
            'satellite': meta.satellite,
            'type': meta.type,
            'datetime_image': image_datetime.isoformat(),
            'bandnumber': len(meta.bandname),
            'bandname': bands
        },
        'Path': {
            'original': meta.path_main,
            'dest': meta.path_dest
        },
        'InputData': {
            'aod': meta.aod,
            'water_vapour': meta.water_vapour,
            'ozone': meta.ozone,
            'altitude': meta.altitude,
            'geometry': geometry,
            'rescale': rescale,
            'sixSV_params': sixsv_params
        },
        'Outputdata': {
            'pixel_value': 'surface_reflectance',
            'NaN_value': '-9999',
            'file_format': 'tif',
            'datetime_processing': current_datetime.isoformat()
        }
    }

    root = ET.Element('root')
    for key, value in metadata_dict.items():
        if isinstance(value, dict):
            dict_to_xml(value, ET.SubElement(root, key))
        else:
            child = ET.SubElement(root, key)
            child.text = str(value)

    tree = ET.ElementTree(root)
    tree.write(dest + '/' + 'MTD.xml', encoding='utf-8', xml_declaration=True)
    return None

def jp2_to_tiff_xarray(path, dest):

    with rasterio.open(path) as src:

        out_meta = src.meta.copy()

        out_image = src.read(1)

        # Save the clipped image
        out_meta.update({"driver": "GTiff",
                         "height": out_image.shape[0],
                         "width": out_image.shape[1],
                         "transform": src.transform})

        with rasterio.open(dest, "w", **out_meta) as dest:
            dest.write(out_image, 1)

def return_water(image_path, rescale_factors, msi_tile):

    '''
    This function returns the water mask (as geodataframe) from a given image path.
    '''

    # Search for msi tile
    sentinel_tiles_path = r"C:\Users\tml411\Documents\Python Scripts\hls_water\src\satwater\auxfiles\tiles\MGRS_tiles.shp"
    msi_tiles_gdf = gpd.read_file(sentinel_tiles_path)

    tile_gdf = msi_tiles_gdf[msi_tiles_gdf['Name'] == msi_tile]

    band_path = [i for i in os.listdir(image_path)]
    green_band = next((band for band in band_path if "B3" in band or "B03" in band), None)
    swir_band = next((band for band in band_path if "B11" in band or "B6" in band or "B06" in band), None)

    xda_green = rxr.open_rasterio(os.path.join(image_path, green_band))
    xda_swir = rxr.open_rasterio(os.path.join(image_path, swir_band))

    xda_green_matched = xda_green.rio.reproject_match(xda_swir)

    # Convert to TOA reflectance
    try:
        xda_swir = xda_swir * rescale_factors[1]['qvalue'] + rescale_factors[1]['offset']
        xda_green_matched = xda_green_matched * rescale_factors[0]['qvalue'] + rescale_factors[0]['offset']
    except:
        xda_swir = xda_swir * rescale_factors[1]['mult'] + rescale_factors[1]['add']
        xda_green_matched = xda_green_matched * rescale_factors[0]['mult'] + rescale_factors[0]['add']

    mndwi = (xda_green_matched - xda_swir) / (xda_green_matched + xda_swir)
    mndwi_mask = mndwi > 0
    water_mask = mndwi_mask

    with rasterio.open(os.path.join(image_path, green_band)) as src:
        _transform = src.transform
        _crs = src.crs

    results = list(
        {"properties": {"raster_val": v}, "geometry": s}
        for s, v in shapes(np.asarray(water_mask, dtype=np.int16), transform=_transform)
        if v  # Only take shapes with raster_val = True (i.e., v=1)
    )

    geometries = [shape(feature['geometry']) for feature in results]

    gdf = gpd.GeoDataFrame(geometry=geometries, crs=_crs).to_crs(4326)

    # Return the gdf_largest polygons that intersect with the tile
    gdf_tile = gpd.overlay(gdf, tile_gdf, how='intersection')

    gdf_tile['area'] = gdf_tile['geometry'].area
    gdf_largest = gdf_tile.sort_values(by='area', ascending=False).head(1)

    return gdf_largest

def return_bbox(image_path):

    with rasterio.open(image_path) as _bandl:

         band_crs = _bandl.crs.to_epsg()
         image_data = _bandl.read(1)
         valid_data_mask = image_data != 0

    results = (
         {'properties': {'raster_val': v}, 'geometry': s}
         for s, v in shapes(valid_data_mask.astype(np.int16), transform=_bandl.transform)
         if v  # Only take shapes with raster_val = True (i.e., v=1)
    )

    geometries = [shape(feature['geometry']) for feature in results]

    gdf = gpd.GeoDataFrame(geometry=geometries, crs=band_crs).to_crs(4326)

    return gdf