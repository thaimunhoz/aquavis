import os
import json
import pathlib
import datetime
import rasterio
import xmltodict
import numpy as np
import geopandas as gpd
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


def export_meta(meta, dest: str) -> None:

    """
    It writes and exports the metadata as a file .xml.
    """

    current_datetime = datetime.datetime.now()
    image_datetime = datetime.datetime(meta.datetime.year, meta.datetime.month, meta.datetime.day, int(meta.datetime.time_hh), int((meta.datetime.time_hh - int(meta.datetime.time_hh)) * 60))

    geometry = {'B' + str(key): value for key, value in meta.geometry.items()}
    rescale = {'B' + str(key): value for key, value in meta.rescale.items()}
    bands = {'B' + str(i): meta.bandname[i] for i in range(len(meta.bandname))}

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
            'rescale': rescale
        },
        'Outputdata': {
            'pixel_value': 'surface_reflectance',
            'NaN_value': '-9999',
            'file_format': 'TIFF',
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