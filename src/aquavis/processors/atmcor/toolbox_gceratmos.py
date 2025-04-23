import os
import glob
import pathlib
import rasterio
import numpy as np
import xarray as xr
import rioxarray as rxr
from typing import Optional
from datetime import datetime
import xml.etree.ElementTree as ET

from src.aquavis.processors.data_class import AquaVisDataLoader
from src.aquavis.processors.atmcor.atm.atmosphere import Atmosphere
from src.aquavis.processors.atmcor.atm.correction import Correction

""" Pack of function used inside the atm_aquavis.py """

def _perform_atmospheric_correction(path_main: str, dest_dir: str, temp_dir: Optional[str] = None) -> None:
    """Perform atmospheric correction on all bands."""

    loader = AquaVisDataLoader()
    params = loader.load_aquavis_data()

    Atmosphere().run()

    for index, band in enumerate(params.bandname):

        if temp_dir is not None:
            input_path = os.path.join(temp_dir, f"{band[:-4]}.tif")
        else:
            input_path = os.path.join(path_main, f"{band[:-4]}.tif")

        output_path = os.path.join(dest_dir, f"{band[:-4]}.tif")

        arr = rxr.open_rasterio(input_path).squeeze().values.astype(float)
        corr = Correction(arr, index)
        corr.run()

        arr_corrected = np.where(corr.arr_sr < 0, -9999, corr.arr_sr)
        export(arr_corrected, input_path, output_path)

    params = loader.load_aquavis_data()
    export_meta(params, dest_dir)

def _convert_jp2_to_tiff(path_main: str, temp_dir: str, granule_dir: Optional[str] = '/GRANULE', img_data_dir: Optional[str] = '/IMG_DATA') -> None:
    """Convert JP2 files to TIFF format."""
    granule_path = glob.glob(os.path.join(path_main + granule_dir, '*L1C_*'))[0]
    band_jp2_files = glob.glob(os.path.join(granule_path + img_data_dir, '*_B*.jp2'))

    for jp2_file in band_jp2_files:
        output_name = os.path.join(temp_dir, f"{os.path.basename(jp2_file)[:-4]}.tif")
        jp2_to_tiff_xarray(jp2_file, output_name)

def _parse_date_from_filename(filename: str, satellite_type: str) -> datetime:
    """Extract and parse date from satellite image filename."""
    basename = os.path.basename(filename)

    if satellite_type.upper() == 'SENTINEL':
        date_str = basename.split('_')[2].split('T')[0]
    else:  # Landsat
        date_str = basename.split('_')[3]

    return datetime.strptime(date_str, '%Y%m%d')

def export(array: np.array, reference: str, dest: str) -> None:

    """Exports a single-band raster to the specified destination using rioxarray."""

    filename_out_factor = dest
    reference_raster = rxr.open_rasterio(reference, masked=True)

    output_raster = xr.DataArray(
        array,
        dims=("y", "x"),
        coords={"y": reference_raster.y, "x": reference_raster.x},
        attrs=reference_raster.attrs  # Copy metadata including CRS
    )

    output_raster.rio.write_crs(reference_raster.rio.crs, inplace=True)
    output_raster.rio.write_transform(reference_raster.rio.transform(), inplace=True)
    output_raster.rio.to_raster(filename_out_factor, dtype="float32", compress="lzw")

def newdirectory(path: str, name: str) -> str:
    """Creates a new directory in the specified path."""
    saved_path = path + '/' + name
    pathlib.Path(saved_path).mkdir(parents=True, exist_ok=True)
    return saved_path

def dict_to_xml(dictionary, parent=None):

    """Converts the dict into .xml."""

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
    """It writes and exports the metadata as a file .xml."""

    current_datetime = datetime.now()
    image_datetime = datetime(meta.datetime.year, meta.datetime.month, meta.datetime.day, int(meta.datetime.time_hh), int((meta.datetime.time_hh - int(meta.datetime.time_hh)) * 60))

    geometry = {'B' + str(key): value for key, value in meta.geometry.items()}
    rescale = {'B' + str(key): value for key, value in meta.rescale.items()}
    bands = {'B' + str(i): meta.bandname[i] for i in range(len(meta.bandname))}
    sixsv_params = {'B' + str(key): value for key, value in meta.values_adjcorr.items()}

    metadata_dict = {
        'General_Info': {
            'software': 'GCERatmos',
            'version': '1',
            'satellite': meta.sensor,
            'type': meta.type,
            'datetime_image': image_datetime.isoformat(),
            'bandnumber': len(meta.bandname),
            'bandname': bands,
            'atmcor_folder': os.path.dirname(dest)
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
    """Convert JP2 file to TIFF using xarray and rasterio."""
    with rasterio.open(path) as src:
        out_meta = src.meta.copy()
        out_image = src.read(1)
        out_meta.update({"driver": "GTiff",
                         "height": out_image.shape[0],
                         "width": out_image.shape[1],
                         "transform": src.transform})
        with rasterio.open(dest, "w", **out_meta) as dest:
            dest.write(out_image, 1)