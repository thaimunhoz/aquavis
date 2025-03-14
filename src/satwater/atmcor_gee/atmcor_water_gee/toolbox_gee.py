import os
import json
import pathlib
import datetime
import rasterio
import xmltodict
import numpy as np
from osgeo import gdal
import geopandas as gpd
import subprocess as sub
from netCDF4 import Dataset
from shapely.geometry import shape
import xml.etree.ElementTree as ET
from rasterio.features import shapes

gdal.UseExceptions()

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


def loadarray(path: str):
    """
    Loads a single band and returns an array.
    """
    # dataset = gdal.Open(path, GA_ReadOnly)
    dataset = gdal.Open(path)
    return dataset.ReadAsArray().astype(float)

def export(array: float, index: str, reference: str, dest: str) -> None:

    """
    Exports a single band to dest.
    """
    filename_reference = reference
    filename_out_factor = dest + '/' + index[0:-4] + '.tif'
    dataset_reference = gdal.Open(filename_reference)

    line = dataset_reference.RasterYSize
    column = dataset_reference.RasterXSize
    bands = 1

    # defining drive
    driver = gdal.GetDriverByName('GTiff')
    # copying the bands data type pre-existing
    data_type = gdal.GetDataTypeByName('Float32')
    # create new dataset
    dataset_output = driver.Create(filename_out_factor, column, line, bands, data_type)
    # copying the spatial information pre-existing
    dataset_output.SetGeoTransform(dataset_reference.GetGeoTransform())
    # copying the projection information pre-existing
    dataset_output.SetProjection(dataset_reference.GetProjectionRef())
    # writing array data in band
    dataset_output.GetRasterBand(1).WriteArray(array)
    dataset_output=None
    return None


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

def jp2_to_geotiff(path: str, dest: str) -> None:
    """
    Converts the .jp2 into .GeoTIFF.
    """
    image = gdal.Open(path)
    gdal.Translate(dest, image, format='GTiff')
    return None


def netcdf_to_geotiff(path_main: str, dest: str) -> None:
    """
    Converts .netcdf into .GeoTIFF --OLCI/Sentinel-3 images.
    """
    # Recovering the spectral bands and coordinates:
    BANDNAMES = ['Oa{0}_radiance'.format(str(i).zfill(2)) for i in range(1, 22)] # Including the bands from Oa1 to Oa21.
    nc_coords = os.path.join(path_main, 'geo_coordinates.nc')
    ds_nc = Dataset(nc_coords, 'r')
    coords = (ds_nc.variables['latitude'], ds_nc.variables['longitude'])
    # Getting the lat/long:
    lat_tif = os.path.join(dest, coords[0].name + '.tif')
    lon_tif = os.path.join(dest, coords[1].name + '.tif')
    rad_tifs = []
    nc_paths = [os.path.join(path_main, band + '.nc') for band in BANDNAMES]
    for v, var in enumerate(coords):
        nodata = var._FillValue
        scale = var.scale_factor
        var_vrt = os.path.join(dest, var.name + '.vrt')
        var_tif = os.path.join(dest, var.name + '.tif')
        # Creating the lat/long in the .VRT:
        cmd = ['gdalbuildvrt', '-sd', str(1 + var._varid), '-separate', '-overwrite', var_vrt, nc_coords]
        sub.call(cmd, stdout=sub.DEVNULL, stderr=sub.STDOUT)
        # Editing the .VRT from lat/long:
        # It is necessary to sure that the data are well dimensioned and located.
        with open(var_vrt, 'r') as f:
            xml = f.readlines()
        for line in xml:
            if '<VRTRasterBand ' in line:
                head_index = xml.index(line) + 1
            if '<DstRect xOff' in line:
                tail_index = xml.index(line) + 1
        xml.insert(head_index, '    <NoDataValue>{nd}</NoDataValue>\n'.format(nd=nodata))
        xml.insert(head_index + 1, '    <Scale>{sc}</Scale>\n'.format(sc=scale))
        tail_index = tail_index + 2
        xml.insert(tail_index, '      <NODATA>{nd}</NODATA>\n'.format(nd=nodata))
        xml.insert(tail_index + 2, '    <Offset>0.0</Offset>\n')
        xml.insert(tail_index + 3, '    <Scale>{sc}</Scale>\n'.format(sc=scale))
        xml = [line.replace('="Int32', '="Float32') for line in xml]
        with open(var_vrt, 'w') as f:
            f.writelines(xml)
        # Writing the lat/long .tifF --temporary:
        cmd = ['gdal_translate', '-unscale', var_vrt, var_tif]
        sub.call(cmd, stdout=sub.DEVNULL, stderr=sub.STDOUT)
    ds_nc.close()
    # Creating .VRT to bands:
    for n, nc in enumerate(nc_paths):
        ds_nc = Dataset(nc, 'r')
        var = ds_nc.variables[os.path.basename(nc)[:-3]]
        nodata = var._FillValue
        offset = var.add_offset
        rows = var.shape[0]
        scale = var.scale_factor
        ds_nc.close()
        data_vrt = os.path.join(dest, 'data.vrt')
        data_vrt_tif = data_vrt.replace('.vrt', '.tif')
        out_vrt = os.path.join(dest, os.path.basename(nc)[:-3] + '.vrt')
        out_tif = out_vrt.replace('.vrt', '.tif')
        cmd = ['gdalbuildvrt', '-sd', '1', '-separate', '-overwrite', data_vrt, nc]
        sub.call(cmd, stdout=sub.DEVNULL, stderr=sub.STDOUT)
        # Editing the band's .VRT:
        with open(data_vrt, 'r') as f:
            xml = f.readlines()
        for line in xml:
            if '<VRTRasterBand ' in line:
                head_index = xml.index(line)
            if '<DstRect xOff' in line:
                tail_index = xml.index(line) + 1
        xml[head_index] = '  <VRTRasterBand dataType="Float32" band="1">\n'
        xml.insert(head_index + 1, '    <NoDataValue>{nd}</NoDataValue>\n'.format(nd=nodata))
        xml[head_index + 2] = '    <ComplexSource>\n'
        xml[head_index + 5] = xml[head_index + 5].replace('DataType="UInt16"', 'DataType="Float32"')
        tail_index = tail_index + 1
        xml.insert(tail_index, '      <NODATA>{nd}</NODATA>\n'.format(nd=nodata))
        xml[tail_index + 1] = '    </ComplexSource>\n'
        xml.insert(tail_index + 2, '    <Offset>{off}</Offset>\n'.format(off=offset))
        xml.insert(tail_index + 3, '    <Scale>{sc}</Scale>\n'.format(sc=scale))
        with open(data_vrt, 'w') as f:
            f.writelines(xml)
        # Writing .tifF from band:
        cmd = ['gdal_translate', '-unscale', data_vrt, data_vrt_tif]
        sub.call(cmd, stdout=sub.DEVNULL, stderr=sub.STDOUT)
        # Updating the GeoTransform:
        ds = gdal.Open(data_vrt_tif, gdal.GA_Update)
        ds.SetGeoTransform((0.0, 1.0, 0.0, float(rows), 0.0, -1.0))
        ds.FlushCache()
        # Building the new band's .VRT -- inserting the geolocation information:
        cmd = ['gdalbuildvrt', '-sd', '1', '-separate', '-overwrite', out_vrt, data_vrt_tif]
        sub.call(cmd, stdout=sub.DEVNULL, stderr=sub.STDOUT)
        # Editing the new .VRT:
        with open(out_vrt, 'r') as f:
            xml = f.readlines()
        for line in xml:
            if '<VRTRasterBand ' in line:
                head_index = xml.index(line)
                break
        xml[head_index] = '  <VRTRasterBand dataType="Float32" band="1">\n'
        xml.insert(-1, '''  <metadata domain="GEOLOCATION">
            <mdi key="X_DATASET">{lon}</mdi>
            <mdi key="X_BAND">1</mdi>
            <mdi key="Y_DATASET">{lat}</mdi>
            <mdi key="Y_BAND">1</mdi>
            <mdi key="PIXEL_OFFSET">0</mdi>
            <mdi key="LINE_OFFSET">0</mdi>
            <mdi key="PIXEL_STEP">1</mdi>
            <mdi key="LINE_STEP">1</mdi>
          </metadata>\n'''.format(lon=lon_tif, lat=lat_tif))
        for line in xml:
            if os.sep in line:
                xml[xml.index(line)] = line.replace(os.sep, '/')
        with open(out_vrt, 'w') as f:
            f.writelines(xml)
        # Converting the .VRT in .tifF:
        cmd = ['gdalwarp', '-t_srs', 'epsg:4326', '-geoloc', '-srcnodata', str(nodata), '-dstnodata', '-9999', out_vrt,
               out_tif]
        sub.call(cmd, stdout=sub.DEVNULL, stderr=sub.STDOUT)
        # Removing the temp files safely:
        os.remove(out_vrt)
        ds = gdal.Open(data_vrt_tif, gdal.GA_ReadOnly)
        ds = None
        os.remove(data_vrt_tif)
        rad_tifs.append(out_tif)
    return None


def writeNetCDF(original_netcdf: str, blue_stack, red_stack, swir_stack, new_nc_file: str) -> None:
    # Open the original NetCDF file:
    with Dataset(original_netcdf, 'r') as ncfile:
        # Reads latitude and longitude:
        latitude = ncfile['geolocation_data']['latitude'][:]
        longitude = ncfile['geolocation_data']['longitude'][:]
        wavelength_blue = ncfile['sensor_band_parameters']['blue_wavelength'][:]
        wavelength_red = ncfile['sensor_band_parameters']['red_wavelength'][:]
        wavelength_swir = ncfile['sensor_band_parameters']['SWIR_wavelength'][:]

    # Create a new NetCDF file:
    with Dataset(new_nc_file, 'w', format='NETCDF4') as ncfile:
        # Create dimensions
        lat_dim = ncfile.createDimension('lat', latitude.shape[0])
        lon_dim = ncfile.createDimension('lon', latitude.shape[1])
        wavelength_blue_dim = ncfile.createDimension('wavelength_blue', wavelength_blue.shape[0])
        wavelength_red_dim = ncfile.createDimension('wavelength_red', wavelength_red.shape[0])
        wavelength_swir_dim = ncfile.createDimension('wavelength_swir', wavelength_swir.shape[0])

        # Create coordinate variables:
        latitudes_var = ncfile.createVariable('latitude', np.float32, ('lat', 'lon'))
        longitudes_var = ncfile.createVariable('longitude', np.float32, ('lat', 'lon'))
        wavelengths_blue_var = ncfile.createVariable('wavelength_blue', np.float32, ('wavelength_blue',))
        wavelengths_red_var = ncfile.createVariable('wavelength_red', np.float32, ('wavelength_red',))
        wavelengths_swir_var = ncfile.createVariable('wavelength_swir', np.float32, ('wavelength_swir',))

        # Assign data to coordinate variables:
        latitudes_var[:, :] = latitude
        longitudes_var[:, :] = longitude
        wavelengths_blue_var[:] = wavelength_blue
        wavelengths_red_var[:] = wavelength_red
        wavelengths_swir_var[:] = wavelength_swir

        # Create a new group:
        group = ncfile.createGroup('observation_data')

        # Create variables within the group for blue band data:
        data_blue_var = group.createVariable('band_blue', np.float32, ('wavelength_blue', 'lat', 'lon'))
        data_red_var = group.createVariable('band_red', np.float32, ('wavelength_red', 'lat', 'lon'))
        data_swir_var = group.createVariable('band_swir', np.float32, ('wavelength_swir', 'lat', 'lon'))

        # Assign blue | red | swir stack data to the variable:
        data_blue_var[:, :] = blue_stack  # Ensure blue_stack is correctly generated/retrieved
        data_red_var[:, :] = red_stack  # Ensure blue_stack is correctly generated/retrieved
        data_swir_var[:, :] = swir_stack  # Ensure blue_stack is correctly generated/retrieved

        # Add attributes to coordinate variables:
        latitudes_var.units = 'degrees_north'
        longitudes_var.units = 'degrees_east'

        # Add attributes to the data variable:
        data_blue_var.units = 'None'
        data_blue_var.long_name = 'Surface Reflectance Blue'

        data_red_var.units = 'None'
        data_red_var.long_name = 'Surface Reflectance Red'

        data_swir_var.units = 'None'
        data_swir_var.long_name = 'Surface Reflectance SWIR'
    return None

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