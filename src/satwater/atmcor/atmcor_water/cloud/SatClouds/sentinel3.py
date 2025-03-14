# -*- mode: python -*- Rejane Paulino 06-05-2024

"""
SatClouds - Sentinel-3/OLCI.
"""


import os
from osgeo import gdal
gdal.UseExceptions()
from netCDF4 import Dataset
import pathlib
import shutil
import subprocess as sub
import glob
import numpy as np


class Sentinel3:

    def __init__(self, path_img, dest):

        self.path_img = path_img
        self.dest = dest # directory with images (*.tifF).

    def run(self):
        """
        Runs and applies the cloud mask derived from qualityFlags.nc in S3X_OL_1_EFR____*.SEN3.
        """
        # Creates a new directory:
        pathxmain = self.newdirectory(self.dest, 'SatClouds')
        pathxtempdir = self.newdirectory(pathxmain, 'temp')
        # Exports the qualityFlags.nc as .tifF:
        self.netcdf_to_geotiff(self.path_img, pathxtempdir)
        # Filters the cloud mask and glint risk:
        arr = self.loadarray(pathxtempdir + '/qualityFlags.tif', np.uint32)
        cloud_mask = np.where(arr & 0x08000000, 0, 1) # 0 == cloud
        # sunglint_risk = np.where(arr & 0x00400000, 0, 1) # 0 == sunglint risk # It is based on Cox & Munk model, but it filters a large part from images.
        # Applies the cloud mask:
        paths = [band for band in glob.glob(os.path.join(self.dest, '*.tif'))]
        paths = sorted(paths)
        BANDNAMES = ['Oa{0}_reflectance.tif'.format(str(i).zfill(2)) for i in range(1, 22)]  # Including the bands from Oa1 to Oa21.
        for path, id in zip(paths, BANDNAMES):
            arr = self.loadarray(path, float)
            masked = arr * cloud_mask
            masked = np.where(masked == 0, -9999, masked)
            self.export(masked, id, path, pathxmain)
        # Removes the tempdir:
        shutil.rmtree(pathxtempdir)


    def newdirectory(self, path: str, name: str) -> str:
        """
        Creates a new directory in the specified path.
        """
        saved_path = path + '/' + name
        pathlib.Path(saved_path).mkdir(parents=True, exist_ok=True)
        return saved_path


    def netcdf_to_geotiff(self, path_main: str, dest: str) -> None:
        """
        Converts .netcdf into .GeoTIFF --OLCI/Sentinel-3 images.
        """
        # Recovering the spectral bands and coordinates:
        BANDNAMES = ['qualityFlags'] # Including the qualityFlags.nc
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
        reference = Dataset(path_main + '/Oa01_radiance.nc', 'r')
        var_reference = reference.variables['Oa01_radiance']
        nc = nc_paths[0]
        ds_nc = Dataset(nc, 'r')
        var = ds_nc.variables['quality_flags']
        nodata = var._FillValue if '_FillValue' in var.ncattrs() else None
        offset = var_reference.add_offset
        rows = var.shape[0]
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
        # xml.insert(tail_index + 3, '    <Scale>{sc}</Scale>\n'.format(sc=scale))
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

    def loadarray(self, path: str, type):
        """
        Loads a single band and returns an array.
        """
        # dataset = gdal.Open(path, GA_ReadOnly)
        dataset = gdal.Open(path)
        return dataset.ReadAsArray().astype(type)

    def export(self, array: float, index: str, reference: str, dest: str) -> None:
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
