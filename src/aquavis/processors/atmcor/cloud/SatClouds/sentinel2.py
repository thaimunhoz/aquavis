# -*- mode: python -*- Rejane Paulino 06-05-2024

"""
SatClouds - Sentinel-2/MSI.
"""

import os
from osgeo import gdal
gdal.UseExceptions()
import pathlib
import shutil
import glob
import numpy as np
import sys

os.environ['GDAL_DATA'] = 'C:/OSGeo4W/bin'
os.environ['PATH'] = 'C:/OSGeo4W/bin' + os.pathsep + os.environ['PATH']

class Sentinel2:

    def __init__(self, fmask_env, path_img, dest):

        self.fmask_env = fmask_env # conda env */bin/python.exe
        self.path_img = path_img # path from safe dir
        self.dest = dest # directory with images (*.tifF)

        self.SHADOW_BUFFER_DISTANCE = 300
        self.CLOUD_BUFFER_DISTANCE = 600
        self.PIXEL_SIZE = 60

    def run(self):

        """
        Runs and applies the cloud mask derived from Fmask-algorithm and S2X_MSIL1C_*.SAFE.
        """

        # Creates a new directory:
        pathxmain = self.newdirectory(self.dest, 'SatClouds')
        pathxtempdir = self.newdirectory(pathxmain, 'temp')

        # Generates the cloud mask based on Fmask:
        fmask_mode = r"C:\FMASK_2\python-fmask-master\bin\fmask_sentinel2Stacked.py"
        fmask_env = r'C:\Users\tml411\AppData\Local\anaconda3\envs\fmask_env\python'

        cmd = rf"{fmask_env} {fmask_mode} --shadowbufferdistance {self.SHADOW_BUFFER_DISTANCE} --cloudbufferdistance " \
              rf"{self.CLOUD_BUFFER_DISTANCE} --pixsize {self.PIXEL_SIZE} -e {pathxtempdir} -o {pathxtempdir + '/cloud.tif'} --safedir {self.path_img} --tempdir {pathxtempdir}"
        os.system(cmd)

        # Resamples the
        # Applies and filters the cloud mask:
        '''arr60 = self.loadarray(pathxtempdir + '/cloud.tif')
        bit_values = {
            0: 'Null',
            1: 'Clear land',
            2: 'Cloud',
            3: 'Shadow',
            4: 'Snow',
            5: 'Clear water'
        }
        self.resample(pathxtempdir + '/cloud.tif', 20, pathxtempdir)
        self.resample(pathxtempdir + '/cloud.tif', 10, pathxtempdir)
        paths = [band for band in glob.glob(os.path.join(self.dest, '*.tif'))]
        print(paths)
        for path in paths:
            if '_B02' in path or '_B03' in path or '_B04' in path or '_B08' in path:
                # For 10-meters:
                arr10 = self.loadarray(pathxtempdir + '/cloud10.tif')
                cloud_mask = np.where(arr10 == 2, 0, 1)
                shadow_mask = np.where(arr10 == 3, 0, 1)
                arr = self.loadarray(path)
                masked = arr * cloud_mask * shadow_mask
                masked_out = np.where(masked == 0, -9999, masked)
                self.export(masked_out, path[-30:], path, pathxmain)
            elif '_B05' in path or '_B06' in path or '_B07' in path or '_B8A' in path or '_B11' in path or '_B12' in path:
                # For 20-meters:
                arr20 = self.loadarray(pathxtempdir + '/cloud20.tif')
                cloud_mask = np.where(arr20 == 2, 0, 1)
                shadow_mask = np.where(arr20 == 3, 0, 1)
                arr = self.loadarray(path)
                masked = arr * cloud_mask * shadow_mask
                masked_out = np.where(masked == 0, -9999, masked)
                self.export(masked_out, path[-30:], path, pathxmain)
            else:
                # For 60-meters:
                cloud_mask = np.where(arr60 == 2, 0, 1)
                shadow_mask = np.where(arr60 == 3, 0, 1)
                arr = self.loadarray(path)
                masked = arr * cloud_mask * shadow_mask
                masked_out = np.where(masked == 0, -9999, masked)
                self.export(masked_out, path[-30:], path, pathxmain)'''
        # Removes the tempdir:
        #shutil.rmtree(pathxtempdir)


    def newdirectory(self, path: str, name: str) -> str:
        """
        Creates a new directory in the specified path.
        """
        saved_path = path + '/' + name
        pathlib.Path(saved_path).mkdir(parents=True, exist_ok=True)
        return saved_path

    def loadarray(self, path: str):
        """
        Loads a single band and returns an array.
        """
        # dataset = gdal.Open(path, GA_ReadOnly)
        dataset = gdal.Open(path)
        return dataset.ReadAsArray().astype(float)

    def resample(self, path: str, size: int, dest: str):
        """
        It resamples the pixel size.
        """
        gdal.Warp(dest + '/cloud' + str(size) + '.tif', path, xRes=size, yRes=size, resampleAlg='near')
        return None

    def export(self, array, index: str, reference: str, dest: str) -> None:
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

# Runs the Sentine-2 procedure:
# Extracts args from command line:
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Error: python sentinel2.py <fmask_env> <path_img> <dest>")
        sys.exit(1)

    fmask_env = sys.argv[1]
    path_img = sys.argv[2]
    dest = sys.argv[3]

    sentinel2 = Sentinel2(fmask_env, path_img, dest)
    sentinel2.run()
