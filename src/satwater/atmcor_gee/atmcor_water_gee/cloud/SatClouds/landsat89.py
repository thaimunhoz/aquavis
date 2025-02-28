# -*- mode: python -*- Rejane Paulino 06-05-2024

"""
SatClouds - Landsat89/OLI.
"""

import os
from osgeo import gdal
gdal.UseExceptions()
import pathlib
import shutil
import glob
import numpy as np
import sys


class Landsat89:

    def __init__(self, fmask_env, path_img, dest):

        self.fmask_env = fmask_env # conda env */bin/python.exe
        self.path_img = path_img # path from safe dir
        self.dest = dest # directory with images (*.TIFF)

        self.SHADOW_BUFFER_DISTANCE = 300
        self.CLOUD_BUFFER_DISTANCE = 600

    def run(self):
        """
        Runs and applies the cloud mask derived from Fmask-algorithm and LC09_L1TP_*.
        """
        # Creates a new directory:
        pathxmain = self.newdirectory(self.dest, 'SatClouds')
        pathxtempdir = self.newdirectory(pathxmain, 'temp')
        # Generates the cloud mask based on Fmask:
        fmask_mode = r"cloud/SatClouds/src/python-fmask-master/bin/fmask_usgsLandsatStacked.py"
        cmd = rf"{self.fmask_env} {fmask_mode} --shadowbufferdistance {self.SHADOW_BUFFER_DISTANCE} --cloudbufferdistance " \
              rf"{self.CLOUD_BUFFER_DISTANCE} -e {pathxtempdir} -o {pathxtempdir + '/cloud.tif'} --scenedir {self.path_img} --tempdir {pathxtempdir}"
        os.system(cmd)

        # Resamples the
        # Applies and filters the cloud mask:
        '''arr30 = self.loadarray(pathxtempdir + '/cloud.tif')
        bit_values = {
            0: 'Null',
            1: 'Clear land',
            2: 'Cloud',
            3: 'Shadow',
            4: 'Snow',
            5: 'Clear water'
        }
        self.resample(pathxtempdir + '/cloud.tif', 15, pathxtempdir)
        paths = [band for band in glob.glob(os.path.join(self.dest, '*.TIF'))]
        for path in paths:
            if '_B8' in path:
                # For 15-meters:
                arr15 = self.loadarray(pathxtempdir + '/cloud15.tif')
                cloud_mask = np.where(arr15 == 2, 1, 0)
                shadow_mask = np.where(arr15 == 3, 1, 0)
                arr = self.loadarray(path)
                masked = arr * cloud_mask * shadow_mask
                masked_out = np.where(masked == 0, -9999, masked)
                self.export(masked_out, path[-43:], path, pathxmain)
            else:
                # For 30-meters:
                cloud_mask = np.where(arr30 == 2, 1, 0)
                shadow_mask = np.where(arr30 == 3, 1, 0)
                arr = self.loadarray(path)
                masked = arr * cloud_mask * shadow_mask
                masked_out = np.where(masked == 0, -9999, masked)
                self.export(masked_out, path[-43:], path, pathxmain)'''
        # Removes the tempdir:
        shutil.rmtree(pathxtempdir)


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

    def export(self, array: float, index: str, reference: str, dest: str) -> None:
        """
        Exports a single band to dest.
        """
        filename_reference = reference
        filename_out_factor = dest + '/' + index[0:-4] + '.TIF'
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

# Runs the Landsat-8/9 procedure:
# Extracts args from command line:
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Error: python landsat89.py <fmask_env> <path_img> <dest>")
        sys.exit(1)

    fmask_env = sys.argv[1]
    path_img = sys.argv[2]
    dest = sys.argv[3]

    landsat89 = Landsat89(fmask_env, path_img, dest)
    landsat89.run()
