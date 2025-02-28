import os
import glob
import shutil
import warnings
import rasterio
import numpy as np

from src.satwater.input_gee import toa_gee
from src.satwater.atmcor_gee.gceratmos_hls import toolbox_gee as tool
from src.satwater.atmcor_gee.gceratmos_hls.metadata_gee import Metadata_MSI_S2
from src.satwater.atmcor_gee.gceratmos_hls.metadata_gee import Metadata_OLI_L89
from src.satwater.atmcor_gee.gceratmos_hls.atm.atmosphere_gee import Atmosphere
from src.satwater.atmcor_gee.gceratmos_hls.atm.correction_gee import Correction

class Gceratmos_gee:

    def __init__(self,
                 path_main: str,
                 path_dest: str,
                 networkdrive_letter: str,
                 satellite: str,
                 aero_type: str,
                 mode=None,
                 path_buffer=None):

        self.GRANULE = '/GRANULE'
        self.IMG_DATA = '/IMG_DATA'

        self.path_main = path_main
        self.path_dest = path_dest
        self.networkdrive_letter = networkdrive_letter
        self.satellite = satellite
        self.aero_type = aero_type
        self.mode = mode
        self.path_buffer = path_buffer

        '''
        path_main (str): Image ID
        path_dest (str): Output directory
        '''

    def run(self):

        """
        Returns the surface reflectance.
        """

        if self.satellite == 'MSI_S2':

            print(self.path_main)

            # Access image from GEE and metadata:
            vnir_xda_10, vnir_xda_20, metadata_xda = toa_gee.sentinel_toa_from_gee(self.path_main)

            # Metadata:
            meta = Metadata_MSI_S2(self.path_main, self.path_dest, self.networkdrive_letter, self.satellite, self.aero_type, vnir_xda_10, metadata_xda, self.mode)
            meta.run()

            # Atmospheric parameters:
            atmos_param = Atmosphere(meta)
            atmos_param.run()

            # Atmospheric correction:
            vnir_xda_10 = vnir_xda_10.rename({'Y': 'y', 'X': 'x'})
            vnir_xda_20 = vnir_xda_20.rename({'Y': 'y', 'X': 'x'})

            bands_index = [1, 2, 3, 8, 11]
            bands_name = ['B2', 'B3', 'B4', 'B8A', 'B11']

            for index, band in enumerate(bands_name):

                if band == 'B8A' or band == 'B11':
                    vnir_xda = vnir_xda_20
                else:
                    vnir_xda = vnir_xda_10

                # Access the band data directly from vnir_xda
                arr = (vnir_xda[band].values)*0.0001  # Extract the band array as a NumPy array
                out_transform = vnir_xda[band].rio.transform()  # Get the geotransform information

                # Perform the correction
                corr = Correction(meta, atmos_param, arr, bands_index[index])
                corr.run()

                # Handle NaN values
                arr_c = np.where(corr.arr_sr < 0, -9999, corr.arr_sr)  # Replace negative values with -9999

                # Prepare metadata for output
                out_meta = {
                    "driver": "GTiff",
                    "height": arr_c.shape[1],
                    "width": arr_c.shape[2],
                    "count": 1,
                    "transform": out_transform,
                    "crs": vnir_xda[band].rio.crs,  # Coordinate Reference System
                    "dtype": 'float64',
                    "compress": 'lzw'
                }

                # Write the corrected array to a new GeoTIFF
                output_path = tool.newdirectory(self.path_dest, self.path_main) + '/' + self.path_main + '_' + band + '.TIF'
                with rasterio.open(output_path, "w", **out_meta) as dest:
                    dest.write(arr_c)

            tool.export_meta(meta, tool.newdirectory(self.path_dest, self.path_main))

        elif self.satellite == 'OLI_L8/9':

            print(self.path_main)

            # Access image from GEE and metadata:
            vnir_xda, angles_xda, metadata_xda = toa_gee.landsat_toa_from_gee(self.path_main)

            # Metadata:
            meta = Metadata_OLI_L89(self.path_main, self.path_dest, self.networkdrive_letter, self.satellite, self.aero_type, vnir_xda, angles_xda, metadata_xda, self.mode)
            meta.run()

            # Atmospheric parameters:
            atmos_param = Atmosphere(meta)
            atmos_param.run()

            # Atmospheric correction:
            vnir_xda = vnir_xda.rename({'Y': 'y', 'X': 'x'})
            for index, band in enumerate(meta.bandname):

                index = index + 1

                # Access the band data directly from vnir_xda
                arr = vnir_xda[band].values  # Extract the band array as a NumPy array
                out_transform = vnir_xda[band].rio.transform()  # Get the geotransform information

                # Perform the correction
                corr = Correction(meta, atmos_param, arr, index)
                corr.run()

                # Handle NaN values
                arr_c = np.where(corr.arr_sr < 0, -9999, corr.arr_sr)  # Replace negative values with -9999

                # Prepare metadata for output
                out_meta = {
                    "driver": "GTiff",
                    "height": arr_c.shape[1],
                    "width": arr_c.shape[2],
                    "count": 1,
                    "transform": out_transform,
                    "crs": vnir_xda[band].rio.crs,  # Coordinate Reference System
                    "dtype": 'float64',
                    "compress": 'lzw'
                }

                # Write the corrected array to a new GeoTIFF
                output_path = tool.newdirectory(self.path_dest, self.path_main[-40:]) + '/' + self.path_main[-40:] + '_' + band + '.TIF'
                with rasterio.open(output_path, "w", **out_meta) as dest:
                    dest.write(arr_c)

            tool.export_meta(meta, tool.newdirectory(self.path_dest, self.path_main[-40:]))

        else:

            warnings.warn("The sensor type was not identified in GCERATMOS.\n", UserWarning)