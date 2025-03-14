import os
import glob
import warnings
import rasterio
import numpy as np
import rioxarray as rxr

from src.satwater.input_gee import toa_gee
from src.satwater.atmcor_gee.atmcor_water_gee import toolbox_gee as tool
from src.satwater.atmcor_gee.atmcor_water_gee.metadata_gee import Metadata_MSI_S2
from src.satwater.atmcor_gee.atmcor_water_gee.metadata_gee import Metadata_OLI_L89
from src.satwater.atmcor_gee.atmcor_water_gee.atm.atmosphere_gee import Atmosphere
from src.satwater.atmcor_gee.atmcor_water_gee.atm.correction_gee import Correction

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

            aux_folder = tool.newdirectory(self.path_dest, 'tempdir2')

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
                with rasterio.open(aux_folder + '/' + self.path_main + '_' + band + '.TIF', "w", **out_meta) as dest:
                    dest.write(arr_c)

            # Glint correction based SWIR subtraction:
            swir_band = [i for i in glob.glob(os.path.join(aux_folder, '*B11.tif'))]
            xda_swir = rxr.open_rasterio(swir_band[0])
            xda_swir = xda_swir.where(xda_swir >= 0, 0)

            for band in bands_name:

                output_path = tool.newdirectory(self.path_dest, self.path_main) + '/' + self.path_main + '_' + band + '.tif'

                if 'B11' in band:
                    xda_band = rxr.open_rasterio(aux_folder + '/' + self.path_main + '_' + band + '.tif')
                    xda_band = xda_band.where(xda_band >= 0, 0)
                    xda_band.rio.to_raster(output_path)
                else:
                    xda_band = rxr.open_rasterio(aux_folder + '/' + self.path_main + '_' + band + '.tif')
                    swir_match = xda_swir.rio.reproject_match(xda_band)
                    glint_corr = np.where((xda_band - swir_match) <= 0, xda_band, (xda_band - swir_match))

                    glint_corr_data = xr.DataArray(glint_corr, dims=xda_band.dims, coords=xda_band.coords, attrs=xda_band.attrs)

                    glint_corr_data.rio.to_raster(output_path)

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

            bands_order = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B8A', 'B9', 'B10', 'B11', 'B12']
            bandnames = sorted(
                meta.bandname,
                key=lambda x: next((i for i, band in enumerate(bands_order) if band in x), float('inf'))
            )

            # Atmospheric correction:
            aux_folder = tool.newdirectory(self.path_main, 'tempdir2')
            vnir_xda = vnir_xda.rename({'Y': 'y', 'X': 'x'})
            for index, band in enumerate(bandnames):

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
                with rasterio.open(tool.newdirectory(aux_folder, self.path_main[-40:]) + '/' + band[0:-4] + '.TIF', "w",**out_meta) as dest:
                    dest.write(arr_c, 1)

            # Glint correction based SWIR subtraction:
            swir_band = [i for i in glob.glob(os.path.join(aux_folder, '*B11.tif'))]
            xda_swir = rxr.open_rasterio(swir_band[0])
            xda_swir = xda_swir.where(xda_swir >= 0, 0)

            for band in bandnames:

                output_path = tool.newdirectory(self.path_dest, self.path_main) + '/' + self.path_main + '_' + band + '.tif'

                if 'B11' in band:
                    xda_band = rxr.open_rasterio(aux_folder + '/' + self.path_main + '_' + band + '.tif')
                    xda_band = xda_band.where(xda_band >= 0, 0)
                    xda_band.rio.to_raster(output_path)
                else:
                    xda_band = rxr.open_rasterio(aux_folder + '/' + self.path_main + '_' + band + '.tif')
                    swir_match = xda_swir.rio.reproject_match(xda_band)
                    glint_corr = np.where((xda_band - swir_match) <= 0, xda_band, (xda_band - swir_match))

                    glint_corr_data = xr.DataArray(glint_corr, dims=xda_band.dims, coords=xda_band.coords,
                                                           attrs=xda_band.attrs)

                    glint_corr_data.rio.to_raster(output_path)

            tool.export_meta(meta, tool.newdirectory(self.path_dest, self.path_main[-40:]))

        else:

            warnings.warn("The sensor type was not identified in GCERATMOS.\n", UserWarning)