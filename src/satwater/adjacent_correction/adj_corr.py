import os
import glob
import shutil
import rasterio
import numpy as np
from rasterio.mask import mask
import xml.etree.ElementTree as ET
from scipy.signal import convolve2d
from src.satwater.utils import satwutils
from scipy.signal import fftconvolve
from src.satwater.adjacent_correction import adjcorr_functions as adjc_fun

class AdjCorrClass:

    def __int__(self):

        pass

    def fast_convolution(self, image, kernel):

        # Handle NaNs
        image = np.nan_to_num(image, nan=0)
        kernel = np.nan_to_num(kernel, nan=0)

        # Get original image and kernel sizes
        image_h, image_w = image.shape
        kernel_h, kernel_w = kernel.shape

        # Compute necessary padding (same as 'fill' boundary in convolve2d)
        pad_h = kernel_h // 2
        pad_w = kernel_w // 2

        # Pad image with zeros (replicates 'boundary=fill', fillvalue=0)
        padded_image = np.pad(image, ((pad_h, pad_h), (pad_w, pad_w)), mode='constant', constant_values=0)

        # Perform convolution in the frequency domain
        convolved = fftconvolve(padded_image, kernel, mode='same')

        # Crop the convolution result to match original image size
        convolved_cropped = convolved[pad_h:pad_h + image_h, pad_w:pad_w + image_w]

        return convolved_cropped

    def apply_adjcorr(self, input_path, output_path, params):

        self.tree = ET.parse(os.path.join(input_path, "MTD.xml"))
        self.output_path = output_path
        self.image_path = input_path

        root = self.tree.getroot()
        metadata_6sv = adjc_fun.xml_to_dict(root)

        params['aux_info']['date_time_info'] = metadata_6sv["General_Info"]["datetime_image"]

        shutil.copy(os.path.join(input_path, "MTD.xml"), os.path.join(output_path, "MTD.xml"))

        # Selecting only those important band for AQUAVis
        if metadata_6sv["General_Info"]["satellite"] == "MSI_S2":

            aquavis_bands = ["B2", "B02", "B3", "B03", "B4", "B04", "B5", "B05", "B8A", "B8A", "B11"]
            filtered_dict = {key: value for key, value in metadata_6sv["General_Info"]["bandname"].items() if
                             any(band in value for band in aquavis_bands)}

            updated_dict = {key: value.replace(".jp2", ".tif") for key, value in filtered_dict.items()}
            filtered_dict = updated_dict
            band_names = list(filtered_dict.values())

        else:

            aquavis_bands = ["B2", "B02", "B3", "B03", "B4", "B04", "B5", "B05", "B6", "B06"]
            filtered_dict = {key: value for key, value in metadata_6sv["General_Info"]["bandname"].items() if
                             any(band in value for band in aquavis_bands)}

            updated_dict = {key: value.replace(".TIF", ".tif") for key, value in filtered_dict.items()}
            filtered_dict = updated_dict

            band_names = list(filtered_dict.values())

        band_index = list(filtered_dict.keys())

        for i in range(len(band_index)):
            band_key = band_index[i]

            # Environmental function
            Fr = adjc_fun.atmospheric_point_scattering_function(metadata_6sv, band_key)

            # Surface Reflectance image without the adjacent correction
            with rasterio.open(os.path.join(self.image_path, filtered_dict[band_key])) as src:
                image_rs_no_adjc = src.read(1)

                image_rs_no_adjc = np.where(image_rs_no_adjc == -9999, np.nan, image_rs_no_adjc)
                image_rs_no_adjc = np.where(image_rs_no_adjc == 0, np.nan, image_rs_no_adjc)

                crs_image = src.crs
                transform = src.transform
                nodata_value = src.nodata

            # Apply the convolution to calculate the <p>
            convolved_image = self.fast_convolution(image_rs_no_adjc, Fr)

            one_value_array = np.full_like(image_rs_no_adjc, 1)
            convolved_Fr = self.fast_convolution(one_value_array, Fr)

            rho_env = convolved_image / convolved_Fr

            # Apply the adjacency correction
            sr_corr = adjc_fun.Adjacency_correction(metadata_6sv, band_key, image_rs_no_adjc, rho_env)

            output_tif = os.path.join(self.output_path, f"{band_names[i]}")

            # Open a new file to write the adjusted raster
            with rasterio.open(output_tif, 'w', driver='GTiff',
                               count=1, dtype=sr_corr.dtype,
                               width=sr_corr.shape[1], height=sr_corr.shape[0],
                               crs=crs_image, transform=transform, compress='lzw', nodata=nodata_value) as dst:
                dst.write(sr_corr, 1)  # Write the data to the file

    def run(self, params):

        if params['aux_info']['sat_name'] == "sentinel":

            sentinel_scene_dir = os.path.join(params['output_dir'], 'atmcor')
            input_path = [os.path.join(sentinel_scene_dir, scene) for scene in os.listdir(sentinel_scene_dir)]

        else:

            landsat_path = os.path.join(params['output_dir'], 'atmcor')
            input_path = [os.path.join(landsat_path, scene) for scene in os.listdir(landsat_path)]

        params['output_dir_adjcorr'] = os.path.join(params['output_dir'], 'adjcorr')
        satwutils.create_dir(params['output_dir_adjcorr'])

        for scene in input_path:

            if params['aux_info']['sat_name'] == "sentinel":
                scene_name = glob.glob(os.path.join(scene, '*.SAFE*'))[0]
            else:
                scene_name = glob.glob(os.path.join(scene, 'LC*'))[0]

            output_path = os.path.join(params['output_dir_adjcorr'], os.path.basename(scene_name))

            if os.path.exists(output_path):  # Check if the path exists
                print(f"Skipping {output_path}, already exists.")
                continue

            satwutils.create_dir(output_path)

            self.apply_adjcorr(scene_name, output_path, params)