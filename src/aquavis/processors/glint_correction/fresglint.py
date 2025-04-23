import os
import logging
import rasterio
import numpy as np
import rioxarray as rxr
import xml.etree.ElementTree as ET

from src.aquavis.utils import io
from src.aquavis.processors.data_class import AquaVisDataLoader
from src.aquavis.processors.glint_correction import glintcorr_toolbox as toolbox

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class FresGLINT:

    def __init__(self):
        self.metadata = AquaVisDataLoader().load_aquavis_data()
        pass

    def set_refraction_index(self):
        if self.atmosphere_parameters["General_Info"]["satellite"] == "MSI_S2":
            return {0: 1.335982329, 1: 1.332901801, 2: 1.329898916, 3: 1.326184102, 4: 1.313277663}
        else:
            return {0: 1.336521293, 1: 1.332967903, 2: 1.330137378, 3: 1.32620903, 4: 1.31328015}

    def paramglint(self, ang: dict, solar_zn: float, view_zn: float, optical_depth_total: float, nw: float) -> dict:

        """Returns the Fresnel's reflectance and direct transmittance from atmosphere."""

        # Incidence angle:
        raa = abs(ang['solar_az'] - ang['view_az'])

        raa = np.where(raa > 180, 360 - raa, raa)

        cosTheta = np.sqrt((((np.cos(ang['solar_zn'] * (np.pi / 180)) * np.cos(ang['view_z'] * (np.pi / 180))) + (
                np.sin(ang['solar_zn'] * (np.pi / 180)) * np.sin(ang['view_z'] * (np.pi / 180)) * np.cos(
            raa * (np.pi / 180)))) + 1) / 2)  # cos2Theta = 2cos^2(Theta) - 1:

        Theta = np.arccos(cosTheta)

        # Transmittance angle:
        Theta_t = np.arcsin((1 / nw) * np.sin(Theta))

        # Fresnel reflectance coefficient:
        rfresnel = 0.5 * (((np.sin(Theta - Theta_t) / np.sin(Theta + Theta_t)) ** 2) + (
                (np.tan(Theta - Theta_t) / np.tan(Theta + Theta_t)) ** 2))

        # Direct Transmittance:
        Tdir = (np.exp(-float(optical_depth_total) / np.cos(float(solar_zn) * (np.pi / 180)))) * (
            np.exp(-float(optical_depth_total) / np.cos(float(view_zn) * (np.pi / 180))))

        return {'rFresnel': rfresnel, 'Tdir': Tdir}

    def reference_band(self, angle):

        w_refIndex = self.set_refraction_index()

        if self.atmosphere_parameters["General_Info"]["satellite"] == "MSI_S2": swir_ref_band = ["B11"]
        else: swir_ref_band = ["B6", "B06"]

        filtered_dict_ref = list({key: value for key, value in self.atmosphere_parameters["General_Info"]["bandname"].items() if
                                  any(band in value for band in swir_ref_band)}.keys())[0]

        return self.paramglint(angle,
                               self.atmosphere_parameters["InputData"]['sixSV_params'][filtered_dict_ref]['solar_zn'],
                                self.atmosphere_parameters["InputData"]['sixSV_params'][filtered_dict_ref]['view_z'],
                                self.atmosphere_parameters["InputData"]['sixSV_params'][filtered_dict_ref]['optical_depth_total'],
                                w_refIndex[4])

    def corr_swir_glint(self, input_path, output_path, band_name, band_key, filtered_dict):

        array_path = os.path.join(input_path, filtered_dict[band_key])

        with rasterio.open(array_path) as src:

            arr = src.read(1).astype(float)
            profile = src.profile

        mask_nodata = np.where((arr < 0) | (arr > 1) | (arr == np.nan), 0, 1)
        arr_integer = np.where(mask_nodata == 0, 0.0001, arr)
        arr_integer = arr_integer.squeeze()

        # Save corrected image
        profile.update(
            dtype=arr.dtype.name,
            count=1,
            compress='lzw',
            driver='GTiff',
            nodata=-9999
        )

        output_path_save = os.path.join(output_path, f"{band_name}")

        with rasterio.open(output_path_save, 'w', **profile) as dest_dataset:
            dest_dataset.write(arr_integer, 1)

    def process_band(self, input_path, output_path, band_key, index, band_name, filtered_dict, swir_band, ang):

        w_refIndex = self.set_refraction_index()

        if (self.atmosphere_parameters["General_Info"]["satellite"] == "MSI_S2") and (band_key == "B8"):

            arr1020 = rxr.open_rasterio(os.path.join(input_path, swir_band))
            arr1020 = np.nan_to_num(arr1020, nan=0.0001)  # Replace NaNs with 0.0001

            arr1020[arr1020 <= 0] = 0.0001

            ang20 = toolbox.calculate_angle_images(self.atmosphere_parameters, arr1020)

            reference_20 = self.reference_band(ang20)

            target = self.paramglint(ang20, self.atmosphere_parameters["InputData"]['sixSV_params'][band_key]['solar_zn'],
                                self.atmosphere_parameters["InputData"]['sixSV_params'][band_key]['view_z'],
                                self.atmosphere_parameters["InputData"]['sixSV_params'][band_key]['optical_depth_total'],
                                w_refIndex[index])

            r_glint = arr1020 * (target['Tdir'] / reference_20['Tdir']) * (target['rFresnel'] / reference_20['rFresnel'])
            r_glint = r_glint

        elif (self.atmosphere_parameters["General_Info"]["satellite"] == "MSI_S2") and (band_key == "B11"):

            self.corr_swir_glint(input_path, output_path, band_name, band_key, filtered_dict)
            return None

        elif (band_key == "B5") or (band_key == "B05"):

            self.corr_swir_glint(input_path, output_path, band_name, band_key, filtered_dict)
            return None

        else:

            target = self.paramglint(ang, self.atmosphere_parameters["InputData"]['sixSV_params'][band_key]['solar_zn'],
                                self.atmosphere_parameters["InputData"]['sixSV_params'][band_key]['view_z'],
                                self.atmosphere_parameters["InputData"]['sixSV_params'][band_key]['optical_depth_total'],
                                w_refIndex[index])

            r_glint = self.arr_ref * (target['Tdir'] / self.reference['Tdir']) * (target['rFresnel'] / self.reference['rFresnel'])
            r_glint = r_glint

        array_path = os.path.join(input_path, filtered_dict[band_key])

        with rasterio.open(array_path) as src:
            arr = src.read(1).astype(float)
            profile = src.profile

        arr_aux = rxr.open_rasterio(array_path)
        self.vnir_values.append(((arr_aux >= -0.1) & (arr_aux < 0)).astype(int))

        r_corr = self.correct_glint(arr, r_glint)

        toolbox.save_corrected_image(r_corr, profile, band_name, output_path)

    def correct_glint(self, arr, r_glint):

        r_corr = arr - r_glint
        r_corr[(r_corr > -0.2) & (r_corr < 0)] = 0.0001
        return r_corr

    def run(self, input_path: str, output_path: str) -> None:

        '''Apply glint correction based on Fresnel Reflectance'''

        tree = ET.parse(os.path.join(input_path, "MTD.xml"))
        root = tree.getroot()
        self.atmosphere_parameters = toolbox.xml_to_dict(root)

        io.validate_file(self.metadata, input_path)

        band_path = [i for i in os.listdir(input_path)]

        # This band is used as reference to resample SWIR band to 10m
        red_band = next((band for band in band_path if "B04" in band or "B4" in band), None)

        # Depending on the sensor, bands name will change
        corr_bands, swir_band = toolbox.return_bands(self.atmosphere_parameters, band_path)

        # Load the SWIR band -> reference band:
        self.arr_ref = toolbox.load_reference_band(input_path, swir_band, red_band)

        #self.glint_mask = (self.arr_ref > 0.0157).astype(int)

        # Angle images -> OAA, OZA, SAA, and SZA:
        ang = toolbox.calculate_angle_images(self.atmosphere_parameters, self.arr_ref)

        # Extracts the fresnel reflectance and transmittance for reference band
        self.reference = self.reference_band(ang)

        # Glint correction for each band
        filtered_dict = {key: value for key, value in self.atmosphere_parameters["General_Info"]["bandname"].items() if
                         any(band in value for band in corr_bands)}

        updated_dict = {key: value.replace(".jp2", ".tif") for key, value in filtered_dict.items()}
        updated_dict = {key: value.replace(".TIF", ".tif") for key, value in updated_dict.items()}

        filtered_dict = updated_dict

        band_names = list(filtered_dict.values())
        band_index = list(filtered_dict.keys())

        self.vnir_values = []

        os.makedirs(output_path, exist_ok=True)

        for i in range(len(band_index)):

            band_key = band_index[i]
            band_name = band_names[i]

            self.process_band(input_path, output_path, band_key, i, band_name, filtered_dict, swir_band, ang)

        logger.info("Glint correction completed successfully")