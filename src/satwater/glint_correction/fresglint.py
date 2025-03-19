import os
import rasterio
import numpy as np
import rioxarray as rxr
import xml.etree.ElementTree as ET

from src.satwater.glint_correction import glintcorr_toolbox as toolbox

class FresGLINT:

    def __init__(self, path_main: str, dest: str):

        self.path_main = path_main
        self.dest = dest

    def read_metadata(self):
        tree = ET.parse(os.path.join(self.path_main, 'MTD.xml'))
        root = tree.getroot()
        return toolbox.xml_to_dict(root)

    def set_refraction_index(self, metadata):
        if metadata["General_Info"]["satellite"] == "MSI_S2":
            return {0: 1.335982329, 1: 1.332901801, 2: 1.329898916, 3: 1.326184102, 4: 1.313277663}
        else:
            return {0: 1.336521293, 1: 1.332967903, 2: 1.330137378, 3: 1.32620903, 4: 1.31328015}

    def calculate_glint_mask(self, red_band, swir_band):

        '''
        Define the region where the glint correction will be applied
        '''

        xda_red = rxr.open_rasterio(os.path.join(self.path_main, red_band))
        swir_20 = rxr.open_rasterio(os.path.join(self.path_main, swir_band))
        arr1020 = swir_20.rio.reproject_match(xda_red).values.astype(float)

        return np.where(arr1020 >= 0.005, 1, 0)

    def calculate_angle_images(self, metadata, arr1020):

        angles = {}

        for i in metadata["InputData"]["geometry"]["B0"].keys():

            angle_value = metadata["InputData"]["geometry"]["B0"][i]

            arr = np.full_like(arr1020, angle_value, dtype=float)

            angles[i] = np.where(arr1020 == -9999, np.nan, arr)

        return angles


    def reference_band(self, angle, metadata):

        if metadata["General_Info"]["satellite"] == "MSI_S2": swir_ref_band = ["B11"]
        else: swir_ref_band = ["B6", "B06"]

        filtered_dict_ref = list({key: value for key, value in metadata["General_Info"]["bandname"].items() if
                                  any(band in value for band in swir_ref_band)}.keys())[0]

        amtcor_input = metadata["InputData"]["sixSV_params"][filtered_dict_ref]

        reference = toolbox.paramglint(angle, metadata["InputData"]["geometry"][filtered_dict_ref]['solar_zn'],
                               metadata["InputData"]["sixSV_params"][filtered_dict_ref]['view_z'],
                               amtcor_input['optical_depth_total'], self.w_refIndex[4])

        return reference

    def corr_swir_glint(self, band_names, band_key, filtered_dict):

        array_path = os.path.join(self.path_main, filtered_dict[band_key])

        with rasterio.open(array_path) as src:

            arr = src.read(1).astype(float)
            profile = src.profile

        mask_nodata = np.where((arr < 0) | (arr > 1) | (arr == np.nan), 0, 1)
        # arr_integer = (arr * 10000).astype(np.int16)
        arr_integer = np.where(mask_nodata == 0, -9999, arr)

        # Save corrected image
        profile.update(
            dtype=arr.dtype.name,
            count=1,
            compress='lzw',
            driver='GTiff',
            nodata=-9999
        )

        output_path_save = os.path.join(self.dest, f"glintcorr_{band_names[i]}")

        with rasterio.open(output_path_save, 'w', **profile) as dest_dataset:
            dest_dataset.write(arr_integer, 1)

    def process_band(self, band_key, index, filtered_dict, metadata, swir_band, ang, w_refIndex):

        band_names = list(filtered_dict.values())
        band_index = list(filtered_dict.keys())

        if (metadata["General_Info"]["satellite"] == "MSI_S2") and (band_key == "B8"):

            arr1020 = rxr.open_rasterio(os.path.join(self.path_main, swir_band))
            gmask20 = np.where(arr1020 >= 0.03, 1, 0)

            ang20 = self.calculate_angle_images(self, metadata, arr1020)

            reference_20 = self.reference_band(ang20, metadata)

            target = toolbox.paramglint(ang20, metadata["InputData"]["geometry"][band_key]['solar_zn'],
                                metadata["InputData"]["sixSV_params"][band_key]['view_z'],
                                metadata["InputData"]["sixSV_params"][band_key]['optical_depth_total'],
                                w_refIndex[index])

            r_glint = arr1020 * (target['Tdir'] / reference_20['Tdir']) * (target['rFresnel'] / reference_20['rFresnel'])
            m_glint = gmask20 == 1
            r_glint = r_glint.values[0, :, :]

        elif (metadata["General_Info"]["satellite"] == "MSI_S2") and (band_key == "B11"):

            self.corr_swir_glint(band_names, band_key, filtered_dict)
            return None

        else:

            target = toolbox.paramglint(ang, metadata["InputData"]["geometry"][band_key]['solar_zn'],
                                metadata["InputData"]["sixSV_params"][band_key]['view_z'],
                                metadata["InputData"]["sixSV_params"][band_key]['optical_depth_total'],
                                w_refIndex[index])

            r_glint = self.arr1020 * (target['Tdir'] / self.reference['Tdir']) * (target['rFresnel'] / self.reference['rFresnel'])
            m_glint = self.gmask == 1
            r_glint = r_glint[0, :, :]

        array_path = os.path.join(self.path_main, filtered_dict[band_key])

        with rasterio.open(array_path) as src:
            arr = src.read(1).astype(float)
            profile = src.profile

        r_corr = self.correct_glint(arr, r_glint, m_glint)

        self.save_corrected_image(r_corr, profile, band_key)

    def correct_glint(self, arr, r_glint, m_glint):

        r_corr = np.copy(arr)
        m_glint = m_glint[0, :, :]
        r_corr[m_glint] = arr[m_glint] - r_glint[m_glint]
        r_corr[(r_corr > -0.2) & (r_corr < 0)] = 0.0001
        return r_corr

    def save_corrected_image(self, r_corr, profile, band_key):

        mask_nodata = np.where((r_corr < 0) | (r_corr > 1) | (r_corr == np.nan), 0, 1)
        # arr_integer = (r_corr_final * 10000).astype(np.int16)
        arr_integer = np.where(mask_nodata == 0, -9999, r_corr)

        profile.update(
            dtype=r_corr.dtype.name,
            count=1,
            compress='lzw',
            driver='GTiff',
            nodata=-9999
        )

        output_path_save = os.path.join(self.dest, f"glintcorr_{band_key}")
        with rasterio.open(output_path_save, 'w', **profile) as dest_dataset:
            dest_dataset.write(arr_integer, 1)

    def run(self):

        self.band_path = [i for i in os.listdir(self.path_main)]

        # Read MTD.xml file:
        metadata = self.read_metadata()

        # Set refraction index for MSI or OLI:
        self.w_refIndex = self.set_refraction_index(metadata)

        red_band = next((band for band in self.band_path if "B04" in band or "B4" in band), None)

        if metadata["General_Info"]["satellite"] == "MSI_S2":
            swir_band = next((band for band in self.band_path if "B11" in band), None)
            corr_bands = ["B02", "B2", "B03", "B3", "B04", "B8A", "B11"]
        else:
            swir_band = next((band for band in self.band_path if "B6" in band or "B06" in band), None)
            corr_bands = ["B02", "B2", "B03", "B3", "B4", "B04", "B5", "B05", "B6", "B06"]

        # Load the 1020nm band -> reference band:
        self.arr1020 = rxr.open_rasterio(os.path.join(self.path_main, swir_band)).rio.reproject_match(
            rxr.open_rasterio(os.path.join(self.path_main, red_band))
        ).values.astype(float)

        # Glint mask
        self.gmask = self.calculate_glint_mask(red_band, swir_band)

        # Angle images -> OAA, OZA, SAA, and SZA:
        self.ang = self.calculate_angle_images(metadata, self.arr1020)

        # Extracts the fresnel reflectance and transmittance for reference band (1020nm):
        self.reference = self.reference_band(metadata)

        # Glint correction for each band
        filtered_dict = {key: value for key, value in metadata["General_Info"]["bandname"].items() if
                         any(band in value for band in corr_bands)}

        updated_dict = {key: value.replace(".jp2", ".tif") for key, value in filtered_dict.items()}

        filtered_dict = updated_dict

        band_names = list(filtered_dict.values())
        band_index = list(filtered_dict.keys())

        for i in range(len(band_index)):

            band_key = band_index[i]

            print(f"Processing band: {band_key}")

            self.process_band(band_key, i, filtered_dict, metadata, swir_band, self.ang, self.w_refIndex, self.dest)