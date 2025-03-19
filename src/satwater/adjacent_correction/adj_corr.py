import os
import rasterio
import numpy as np
from rasterio.mask import mask
from scipy.signal import convolve2d

from src.satwater.adjacent_correction import adjcorr_functions as adjc_fun

class AdjCorrClass:

    def __int__(self, input_image: str, output_path: str, tree:str):

        self.tree = tree
        self.output_path = output_path
        self.image_path = input_image

    def run(self):

        root = self.tree.getroot()
        metadata_6sv = adjc_fun.xml_to_dict(root)

        # Identify the pixels in which we're gonna apply the adjacent correction
        # The first gdf is the water mask, the second is the buffer. We will apply the correction in the buffer, and then clip the water mask to drop pixels in the border.
        water_shp_final = adjc_fun.get_mask_water(self.image_path)

        # Selecting only those important band for AQUAVis
        if metadata_6sv["General_Info"]["satellite"] == "MSI_S2":

            aquavis_bands = ["B2", "B02", "B3", "B03", "B4", "B04", "B5", "B05", "B8A", "B8A", "B11"]
            filtered_dict = {key: value for key, value in metadata_6sv["General_Info"]["bandname"].items() if
                             any(band in value for band in aquavis_bands)}

            updated_dict = {key: value.replace(".jp2", ".TIF") for key, value in filtered_dict.items()}
            filtered_dict = updated_dict
            band_names = list(filtered_dict.values())

        else:

            aquavis_bands = ["B2", "B02", "B3", "B03", "B4", "B04", "B5", "B05", "B6", "B06"]
            filtered_dict = {key: value for key, value in metadata_6sv["General_Info"]["bandname"].items() if
                             any(band in value for band in aquavis_bands)}

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

                if len(water_shp_final) != 0:  # in case of water_shp_final return empty
                    polygon_mask, _ = mask(src, water_shp_final.geometry, crop=False, nodata=src.nodata)
                    polygon_mask = polygon_mask[0, :, :]
                    polygon_mask_binary = np.where(polygon_mask == 0, 0, 1)

                else:
                    polygon_mask_binary = np.zeros_like(image_rs_no_adjc)

            # Apply the convolution to calculate the <p>
            convolved_image = convolve2d(image_rs_no_adjc, Fr, mode='same', boundary='fill', fillvalue=0)

            one_value_array = np.full_like(image_rs_no_adjc, 1)
            convolved_Fr = convolve2d(one_value_array, Fr, mode='same', boundary='fill', fillvalue=0)

            rho_env = convolved_image / convolved_Fr

            # Apply the adjacency correction
            sr_corr = adjc_fun.Adjacency_correction(metadata_6sv, band_key, image_rs_no_adjc, rho_env)

            # Only update pixels where mask == 1
            adjc_sr = np.where(polygon_mask_binary == 1, sr_corr, image_rs_no_adjc)

            output_tif = os.path.join(self.output_path, f"corrected_band_{band_names[i]}")

            # Open a new file to write the adjusted raster
            with rasterio.open(output_tif, 'w', driver='GTiff',
                               count=1, dtype=sr_corr.dtype,
                               width=adjc_sr.shape[1], height=adjc_sr.shape[0],
                               crs=crs_image, transform=transform, nodata=nodata_value) as dst:
                dst.write(adjc_sr, 1)  # Write the data to the file
