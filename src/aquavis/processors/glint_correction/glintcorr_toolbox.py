import os
import rasterio
import numpy as np
import rioxarray as rxr

def xml_to_dict(element):

    if len(element) == 0:
        return element.text

    result = {}

    for child in element:

        child_dict = xml_to_dict(child)

        if child.tag in result:
            if isinstance(result[child.tag], list):
                result[child.tag].append(child_dict)
            else:
                result[child.tag] = [result[child.tag], child_dict]
        else:
            result[child.tag] = child_dict

    return result

def load_reference_band(input_path, swir_band, red_band):
    """Load the reference band for glint correction."""
    glint_ref = rxr.open_rasterio(os.path.join(input_path, swir_band)).rio.reproject_match(
        rxr.open_rasterio(os.path.join(input_path, red_band))).values.astype(float)

    glint_ref = np.nan_to_num(glint_ref, nan=0.0001)  # Replace NaNs with 0.0001
    glint_ref[glint_ref <= 0] = 0.0001

    return glint_ref

def return_bands(metadata, band_path):

    if metadata["General_Info"]["satellite"] == "MSI_S2":
        swir_band = next((band for band in band_path if "B11" in band), None)
        corr_bands = ["B02", "B2", "B03", "B3", "B04", "B8A", "B11"]
    else:
        swir_band = next((band for band in band_path if "B6" in band or "B06" in band), None)
        corr_bands = ["B02", "B2", "B03", "B3", "B4", "B04", "B5", "B05", "B6", "B06"]

    return corr_bands, swir_band

def calculate_angle_images(metadata, arr1020):
    """Calculate angles for glint correction."""
    angles = {}

    for i in metadata["InputData"]['sixSV_params']["B0"].keys():
        angle_value = metadata["InputData"]['sixSV_params']["B0"][i]

        arr = np.full_like(arr1020, angle_value, dtype=float)

        angles[i] = np.where(arr1020 == -9999, np.nan, arr)

    return angles

def save_corrected_image(r_corr, profile, band_name, dest):

    mask_nodata = np.where((r_corr < 0) | (r_corr > 1) | (r_corr == np.nan), 0, 1)
    # arr_integer = (r_corr_final * 10000).astype(np.int16)
    arr_integer = np.where(mask_nodata == 0, 0.0001, r_corr)
    arr_integer = arr_integer.squeeze()

    profile.update(
        dtype=r_corr.dtype.name,
        count=1,
        compress='lzw',
        driver='GTiff',
        nodata=-9999
    )

    output_path_save = os.path.join(dest, f"{band_name}")
    with rasterio.open(output_path_save, 'w', **profile) as dest_dataset:
        dest_dataset.write(arr_integer, 1)

