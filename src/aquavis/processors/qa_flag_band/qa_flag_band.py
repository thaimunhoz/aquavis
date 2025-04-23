import os
import numpy as np
import rioxarray as rxr
from scipy.ndimage import convolve

def surrounded_cloud(cloud_mask):

    kernel = np.ones((3, 3), dtype=int)

    cloud_neighbors = convolve(cloud_mask.values.astype(int), kernel, mode='constant', cval=0)

    surrounded_by_clouds = (cloud_neighbors > 4).astype(int)

    surrounded_by_clouds_raster = cloud_mask.copy()
    surrounded_by_clouds_raster.values = surrounded_by_clouds

    return surrounded_by_clouds_raster

def process_masks(main_path: str) -> dict:

    files = os.listdir(main_path)

    # Find water band
    water_band = next((band for band in files if "_water" in band), None)
    water_mask = rxr.open_rasterio(os.path.join(main_path, water_band)) if water_band else None

    # Initialize masks dictionary
    masks = {
        "water": None,
        "cloud": None,
        "shadow": None,
        "glint": None,
        "swir_negative": None,
        "blue_negative": None,
        "green_negative": None,
        "red_negative": None,
        "nir_negative": None,
    }

    # File mapping for easy handling
    file_mapping = {
        "_cloud_shadow.tif": "cloud_shadow",
        "_glint_mask.tif": "glint",
        "_swir_negative.tif": "swir_negative",
        "_vnir_neg_0.tif": "blue_negative",
        "_vnir_neg_1.tif": "green_negative",
        "_vnir_neg_2.tif": "red_negative",
        "_vnir_neg_3.tif": "nir_negative",
    }

    masks["water"] = water_mask

    # Process files
    for file in files:
        for suffix, key in file_mapping.items():
            if file.endswith(suffix):
                file_path = os.path.join(main_path, file)
                data = rxr.open_rasterio(file_path)
                data = data.rio.reproject_match(water_mask)

                if key == "cloud_shadow":
                    masks["cloud"] = (data == 2).astype(int)
                    masks["shadow"] = (data == 3).astype(int)
                else:
                    masks[key] = data

    return masks

def create_qa_bitmask(water, cloud, shadow, glint_neg, neg_blue, neg_green, neg_red, neg_nir, neg_swir):

    # Encode cloud and shadow together (00 = clear, 01 = shadow, 10 = cloud, 11 = cloud + shadow)
    cloud_shadow = (cloud.astype(int) << 1) | (shadow.astype(int))

    # Construct the QA band with the new flags
    qa_band = (
            (water.astype(int) << 0) |  # Bit 0: Water
            (cloud_shadow << 1) |  # Bits 1-2: Cloud/Shadow Combined
            (glint_neg.astype(int) << 3) |  # Bit 3: Glint Correction
            (neg_blue.astype(int) << 4) |  # Bit 4: Negative Blue
            (neg_green.astype(int) << 5) |  # Bit 5: Negative Green
            (neg_red.astype(int) << 6) |  # Bit 6: Negative Red
            (neg_nir.astype(int) << 7) |  # Bit 7: Negative NIR
            (neg_swir.astype(int) << 8)  # Bit 8: Negative SWIR
    )

    return qa_band

def run_qa_flag(params, scene_path):

    input_path = os.path.join(params.output_dir_qaflag, os.path.basename(scene_path))

    masks = process_masks(input_path)
    qa_band = create_qa_bitmask(masks["water"], masks["cloud"], masks["shadow"], masks["glint"], masks["blue_negative"], masks["green_negative"], masks["red_negative"], masks["nir_negative"], masks["swir_negative"])

    return qa_band