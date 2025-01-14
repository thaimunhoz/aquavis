import os
import rasterio
import numpy as np
import multiprocessing
from scipy.ndimage import median_filter
from src.satwater.utils import satwutils

def calculate_mndwi(green_band, swir_band):

    """
       Calculates the Modified Normalized Difference Water Index (MNDWI).

       Parameters:
           green_band (numpy.ndarray): Array of reflectance values from the green band.
           swir_band (numpy.ndarray): Array of reflectance values from the SWIR band.

       Returns:
           numpy.ndarray: MNDWI values.
       """

    mndwi = (green_band - swir_band) / (green_band + swir_band)
    return mndwi

def apply_glint_correction(input_dir, params, mndwi_threshold=0.3):

    """
        Applies sunglint correction to satellite images based on the MNDWI threshold.

        Parameters:
            scene_dir (str): Path to the directory containing the satellite image bands.
            params (dict): Dictionary containing processing parameters.
            mndwi_threshold (float, optional): Threshold value for MNDWI. Defaults to 0.3.
        """

    sat = params["sat"]

    bands = ['B2', 'B3', 'B4', 'B5', 'B02', 'B03', 'B04', 'B05', 'B8', 'B8A', 'B11', 'B12']

    input_bands = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if any(band in f for band in bands)]

    if sat == 'landsat':
        swir_band = next(os.path.join(input_dir, f) for f in os.listdir(input_dir) if 'B7' in f)
        green_band = next(os.path.join(input_dir, f) for f in os.listdir(input_dir) if 'B3' in f)

        base_dir = os.path.join(params['output_dir_glint'], 'landsat', os.path.basename(os.path.dirname(input_bands[0])))
        satwutils.create_dir(base_dir)

    else:
        green_band = next(os.path.join(input_dir, f) for f in os.listdir(input_dir) if 'B03' in f)
        swir_band = next(os.path.join(input_dir, f) for f in os.listdir(input_dir) if 'B12' in f)
        base_dir = os.path.join(params['output_dir_glint'], 'sentinel', os.path.basename(os.path.dirname(input_bands[0])))
        satwutils.create_dir(base_dir)

    with rasterio.open(green_band) as green_src, rasterio.open(swir_band) as swir_src:

        green = green_src.read(1).astype(np.float32)
        swir = swir_src.read(1).astype(np.float32)
        profile = green_src.profile

        # Apply median filter to SWIR band (3x3 window)
        filtered_swir = median_filter(swir, size=3)

        # Calculate MNDWI
        mndwi = calculate_mndwi(green, swir)

        # Define ROI where MNDWI > threshold
        roi_mask = mndwi > mndwi_threshold

        # Calculate beta (minimum SWIR value in ROI)
        beta = np.percentile(swir[roi_mask], 25)  # np.min(swir[roi_mask])

        for band_path in input_bands:
            with rasterio.open(band_path) as band_src:
                band = band_src.read(1).astype(np.float32)

                # Apply glint correction: reflec_corrected = reflec_original - alpha * (reflec_original(swir) - beta)
                #corrected_band = band - 1 * (filtered_swir - beta)
                corrected_band = (filtered_swir - beta)

                # Save corrected band
                file_name_without_ext = os.path.basename(band_path)
                output_path = f"{base_dir}/{file_name_without_ext}"
                profile.update(dtype=rasterio.float32)

                with rasterio.open(output_path, 'w', **profile) as dst:
                    dst.write(corrected_band, 1)

                print(f"Corrected band saved at: {output_path}")

def run(params):

    """
    Executes the sunglint correction process for a set of satellite scenes.

    Parameters:
        params (dict): Dictionary containing processing parameters.
    """

    sat = params['aux_info']['sat_name']

    if sat == 'landsat':
        sen_tiles = [params[sat]['tiles']]
    else:
        sen_tiles = [params['sentinel']['tiles']]

    for sen_tile_target in sen_tiles:

        params['sen_tile_target'] = sen_tile_target

        # Create the output directory if it does not exist
        params['output_dir_glint'] = fr'{params["output_dir"]}\glint\{sen_tile_target}'
        os.makedirs(fr'{params["output_dir_glint"]}', exist_ok=True)

        if sat == 'landsat':
            ncode = 'L30'
            path_pr = fr'{params["output_dir"]}\tiling\{params["sen_tile_target"]}\{sat}'
            scenes = [fr'{path_pr}\{i}' for i in os.listdir(path_pr)]
        else:
            ncode = 'S30'
            path_pr = fr'{params["output_dir"]}\tiling\{params["sen_tile_target"]}\{sat}'
            scenes = [fr'{path_pr}\{i}' for i in os.listdir(path_pr)]

        params['ncode'] = ncode
        params["sat"] = sat

        n_params = [params] * len(scenes)

        for scene in scenes:
            apply_glint_correction(scene, params)

        #with multiprocessing.Pool(processes=params['aux_info']['n_cores']) as pool:
        #    results = pool.starmap_async(apply_glint_correction, zip(scenes, n_params)).get()
        #    print(results)
        #    pool.close()