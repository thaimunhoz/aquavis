import os
import rasterio
import numpy as np
import xarray as xr
import multiprocessing
import rioxarray as rxr
from scipy.ndimage import median_filter
from src.satwater.utils import satwutils


def apply_glint_correction(input_dir, params):

    """
        Applies sunglint correction to satellite images based on the Wang and Shi, 2007 aproach.

        Parameters:
            scene_dir (str): Path to the directory containing the satellite image bands.
            params (dict): Dictionary containing processing parameters.
            mndwi_threshold (float, optional): Threshold value for MNDWI. Defaults to 0.3.
        """

    sat = params["sat"]

    bands = ['B2', 'B3', 'B4', 'B5', 'B6', 'B02', 'B03', 'B04', 'B8A', 'B11']

    input_bands = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if any(band in f for band in bands)]

    if sat == 'landsat':
        swir_band = next(os.path.join(input_dir, f) for f in os.listdir(input_dir) if 'B6' in f)

        base_dir = os.path.join(params['output_dir_glint'], 'landsat', os.path.basename(os.path.dirname(input_bands[0])))
        satwutils.create_dir(base_dir)

    else:
        swir_band = next(os.path.join(input_dir, f) for f in os.listdir(input_dir) if 'B11' in f)
        base_dir = os.path.join(params['output_dir_glint'], 'sentinel', os.path.basename(os.path.dirname(input_bands[0])))
        satwutils.create_dir(base_dir)

    xda_swir = rxr.open_rasterio(swir_band)

    for band_path in input_bands:

        if 'B6' in band_path or 'B11' in band_path:

            xda_ = rxr.open_rasterio(band_path)
            file_name_without_ext = os.path.basename(band_path)
            output_path = f"{base_dir}/{file_name_without_ext}"
            xda_.rio.to_raster(output_path)

        else:
            xda_ = rxr.open_rasterio(band_path)
            band_glint_corrected = np.where((xda_ - xda_swir) < 0, xda_, (xda_ - xda_swir))
            xda_glint = xr.DataArray(band_glint_corrected, dims=xda_.dims, coords=xda_.coords, attrs=xda_.attrs)

            file_name_without_ext = os.path.basename(band_path)
            output_path = f"{base_dir}/{file_name_without_ext}"

            xda_glint.rio.to_raster(output_path)

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

        # with multiprocessing.Pool(processes=params['aux_info']['n_cores']) as pool:
        #     results = pool.starmap_async(apply_glint_correction, zip(scenes, n_params)).get()
        #     print(results)
        #     pool.close()