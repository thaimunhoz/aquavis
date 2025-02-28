import os
import glob
import shutil
import os.path
import rasterio
import numpy as np
import multiprocessing
from src.satwater.utils import satwutils

def convert_float_to_int16(input_raster, output_raster, params, scale_factor=10000, nodata=-9999):

    '''
    Convert a float raster to int16 by multiplying by a scale factor and setting nodata values.
    Input:
        input_raster: str, path to the input raster
        output_raster: str, path to the output raster
        scale_factor: int, scale factor to multiply the input raster by
        nodata: int, nodata value to set in the output raster
    Output:
        None
    '''

    with rasterio.open(input_raster, 'r') as src:

        # Read the input raster as a float array
        arr = src.read(1).astype('float32')
        arr[arr <= 0] = np.nan

        # Output in rho or rrs:
        if params['aux_info']['output_type'] == 'rho':
            arr = arr # surface reflectance
            arr *= scale_factor
        else:
            arr = arr/np.pi # Rrs

        # Multiply the array by 10000
        #arr *= scale_factor
        arr = np.where(np.isnan(arr), nodata, arr)

        # Copy the metadata from the input raster
        kwargs = src.meta.copy()

        # Update the data type, nodata value, and compression for the output raster
        kwargs.update({
            'dtype': rasterio.int16,
            'nodata': nodata,
            'compress': 'lzw'
        })

        # Create the output raster in write mode
        with rasterio.open(output_raster, 'w', **kwargs) as dst:

            # Write the modified array to the output raster
            dst.write(arr.astype(rasterio.int16), 1)

def get_scene_details(scene_path, sat='landsat'):

    """
    Returns a dictionary of scene details for a given scene path
    Input:
        scene_path: str, path to the scene
        sat: str, satellite name (landsat or sentinel)
    """

    if sat == 'sentinel':

        scene_name = os.path.basename(scene_path)
        date_time = scene_name.split('_')[2]

    else:

        scene_name = os.path.basename(scene_path)

        img = glob.glob(fr'{scene_path}\*.tif')[0]
        hh = os.path.basename(img).split("_")[5]
        mm = os.path.basename(img).split("_")[6]
        ss = os.path.basename(img).split("_")[7]

        date_time = scene_name.split('_')[3]
        #date_time = f"{date_time}T{hh}{mm}{ss}"

    return date_time

def gen_hls(scene_path, params):

    '''
    Generate HLS products from a given scene
    Input:
        scene_path: str, path to the scene
        params: dict, parameters for the HLS generation
    '''

    sentinel_bands = ['B02', 'B03', 'B04', 'B8A', 'B11', 'B12']

    date_time = get_scene_details(scene_path, params["sat"])

    if params["sat"] == 'landsat':
        landsat_tile = os.path.basename(scene_path).split('_')[2]
        out_dir_hls = fr"{params['output_dir_hls']}\HLS_T{params['sen_tile_target']}_{date_time}_{landsat_tile}_{params['ncode']}_v1.0"
    else:
        out_dir_hls = fr"{params['output_dir_hls']}\HLS_T{params['sen_tile_target']}_{date_time}_{params['ncode']}_v1.0"

    os.makedirs(out_dir_hls, exist_ok=True)

    scene_all_bands = glob.glob(f"{scene_path}\**.tif")
    scene_all_bands = sorted(scene_all_bands, key=lambda x: next((i for i, band in enumerate(sentinel_bands) if band in x), float('inf')))

    for file_bandi in scene_all_bands:

        band_nmi = os.path.basename(file_bandi).split("_")[-1].split(".")[0]

        # Fix band names:
        if band_nmi == 'B8A': band_nmi = 'B05'

        band_nmi = satwutils.fix_band_name(band_nmi, params['aux_info']['sat_name'])

        if params["sat"] == 'landsat':
            basefilename = f'HLS_T{params["sen_tile_target"]}_{date_time}_{landsat_tile}_{params["ncode"]}_v1.0_{band_nmi}.tif'
            cloudname = f'HLS_T{params["sen_tile_target"]}_{date_time}_{landsat_tile}_{params["ncode"]}_v1.0_CLOUD.tif'
        else:
            basefilename = f'HLS_T{params["sen_tile_target"]}_{date_time}_{params["ncode"]}_v1.0_{band_nmi}.tif'  ## HLS.T17SLU.2020209T155956.L30.v1.5.B01.tif
            cloudname = f'HLS_T{params["sen_tile_target"]}_{date_time}_{params["ncode"]}_v1.0_CLOUD.tif'

        out_base_dir_nm = fr"{out_dir_hls}\{basefilename}"

        # convert from float to int16, change range to 0-10000, nodata = -9999
        convert_float_to_int16(file_bandi, out_base_dir_nm, params, scale_factor=10000, nodata=-9999)

    # Save cloud information as a band
    cloud_path = os.path.join(params["output_dir"], "atmcor", "landsat", os.path.basename(params["output_dir"]), os.path.basename(scene_path), "SatClouds", "temp", "cloud.tif")
    destination_path = os.path.join(out_dir_hls, cloudname)

    try:
        shutil.copy(cloud_path, destination_path)
    except:
        print(f"No cloud information found for scene {scene_path}")

def run(params):

    """
    Runs a given set of parameters to initiate a selection process for satellite data.
    The 'select_sat' parameter defines the satellite to select (defaults to 'landsat')
    """

    sat = params['aux_info']['sat_name']

    if sat == 'landsat':
        sen_tiles = [params[sat]['tiles']]
    else:
        sen_tiles = [params['sentinel']['tiles']]

    for sen_tile_target in sen_tiles:

        params['sen_tile_target'] = sen_tile_target

        # Create the output directory if it does not exist
        params['output_dir_hls'] = fr'{params["output_dir"]}\hlswater'
        os.makedirs(fr'{params["output_dir_hls"]}', exist_ok=True)

        if sat == 'landsat':
            ncode = 'L30'
            path_pr = fr'{params["output_dir"]}\tiling\{params["sen_tile_target"]}\landsat'
            scenes = [fr'{path_pr}\{i}' for i in os.listdir(path_pr)]
        else:
            ncode = 'S30'
            path_pr = fr'{params["output_dir"]}\tiling\{params["sen_tile_target"]}\sentinel'
            scenes = [fr'{path_pr}\{i}' for i in os.listdir(path_pr)]

        params['ncode'] = ncode
        params["sat"] = sat

        n_params = [params] * len(scenes)

        # for scene in scenes:
        #     gen_hls(scene, params)

        with multiprocessing.Pool(processes=params['aux_info']['n_cores']) as pool:
            results = pool.starmap_async(gen_hls, zip(scenes, n_params)).get()
            print(results)
            pool.close()

        atmcor_folder = fr'{params["output_dir"]}\atmcor'
        tiling_folder = fr'{params["output_dir"]}\tiling'
        temp_folder = fr'{params["output_dir"]}\temp'

        # Clean folders
        try:
            if params['aux_info']["keep_atmor"] == True:

                os.chmod(tiling_folder, 0o777)
                os.chmod(temp_folder, 0o777)

                shutil.rmtree(tiling_folder)
                shutil.rmtree(temp_folder)

            else:

                os.chmod(atmcor_folder, 0o777)
                os.chmod(tiling_folder, 0o777)
                os.chmod(temp_folder, 0o777)

                shutil.rmtree(atmcor_folder)
                shutil.rmtree(tiling_folder)
                shutil.rmtree(temp_folder)
        except:
            pass