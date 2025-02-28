import os
import glob
import pandas as pd
import multiprocessing
from src.satwater.utils import satwutils
from src.satwater.tiling import brdf_lee11_QAA_RGB as brdf

def gen_tiles(landsat_scene: str, params: dict) -> None:

    """
    Generate tiles for a given Landsat scene by reprojecting, resampling, and clipping it to Sentinel tile geometry.

    Args:
        landsat_scene (str): Path to the Landsat scene directory.
        params (dict): Dictionary containing necessary parameters for tiling, including Sentinel tile geometry and EPSG code.

    Returns:
        None
    """

    print(f"Processing Landsat scene: {landsat_scene}")

    landsat_bands = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6']
    landsat_scene_all_bands = [
        f for f in glob.glob(os.path.join(landsat_scene, '*LC*', '*_B*.TIF')) if
        any(band in f for band in landsat_bands)
    ]

    landsat_scene_all_bands = sorted(
        landsat_scene_all_bands,
        key=lambda x: next((i for i, band in enumerate(landsat_bands) if band in x), float('inf'))
    )

    sen2_epsg_code = params['sen2_epsg_code']

    imgtemp_dir = os.path.join(params['output_dir'], 'temp', f"temp_{os.path.basename(landsat_scene)}")
    satwutils.create_dir(imgtemp_dir)

    for landsat_band in landsat_scene_all_bands:

        base_dir = os.path.join(params['output_dir_tiling'], 'landsat', os.path.basename(os.path.dirname(landsat_band)))
        satwutils.create_dir(base_dir)

        out_band = os.path.join(base_dir, os.path.basename(landsat_band))

        # Temporary file for reprojected image
        temp_img = os.path.join(imgtemp_dir, f"temp_{os.path.basename(landsat_band)}")

        # Reproject Landsat band to Sentinel projection
        satwutils.reproject(landsat_band, temp_img, sen2_epsg_code)

    # Apply BRDF correction
    if params['aux_info']['brdf_corr']:
        brdf.call_brdf_correction(imgtemp_dir, imgtemp_dir, 'landsat')

        path_out = imgtemp_dir
        images_brdf = [f for f in glob.glob(fr'{path_out}\*.TIF') if "brdf_corrected" in f]

    else:
        path_out = imgtemp_dir
        images_brdf = [f for f in glob.glob(fr'{path_out}\*.TIF') if "temp" in f]

    # Clip and resample image to match Sentinel tile
    i = 0

    for img in images_brdf:

        base_dir = os.path.join(params['output_dir_tiling'], 'landsat', os.path.basename(os.path.dirname(landsat_band)))
        satwutils.create_dir(base_dir)

        out_band = os.path.join(base_dir, os.path.basename(landsat_scene_all_bands[i]))

        satwutils.cut_images_res(img, params['sen_tile_target_shp'], out_band, 30)
        i += 1

def run(select_sat: str, params: dict) -> None:

    """
    Main function to generate tiles for all Landsat scenes corresponding to specified Sentinel tiles.

    Args:
        select_sat (str): Satellite type (e.g., 'landsat').
        params (dict): Dictionary containing input/output directories, tile information, and processing settings.

    Returns:
        None
    """

    os.makedirs(os.path.join(params['output_dir'], 'temp'), exist_ok=True)

    sentinel_tiles_names = pd.read_csv(r'C:\Users\tml411\Documents\Python Scripts\hls_water\src\satwater\auxfiles\tiles\sentinel_landsat_intersections.csv')

    sentinel_tiles = [params[select_sat]['tiles']]

    for sentinel_tile in sentinel_tiles:
        try:
            sent_tile_info = sentinel_tiles_names[sentinel_tiles_names['sentinel_tile'] == sentinel_tile]

            params['sen_tiles_target'] = sent_tile_info['sentinel_tile'].values[0]
            params['sen2_epsg_code'] = sent_tile_info['sentinel_epsg'].values[0]

            # Get Sentinel tile geometry
            params['sen_tile_target_shp'] = satwutils.get_tile_shp(
                params['sen_tiles_target'], params, params['sen2_epsg_code']
            )

            # Create output directory for tiling
            params['output_dir_tiling'] = os.path.join(params['output_dir'], 'tiling', sentinel_tile)
            satwutils.create_dir(os.path.join(params['output_dir_tiling'], 'landsat'))

            # Gather Landsat scenes for the corresponding Sentinel tile
            landsat_path = os.path.join(params['output_dir'], 'atmcor', 'landsat', sentinel_tile)
            landsat_scenes = [os.path.join(landsat_path, scene) for scene in os.listdir(landsat_path)]

            # Process Landsat scenes in parallel
            with multiprocessing.Pool(processes=params['aux_info']['n_cores']) as pool:
                pool.starmap(gen_tiles, [(scene, params) for scene in landsat_scenes])
                pool.close()

            print(f"Completed processing for Sentinel tile: {sentinel_tile}")

        except Exception as e:
            print(f"Error processing Sentinel tile {sentinel_tile}: {e}")