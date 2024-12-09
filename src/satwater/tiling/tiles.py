import os
import ast
import glob
import os.path
import multiprocessing

import pandas as pd

from src.satwater.utils import satwutils

def gen_tiles(landsat_scene, params):

    '''
    Generate the tiles for the Landsat scene. The Landsat bands are reprojected, resampled, and clipped to the Sentinel tile (MGRS).
    For each Landsat tile, we look for all intersections with Sentinel tiles and generate the corresponding tiles.
    '''

    print(landsat_scene)

    landsat_bands = ['B2', 'B3', 'B4', 'B5']

    sen_tile = params['sen_tiles_target']

    #path = fr"{landsat_scene}\{os.listdir(landsat_scene)[1]}"

    landsat_scene_all_bands = [f for f in glob.glob(os.path.join(landsat_scene, '*')) if any(band in f for band in landsat_bands)]

    #sen_gdf = sen_tile[sen_tile['Name'] == sen_tile_target]

    sen2_epsg_code = params['sen2_epsg_code']

    for landsat_band in landsat_scene_all_bands:

        base_dir_nm = fr"{params['output_dir_tiling']}\landsat\{os.path.basename(os.path.dirname(landsat_band))}"
        satwutils.create_dir(base_dir_nm)
        out_band = fr"{base_dir_nm}\{os.path.basename(landsat_band)}"

        if os.path.exists(out_band):
            continue

        imgtemp = fr"{params['output_dir']}\temp\temp_{os.path.basename(landsat_band)}"
        satwutils.reproject(landsat_band, imgtemp, sen2_epsg_code)  # reprojecting the landsat to sentienl projection

        satwutils.cut_images_res(imgtemp, params['sen_tile_target_shp'], out_band, 30)

def run(select_sat, params):

    os.makedirs(fr"{params['output_dir']}\temp", exist_ok=True)

    sentinel_tiles_names = pd.read_csv(r'C:\Users\tml411\Documents\Python Scripts\hls_water\src\satwater\auxfiles\tiles\sentinel_landsat_intersections.csv')

    sen_tiles = [params[select_sat]['tiles']]

    for sen_tile_target in sen_tiles:

        sent_tile = sentinel_tiles_names[sentinel_tiles_names['sentinel_tile'] == sen_tile_target]

        params['sen_tiles_target'] = sent_tile['sentinel_tile'].values[0]

        # sentinel tile and projection are my reference
        #src_dir_tile = fr'{params["sentinel"]["input_dir"]}\{sen_tile_target}'
        #sentinel_path = \
        #[fr'{src_dir_tile}\{i}\{os.listdir(os.path.join(src_dir_tile, i))[0]}' for i in os.listdir(src_dir_tile)][0]
        #sentinel_img = glob.glob(f'{sentinel_path}\**\*B1*.jp2', recursive=True)[0]

        params['sen2_epsg_code'] = sent_tile["sentinel_epsg"].values[0]
        params['sen_tile_target_shp'] = satwutils.get_tile_shp(params['sen_tiles_target'], params, params['sen2_epsg_code'])

        # Create the output directory if it does not exist
        params['output_dir_tiling'] = fr'{params["output_dir"]}\tiling\{sen_tile_target}'
        satwutils.create_dir(fr'{params["output_dir_tiling"]}\landsat')

        #pathrows = ast.literal_eval(sent_tile["landsat_tiles"].values[0])
        #landsat_pathrows = [f"{int(tile.split('_')[0]):03d}_{int(tile.split('_')[1]):03d}" for tile in pathrows]

        path_pr = fr'{params["output_dir"]}\atmcor\landsat\{sent_tile["sentinel_tile"].values[0]}'
        landsat_scenes_aux = [fr'{path_pr}\{i}' for i in os.listdir(path_pr)]
        landsat_scenes = [fr"{landsat_scene}\{os.listdir(landsat_scene)[1]}" for landsat_scene in landsat_scenes_aux]

        n_params = [params] * len(landsat_scenes)

        with multiprocessing.Pool(processes=params['aux_info']['n_cores']) as pool:
            results = pool.starmap_async(gen_tiles, zip(landsat_scenes, n_params)).get()
            print(results)
            pool.close()