import os
import glob
import os.path
import multiprocessing

from src.satwater.utils import satwutils
import bandpass

def gen_resample(sentinel_scene, params):

    '''
    Resample the Sentinel-2 bands to the Landsat spatial resolution.
    '''

    sentinel_bands = ['B02', 'B03', 'B04', 'B8A']

    imgtemp_dir = fr"{params['output_dir']}\temp\temp_{os.path.basename(sentinel_scene)}"
    satwutils.create_dir(imgtemp_dir)

    sentinel_scene_all_bands = [f for f in glob.glob((fr'{sentinel_scene}\*.SAFE*\*_B*.tif')) if any(band in f for band in sentinel_bands)]

    for sentinel_band in sentinel_scene_all_bands:

        base_dir_nm = fr"{params['output_dir_tiling']}\sentinel\{os.path.basename(os.path.dirname(sentinel_band))}"
        satwutils.create_dir(base_dir_nm)
        out_scene = fr"{base_dir_nm}\{os.path.basename(sentinel_band)}"

        if os.path.exists(out_scene):
            continue

        resampling_imgtemp = fr"{imgtemp_dir}\{os.path.basename(sentinel_band)}"
        satwutils.cut_images_res(sentinel_band, params['sen_tile_target_shp'], resampling_imgtemp, 30)

        bandpass.apply_bandpass(resampling_imgtemp, out_scene)

def run(params):

    sentinel_bands = ['B02', 'B03', 'B04', 'B8A']

    _imgtemp_dir = fr"{params['output_dir']}\temp"
    os.makedirs(_imgtemp_dir, exist_ok=True)

    sen_tiles = params['sentinel']['tiles']

    for sen_tile_target in sen_tiles:

        # sentinel tile and projection are my reference
        sentinel_img = [f for f in glob.glob((fr'{params["output_dir"]}\atmcor\sentinel\{sen_tile_target}\**\*.SAFE*\*_B*.tif')) if any(band in f for band in sentinel_bands)][0]
        params['sen2_epsg_code'] = satwutils.raster2meta(sentinel_img)
        params['sen_tile_target_shp'] = satwutils.get_tile_shp(sen_tile_target, params, params['sen2_epsg_code'])

        # Create the output directory if it does not exist
        params['output_dir_tiling'] = fr'{params["output_dir"]}\tiling\{sen_tile_target}'
        satwutils.create_dir(fr'{params["output_dir_tiling"]}\sentinel')

        # Landsat spatial resolution is my reference
        sen2_scene_dir = fr'{params["output_dir"]}\atmcor\sentinel\{sen_tile_target}'
        sentinel_scenes = [fr'{sen2_scene_dir}\{i}' for i in os.listdir(sen2_scene_dir)]

        n_params = [params] * len(sentinel_scenes)

        with multiprocessing.Pool(processes=params['aux_info']['n_cores']) as pool:
            results = pool.starmap_async(gen_resample, zip(sentinel_scenes, n_params)).get()
            print(results)
            pool.close()