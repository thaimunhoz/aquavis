import os
import glob
import multiprocessing
from src.satwater.utils import satwutils
from src.satwater.tiling import bandpass

def gen_resample(sentinel_scene, params):

    """
    Resamples Sentinel-2 bands to match Landsat spatial resolution and applies bandpass correction.

    Args:
        sentinel_scene (str): Path to the Sentinel-2 scene directory.
        params (dict): Dictionary of parameters containing output directories, tile shapefiles, etc.
    """

    sentinel_bands = ['B02', 'B03', 'B04', 'B8A', 'B11', 'B12']
    imgtemp_dir = os.path.join(params['output_dir'], 'temp', f"temp_{os.path.basename(sentinel_scene)}")
    satwutils.create_dir(imgtemp_dir)

    sentinel_scene_bands = [
        f for f in glob.glob(os.path.join(sentinel_scene, '*.SAFE*', '*_B*.tif'))
        if any(band in f for band in sentinel_bands)
    ]
    sentinel_scene_bands = sorted(
        sentinel_scene_bands, key=lambda x: next((i for i, band in enumerate(sentinel_bands) if band in x), float('inf'))
    )

    for sentinel_band in sentinel_scene_bands:
        output_dir = os.path.join(params['output_dir_tiling'], 'sentinel', os.path.basename(os.path.dirname(sentinel_band)))
        satwutils.create_dir(output_dir)

        output_path = os.path.join(output_dir, os.path.basename(sentinel_band))

        if os.path.exists(output_path):
            continue

        temp_path = os.path.join(imgtemp_dir, os.path.basename(sentinel_band))
        satwutils.cut_images_res(sentinel_band, params['sen_tile_target_shp'], temp_path, 30)

        bandpass.apply_bandpass(temp_path, output_path)

def run(params):

    """
    Coordinates the resampling of Sentinel-2 bands for multiple tiles.

    Args:
        params (dict): Dictionary of parameters containing output directories, tile shapefiles, and settings.
    """

    sentinel_bands = ['B02', 'B03', 'B04', 'B8A', 'B11', 'B12']
    temp_dir = os.path.join(params['output_dir'], 'temp')
    os.makedirs(temp_dir, exist_ok=True)

    tiles = [params['sentinel']['tiles']]

    for tile in tiles:

        # Locate the reference Sentinel-2 image
        sentinel_images = [
            f for f in glob.glob(
                os.path.join(params['output_dir'], 'atmcor', 'sentinel', tile, '**', '*.SAFE*', '*_B*.tif')
            )
            if any(band in f for band in sentinel_bands)
        ]
        if not sentinel_images:
            raise FileNotFoundError(f"No Sentinel-2 images found for tile: {tile}")

        sentinel_img = sentinel_images[0]

        # Set Sentinel tile projection and shapefile
        params['sen2_epsg_code'] = satwutils.raster2meta(sentinel_img)
        params['sen_tile_target_shp'] = satwutils.get_tile_shp(tile, params, params['sen2_epsg_code'])

        # Create output directories
        params['output_dir_tiling'] = os.path.join(params['output_dir'], 'tiling', tile)
        satwutils.create_dir(os.path.join(params['output_dir_tiling'], 'sentinel'))

        # Identify Sentinel-2 scenes for processing
        sentinel_scene_dir = os.path.join(params['output_dir'], 'atmcor', 'sentinel', tile)
        sentinel_scenes = [os.path.join(sentinel_scene_dir, scene) for scene in os.listdir(sentinel_scene_dir)]

        # Run resampling in parallel
        with multiprocessing.Pool(processes=params['aux_info']['n_cores']) as pool:
            results = pool.starmap_async(gen_resample, [(scene, params) for scene in sentinel_scenes]).get()
            print(f"Processing results for tile {tile}: {results}")