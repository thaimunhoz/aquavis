import os
import glob
import os.path
import multiprocessing
from satwater.utils import satwutils

#sys.path.append(r'Z:\projects\funded\nasarid\sat-water')

def gen_tiles(landsat_scene, params):

    landsat_bands = ['B2', 'B3', 'B4', 'B5']

    # Tiles do Sentinel-2 que interceptam o tile do Landsat. Para cada Sentinel tile, a gente precisa clipar e reprojetar a imagem do Landsat.
    sen_tile = params['sen_tiles_target']

    landsat_scene_all_bands = [f for f in glob.glob(os.path.join(landsat_scene, '*')) if any(band in f for band in landsat_bands)]

    for sen_tile_target in sen_tile['Name']: # Loop through the Sentinel tiles

        # First we need to check if we have the image tile in our database
        #sen_tile_path = glob.glob(fr'{params["sentinel"]["input_dir"]}\{sen_tile_target}')
        #if sen_tile_path == []:
        #    continue

        #sen_tile_ref = glob.glob(fr'{sen_tile_path[0]}\{os.listdir(sen_tile_path[0])[0]}\**\GRANULE\**\IMG_DATA\*_B04*.jp2')

        sen_gdf = sen_tile[sen_tile['Name'] == sen_tile_target]

        sen2_epsg_code = sen_gdf['ESPG'].values[0]

        for landsat_band in landsat_scene_all_bands:

            base_dir_nm = fr"{params['output_dir_tiling']}\landsat\{sen_tile_target}\{os.path.basename(os.path.dirname(landsat_band))}"
            satwutils.create_dir(base_dir_nm)
            out_band = fr"{base_dir_nm}\{os.path.basename(landsat_band)}"

            if os.path.exists(out_band):
                continue

            landsat_epsg_code = satwutils.raster2meta(landsat_band)

            imgtemp = fr"{params['output_dir']}\temp\temp_{os.path.basename(landsat_band)}"
            satwutils.reproject(landsat_band, imgtemp, sen2_epsg_code)  # reprojecting the landsat to sentienl projection

            #imgtemp2 = fr"{params['output_dir']}\temp\temp2_{os.path.basename(landsat_band)}"
            #satwutils.resampling_landsat(imgtemp, imgtemp2, sen_tile_ref[0])

            satwutils.cut_images_res(imgtemp, sen_gdf, out_band, 30)

            #os.remove(imgtemp)

def run(params):

    os.makedirs(fr"{params['output_dir']}\temp", exist_ok=True)

    sentinel_pathrows = satwutils.intersection_landsat(params)

    params['sen_tiles_target'] = sentinel_pathrows

    # Create the output directory if it does not exist
    params['output_dir_tiling'] = fr'{params["output_dir"]}\tiling'
    satwutils.create_dir(fr'{params["output_dir_tiling"]}\landsat')

    pathrows = [f"{i[0]:03d}_{i[1]:03d}" for i in zip(sentinel_pathrows['PATH'].to_list(), sentinel_pathrows['ROW'].to_list())]
    landsat_pathrows = list(set(pathrows))

    for landsat_pathrow in landsat_pathrows:

        landsat_scenes = glob.glob(fr'{params["output_dir"]}\atmcor\landsat\{landsat_pathrow}\**\*LC*')

        n_params = [params]*len(landsat_scenes)

        with multiprocessing.Pool(processes=params['aux_info']['n_cores']) as pool:
            results = pool.starmap_async(gen_tiles, zip(landsat_scenes, n_params)).get()
            print(results)
            pool.close()