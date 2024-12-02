import os
import glob
import rasterio
import matplotlib
import numpy as np
from osgeo import gdal
import multiprocessing
from matplotlib import pyplot as plt

matplotlib.use('Agg')  # this "mode" does not show the results and allows the sequential plots

def getimgsize(path_geotiff):

    # Open the raster file in read mode
    with rasterio.open(path_geotiff, 'r') as src:
        columns = src.width
        rows = src.height
    return rows, columns

def plot3B(image_b, image_g, image_r, output_png):

    ds = gdal.Open(image_b)
    arr_b = ds.ReadAsArray().astype(np.int16)  # uint8
    ds = gdal.Open(image_g)
    arr_g = ds.ReadAsArray().astype(np.int16)  # uint8
    ds = gdal.Open(image_r)
    arr_r = ds.ReadAsArray().astype(np.int16)  # uint8

    arr_b = np.where(arr_b == -9999, np.nan, arr_b)
    arr_g = np.where(arr_g == -9999, np.nan, arr_g)
    arr_r = np.where(arr_r == -9999, np.nan, arr_r)

    try:
        stacked_array = np.stack((arr_b, arr_g, arr_r), axis=2)
    except:
        min_height = min(arr_b.shape[0], arr_g.shape[0], arr_r.shape[0])
        min_width = min(arr_b.shape[1], arr_g.shape[1], arr_r.shape[1])

        arr_b_clipped = arr_b[:min_height, :min_width]
        arr_g_clipped = arr_g[:min_height, :min_width]
        arr_r_clipped = arr_r[:min_height, :min_width]

        stacked_array = np.stack((arr_b_clipped, arr_g_clipped, arr_r_clipped), axis=2)

    stacked_array/=10000.0

    brightness = 5.0
    stacked_array = np.clip(stacked_array * brightness, 0.0, 1.0)

    plt.ioff()
    fig = plt.figure()
    fig.set_size_inches(8.5, 6.5)

    plt.imshow(stacked_array, cmap='gist_earth')

    tile = os.path.basename(output_png)

    plt.title(tile.split(".")[1] + '_' + tile)
    plt.tight_layout()
    plt.xticks([])
    plt.yticks([])

    plt.savefig(output_png, bbox_inches='tight', dpi=250)
    plt.close()

def run_plots(path_hls, hlsscene, output_dir_png):

    image_b = glob.glob(fr'{path_hls}\{hlsscene}\*B02.tif', recursive=True)[0]
    image_g = glob.glob(fr'{path_hls}\{hlsscene}\*B03.tif', recursive=True)[0]
    image_r = glob.glob(fr'{path_hls}\{hlsscene}\*B04.tif', recursive=True)[0]

    output_png = rf'{output_dir_png}\{hlsscene}.png'

    plot3B(image_b, image_g, image_r, output_png)


def run(params):

    sat = params['aux_info']['sat_name']

    if sat == 'landsat':
        tiles = os.listdir(fr'{params["output_dir_tiling"]}\landsat')
    else:
        tiles = params['sentinel']['tiles']

    for tile in tiles:

        output_dir_png = fr"{params['output_dir']}\plots\{tile}"
        os.makedirs(output_dir_png, exist_ok=True)

        path_hls = fr"{params['output_dir']}\hlswater\{tile}"
        n_hlsscene = os.listdir(path_hls)

        n_path_hls = [path_hls]*len(n_hlsscene)
        n_output_dir_png = [output_dir_png]*len(n_hlsscene)

        with multiprocessing.Pool(processes=4) as pool:
            results = pool.starmap_async(run_plots, zip(n_path_hls, n_hlsscene, n_output_dir_png)).get()
            print(results)
            pool.close()