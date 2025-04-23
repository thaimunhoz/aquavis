import os
import glob
import rasterio
import matplotlib
import numpy as np
import xarray as xr
import multiprocessing
import rioxarray as rxr
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap

matplotlib.use('Agg')  # this "mode" does not show the results and allows the sequential plots

def getimgsize(path_geotiff):

    # Open the raster file in read mode
    with rasterio.open(path_geotiff, 'r') as src:
        columns = src.width
        rows = src.height
    return rows, columns

def contrast_stretch(image, low_perc=0.5, high_perc=99.5):
    """
    Perform contrast stretching on an image using specified percentiles.

    Parameters:
    - image: xarray.DataArray, the image to be adjusted
    - low_perc: float, the lower percentile to clip
    - high_perc: float, the upper percentile to clip

    Returns:
    - contrast_stretched_image: xarray.DataArray, the contrast-stretched image
    """

    # Compute the lower and upper percentiles
    p_low = np.percentile(image, low_perc)
    p_high = np.percentile(image, high_perc)

    # Stretch the contrast
    contrast_stretched_image = (image - p_low) / (p_high - p_low)
    contrast_stretched_image = np.clip(contrast_stretched_image, 0, 1)

    return contrast_stretched_image

def plot3B(image_b, image_g, image_r, output_png):

    arr_b = rxr.open_rasterio(image_b)/10000.0
    arr_g = rxr.open_rasterio(image_g)/10000.0
    arr_r = rxr.open_rasterio(image_r)/10000.0

    try:
       stacked_array = np.stack((arr_b, arr_g, arr_r))

    except:
       min_height = min(arr_b.shape[0], arr_g.shape[0], arr_r.shape[0])
       min_width = min(arr_b.shape[1], arr_g.shape[1], arr_r.shape[1])

       arr_b_clipped = arr_b[:min_height, :min_width]
       arr_g_clipped = arr_g[:min_height, :min_width]
       arr_r_clipped = arr_r[:min_height, :min_width]

       stacked_array = np.stack((arr_b_clipped, arr_g_clipped, arr_r_clipped))

    brightness = 5.0
    stacked_array = np.clip(stacked_array * brightness, 0.0, 1.0)
    stacked_array = np.squeeze(stacked_array)
    stacked_array = np.transpose(stacked_array, (1, 2, 0)) if stacked_array.shape[0] == 3 else stacked_array

    cmap = LinearSegmentedColormap.from_list('bright colormap', [(0, 0, 0), (1, 1, 1)], N=256)

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

    # rgb = xr.concat([arr_r, arr_g, arr_b], dim='band')
    #
    # # Normalize the data to [0, 1] for plotting
    # rgb_scaled = (rgb * (255 / rgb.max().values)).astype(np.uint8)
    #
    # rgb_contrast = xr.concat([contrast_stretch(rgb_scaled[0]), contrast_stretch(rgb_scaled[1]), contrast_stretch(rgb_scaled[2])], dim='band')
    #
    # plt.ioff()
    # fig = plt.figure()
    # fig.set_size_inches(8.5, 6.5)
    # plt.imshow(rgb_contrast.transpose('y', 'x', 'band'))
    # plt.axis('off')
    # tile = os.path.basename(output_png)
    # plt.title(tile.split(".")[1] + '_' + tile)
    # plt.tight_layout()
    # plt.xticks([])
    # plt.yticks([])
    #
    # plt.savefig(output_png, bbox_inches='tight', dpi=250)
    # plt.close()

def run_plots(path_hls, hlsscene, output_dir_png):

    image_b = glob.glob(fr'{path_hls}\{hlsscene}\*B02.tif', recursive=True)[0]
    image_g = glob.glob(fr'{path_hls}\{hlsscene}\*B03.tif', recursive=True)[0]
    image_r = glob.glob(fr'{path_hls}\{hlsscene}\*B04.tif', recursive=True)[0]

    output_png = rf'{output_dir_png}\{hlsscene}.png'

    plot3B(image_b, image_g, image_r, output_png)

def run(params):

    output_dir_png = fr"{params['output_dir']}\plots"
    os.makedirs(output_dir_png, exist_ok=True)

    path_hls = fr"{params['output_dir']}\AQUAVis_Product"
    n_hlsscene = os.listdir(path_hls)

    n_path_hls = [path_hls]*len(n_hlsscene)
    n_output_dir_png = [output_dir_png]*len(n_hlsscene)

    with multiprocessing.Pool(processes=4) as pool:
        results = pool.starmap_async(run_plots, zip(n_path_hls, n_hlsscene, n_output_dir_png)).get()
        print(results)
        pool.close()