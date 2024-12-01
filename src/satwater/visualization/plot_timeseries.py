
import os
import matplotlib.pyplot as plt
import numpy as np
import glob
import rasterio
import shutil

import glob
from osgeo import gdal
import matplotlib
matplotlib.use('Agg')  # this "mode" does not show the results and allows the sequential plots
from matplotlib import pyplot as plt
from matplotlib import colors
import multiprocessing
import numpy as np
import os
from osgeo import gdal


def raster2meta(geotif_file):
    import rasterio

    # Open the GeoTIFF file
    with rasterio.open(geotif_file) as src:
        proj = src.crs
        # Get the EPSG code
        epsg_code = int(proj.to_epsg())

    return epsg_code

def epsg2utmzone(epsg_code):
    import pyproj
     # Define the EPSG identifier
    epsg_code = f'EPSG:{epsg_code}' # Replace with your own EPSG identifier
     # Extract the UTM zone from the EPSG identifier
    epsg_num = int(epsg_code.split(':')[1])
    zone = (epsg_num - 32600) % 100
     # Print the resulting UTM zone
    print(f'UTM zone: {zone}')
    return zone

def getimgsize(path_geotiff):
    import rasterio
    # Open the raster file in read mode
    with rasterio.open(path_geotiff, 'r') as src:
        # Get the dimensions of the raster
        columns = src.width
        rows = src.height
        # Print the dimensions of the raster
    return rows, columns

import pyproj
import rasterio
import datetime
img_path = "path/to/img.tif"

central_lon_wgs = -88.68
central_lat_wgs = 30.32
epsg_code = raster2meta(img_path)  # get the EPSG code from the GeoTIFF file
zone = epsg2utmzone(epsg_code)  # get the UTM zone from the EPSG code

# Define the WGS84 and target UTM projections
wgs84 = pyproj.Proj(proj='latlong', datum='WGS84')
utm_target = pyproj.Proj(proj='utm', zone=zone, ellps='WGS84', north=True)
# Convert the WGS84 coordinates to the target UTM zone
lon, lat = (central_lon_wgs, central_lat_wgs)  # Replace with your own WGS84 coordinates
central_long_utm, central_lat_utm = pyproj.transform(wgs84, utm_target, central_lon_wgs, central_lat_wgs)

# Open the geotiff file
with rasterio.open(img_path) as src:
    # Extract the row and column pixel indices corresponding to the (x, y) coordinates
    y_row, x_col = src.index(central_long_utm, central_lat_utm)

    # Print the resulting row and column indices
    print(f"Row: {y_row}, Column: {x_col}")


if __name__ == '__main__':

    # what satellite you want to process?
    period_ini = '20130101'
    period_end = '20231201'
    hls_dir = r'V:\processed\funded_projects\nasarid2023\out_satwater\hlswater'
    output_dir = r'V:\processed\funded_projects\nasarid2023\out_satwater\plots'
    tiles = '16RCU'

    period = (period_ini, period_end)
    path_ = fr'{hls_dir}\{tiles}'
    path_lists = os.listdir(path_)

    hlspath_img = glob.glob(fr'{hls_dir}\{tiles}\{path_lists[0]}\*B04.tif', recursive=True)[0]
    rows, columns = getimgsize(hlspath_img)

    dates = []
    stack_time = np.zeros((rows, columns, len(path_lists)), dtype=np.float32)
    for i, hlspath in enumerate(path_lists):
        hlspath_img = glob.glob(fr'{hls_dir}\{tiles}\{hlspath}\*B04.tif', recursive=True)[0]
        dates.append(hlspath.split('.')[1].split('T')[0])
        with rasterio.open(hlspath_img, 'r') as src:
            raster_array = src.read()
            stack_time[:,:,i] = raster_array


    # Extract the time series of the specified row and column
    time_series = stack_time[y_row, x_col, :]

    time_series=np.where(time_series == -9999, np.nan, time_series)
    time_series/= 10000.0

     # Create a list of date strings
     # Convert the date strings to datetime objects
    dates_dt = [datetime.strptime(date, '%Y%m%d').date() for date in dates]
     # Create the plot
    plt.plot(dates_dt, time_series)
    plt.title('Time Series with Dates')
    plt.xlabel('Date')
    plt.ylabel('Value')
    plt.show()



