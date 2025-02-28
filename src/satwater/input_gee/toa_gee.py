import ee
import geemap
import xarray as xr

ee.Initialize(project = "ee-thaimunhoz98")

import logging

class IgnoreSpecificMessageFilter(logging.Filter):
    def filter(self, record):
        # Ignore messages containing "Connection pool is full"
        return "Connection pool is full" not in record.getMessage()

# Apply the filter to the root logger
logger = logging.getLogger()
logger.addFilter(IgnoreSpecificMessageFilter())

def get_available_dates_sentinel(tile, period_ini, period_end):
    """
    Retrieves available dates for Sentinel-2 within a given date range and cloud threshold.

    Args:
    - tile: str, tile in the format "path_row" (e.g., "089_120")
    - period_ini: str, start date in the format "YYYY-MM-DD"
    - period_end: str, end date in the format "YYYY-MM-DD"
    - cloud_threshold: int, minimum cloud coverage percentage (default is 30)

    Returns:
    - list: Available dates within the period with cloud coverage below the threshold
    """

    # Define Landsat 8 and Landsat 9 collections
    sentinel_coll = ee.ImageCollection("COPERNICUS/S2_HARMONIZED") \
        .filter(ee.Filter.eq("MGRS_TILE", tile)) \
        .filterDate(period_ini, period_end) \
        .filter(ee.Filter.lt("CLOUDY_PIXEL_PERCENTAGE", 30))

    # Extract dates and convert to a list
    available_images = sentinel_coll.aggregate_array('PRODUCT_ID').getInfo()

    return available_images

def get_available_dates_landsat(tile, period_ini, period_end, cloud_threshold=30):
    """
    Retrieves available dates for Landsat 8 and 9 within a given date range and cloud threshold.

    Args:
    - tile: str, tile in the format "path_row" (e.g., "089_120")
    - period_ini: str, start date in the format "YYYY-MM-DD"
    - period_end: str, end date in the format "YYYY-MM-DD"
    - cloud_threshold: int, minimum cloud coverage percentage (default is 30)

    Returns:
    - list: Available dates within the period with cloud coverage below the threshold
    """

    # Split tile input into path and row
    path, row = map(int, tile.split('_'))

    #period_ini = ee.Date(period_ini)
    #period_end = ee.Date(period_end)

    # Define Landsat 8 and Landsat 9 collections
    landsat_8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_TOA") \
        .filter(ee.Filter.eq("WRS_PATH", path)) \
        .filter(ee.Filter.eq("WRS_ROW", row)) \
        .filterDate(period_ini, period_end) \
        .filter(ee.Filter.lt("CLOUD_COVER", 30))

    landsat_9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_TOA") \
        .filter(ee.Filter.eq("WRS_PATH", path)) \
        .filter(ee.Filter.eq("WRS_ROW", row)) \
        .filterDate(period_ini, period_end) \
        .filter(ee.Filter.lt("CLOUD_COVER", 30))

    # Merge both collections
    merged_collection = landsat_8.merge(landsat_9)

    # Extract dates and convert to a list
    available_images = merged_collection.aggregate_array('LANDSAT_PRODUCT_ID').getInfo()

    # Convert to readable date format
    #date_list = [ee.Date(date).format('YYYY-MM-dd').getInfo() for date in available_dates]

    return available_images

def sentinel_toa_from_gee(img_id):

    """
    Fetches an image from GEE, converts it to a DataArray, and includes metadata.

    Parameters:
        tile (ee.Geometry): The region of interest (geometry for the tile).
        date (str): The date for the image (format: 'YYYY-MM-DD').
        collection_name (str): GEE collection name ('LANDSAT/LC08/C02/T1_L2', 'COPERNICUS/S2', etc.).

    Returns:
        xr.DataArray: Image data as a stacked DataArray with metadata.
    """

    # start_date = ee.Date(date)
    # end_date = start_date.advance(1, 'days')
    #
    # path, row = map(int, path_row.split("_"))

    # Load the collection and filter
    sentinel_collection = (
        ee.ImageCollection("COPERNICUS/S2_HARMONIZED")
        .filter(ee.Filter.eq("PRODUCT_ID", img_id))
    )

    # Retrieve metadata information
    metadata_dict = sentinel_collection.getInfo()

    if sentinel_collection.size().getInfo() == 0:
        print(f"No images found for ID {img_id}.")
        return None

    image = sentinel_collection.first()

    selected_bands = ["B2", "B3", "B4", "B8A", "B11"]

    # Return the region of interest
    geom = sentinel_collection.first().select(1).geometry().bounds().getInfo()['coordinates'][0]
    longitudes = [coord[0] for coord in geom]
    latitudes = [coord[1] for coord in geom]
    bounding_box = [min(longitudes), min(latitudes), max(longitudes), max(latitudes)]

    images_per_band_10 = []
    images_per_band_20 = []
    for i in selected_bands:

        image_vnir = image.select(i)
        if i == "B8A" or i == "B11":
            scale = 20
            band_xda = geemap.ee_to_xarray(image_vnir, projection=sentinel_collection.first().select(i).projection(), scale=scale, geometry=bounding_box).transpose("time", "Y", "X")
            #band_xda.assign_coords(band=[i])
            images_per_band_20.append(band_xda)
        else:
            scale = 10
            band_xda = geemap.ee_to_xarray(image_vnir, projection=sentinel_collection.first().select(i).projection(), scale = scale, geometry = bounding_box).transpose("time", "Y", "X")
            #band_xda = band_xda.assign_coords(band=[i])
            images_per_band_10.append(band_xda)

    #vnir_xda_10 = xr.concat(images_per_band_10, dim="band")
    #vnir_xda_20 = xr.concat(images_per_band_20, dim="band")
    vnir_xda_10 = xr.merge(images_per_band_10)
    vnir_xda_20 = xr.merge(images_per_band_20)

    return vnir_xda_10, vnir_xda_20, metadata_dict

def landsat_toa_from_gee(img_id):

    """
    Fetches an image from GEE, converts it to a DataArray, and includes metadata.

    Parameters:
        tile (ee.Geometry): The region of interest (geometry for the tile).
        date (str): The date for the image (format: 'YYYY-MM-DD').
        collection_name (str): GEE collection name ('LANDSAT/LC08/C02/T1_L2', 'COPERNICUS/S2', etc.).

    Returns:
        xr.DataArray: Image data as a stacked DataArray with metadata.
    """

    # start_date = ee.Date(date)
    # end_date = start_date.advance(1, 'days')
    #
    # path, row = map(int, path_row.split("_"))

    # Load the collection and filter
    landsat8_collection = (
        ee.ImageCollection("LANDSAT/LC08/C02/T1_TOA")
        .filter(ee.Filter.eq("LANDSAT_PRODUCT_ID", img_id))
    )

    landsat9_collection = (
        ee.ImageCollection("LANDSAT/LC09/C02/T1_TOA")
        .filter(ee.Filter.eq("LANDSAT_PRODUCT_ID", img_id))
    )

    # Merge both collections
    combined_collection = landsat8_collection.merge(landsat9_collection)

    # Retrieve metadata information
    metadata_dict = combined_collection.getInfo()

    if combined_collection.size().getInfo() == 0:
        print(f"No images found for ID {img_id}.")
        return None

    image = combined_collection.first()

    # Get common band names
    #band_names = image.bandNames().getInfo()
    selected_bands = ["B2", "B3", "B4", "B5", "B6"]

    image_vnir = image.select(selected_bands)

    # Return the region of interest
    geom = combined_collection.first().select(0).geometry().bounds().getInfo()['coordinates'][0]
    longitudes = [coord[0] for coord in geom]
    latitudes = [coord[1] for coord in geom]
    bounding_box = [min(longitudes), min(latitudes), max(longitudes), max(latitudes)]

    vnir_xda = geemap.ee_to_xarray(image_vnir, projection=combined_collection.first().select(0).projection(), scale = 30, geometry = bounding_box).transpose("time", "Y", "X")

    # Extract angle data
    angle_bands = ["SAA", "SZA", "VAA", "VZA"]
    image_angles = image.select(angle_bands)
    angles_xda = geemap.ee_to_xarray(image_angles, projection=combined_collection.first().select(0).projection(), scale = 30, geometry = bounding_box).transpose("time", "Y", "X")

    # Applying the scale and offset factors
    # Apply condition and scaling to all bands
    for band in angles_xda.data_vars:
        angles_xda[band] = angles_xda[band].where(angles_xda[band] != -32768, float('nan'))
        angles_xda[band] = angles_xda[band].where(angles_xda[band] != 32768, float('nan'))
        angles_xda[band] = angles_xda[band] / 100  # Apply scale factor (if necessary)

    return vnir_xda, angles_xda, metadata_dict

# Example Usage
# tile = "21HVB"
# period_ini = "2023-01-01"
# period_end = "2023-12-31"
#
# available_images = get_available_dates_sentinel(tile, period_ini, period_end)
# print(available_images)