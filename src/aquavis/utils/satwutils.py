import os
import shutil
import os.path
import rasterio
import rioxarray
import numpy as np
import xarray as xr
import rioxarray as rxr
import geopandas as gpd
from shapely.geometry import shape
from rasterio.features import shapes
from shapely.geometry import mapping
from typing import List, Dict, Any, Optional
from scipy.ndimage import binary_dilation as bn
from rasterio.warp import calculate_default_transform, reproject, Resampling

def raster2shp(raster_numpy, raster_transform, raster_crs):

    results = list(
        {"properties": {"raster_val": v}, "geometry": s}
        for s, v in shapes(np.asarray(raster_numpy, dtype=np.int16), transform=raster_transform)
        if v  # Only take shapes with raster_val = True (i.e., v=1)
    )
    geometries = [shape(feature["geometry"]) for feature in results]

    gdf = gpd.GeoDataFrame(geometry=geometries, crs=raster_crs)

    # Keep only the polygons with the biggest area
    # gdf['area'] = gdf['geometry'].area
    # gdf_largest = gdf.sort_values(by='area', ascending=False).head(1)

    return gdf

def buffer_into_polygon(gdf, buffer_size):

    gdf['geometry'] = gdf.geometry.buffer(buffer_size)
    gdf_cleaned = gdf[~gdf['geometry'].is_empty & gdf['geometry'].is_valid]

    # Keep only the polygons with the biggest area
    # gdf_cleaned['area'] = gdf_cleaned['geometry'].area
    # gdf_largest = gdf_cleaned.sort_values(by='area', ascending=False).head(1)

    return gdf_cleaned

def create_inward_buffer(input_raster_path, buffer_distance_meters=3000, pixel_resolution=30):
    band_path = [i for i in os.listdir(input_raster_path)]
    green_band = next((band for band in band_path if "B3" in band or "B03" in band), None)
    swir_band = next((band for band in band_path if "B11" in band or "B6" in band or "B06" in band), None)

    xda_green = rxr.open_rasterio(os.path.join(input_raster_path, green_band))
    xda_swir = rxr.open_rasterio(os.path.join(input_raster_path, swir_band))

    xda_green_matched = xda_green.rio.reproject_match(xda_swir)

    mndwi = (xda_green_matched - xda_swir) / (xda_green_matched + xda_swir)
    mndwi_mask = mndwi > 0.2
    water_mask = mndwi_mask.astype(int).values[0, :, :]

    # Calculate buffer in pixels (buffer distance converted to pixels)
    buffer_pixels = int(abs(buffer_distance_meters) / pixel_resolution)

    # Invert the water mask (0 -> water, 1 -> non-water)
    inverted_mask = np.where(water_mask == 1, 0, 1)

    # Apply dilation to the inverted mask to expand the non-water area (expand the outer non-water body)
    dilated_mask = bn(inverted_mask, structure=np.ones((buffer_pixels, buffer_pixels)))

    # Invert back the dilated mask to create the inward buffer (1 for inside the buffer, 0 for outside)
    buffer_mask = np.where(np.logical_and(dilated_mask == 1, water_mask == 1), 1, 0)

    return

def get_mask_water(img_path):

    band_path = [i for i in os.listdir(img_path)]
    green_band = next((band for band in band_path if "B3" in band or "B03" in band), None)
    swir_band = next((band for band in band_path if "B11" in band or "B6" in band or "B06" in band), None)

    xda_green = rxr.open_rasterio(os.path.join(img_path, green_band))
    xda_swir = rxr.open_rasterio(os.path.join(img_path, swir_band))

    xda_green_matched = xda_green.rio.reproject_match(xda_swir)

    mndwi = (xda_green_matched - xda_swir) / (xda_green_matched + xda_swir)
    mndwi_mask = mndwi > 0.2
    water_mask = mndwi_mask

    with rasterio.open(os.path.join(img_path, green_band)) as src:
        _transform = src.transform
        _crs = src.crs

    # Include adjacency buffer
    water_shp = raster2shp(water_mask, _transform, _crs)
    water_copy = water_shp.copy()

    water_shp_buffered = buffer_into_polygon(water_copy, -3000)

    # Get only the difference between the original and the buffered
    water_shp_final = gpd.GeoDataFrame(geometry=water_shp.difference(water_shp_buffered), crs=water_shp.crs)

    # water_shp_final_copy = water_shp_final.copy()
    # water_adjc = buffer_into_polygon(water_shp_final_copy, 500)

    return water_shp_final

def raster2meta(geotif_file):

    with rasterio.open(geotif_file) as src:

        proj = src.crs
        epsg_code = int(proj.to_epsg())

    return epsg_code

def resampling_landsat(src_file, dst_file, ref_path):

    with rasterio.open(ref_path) as ref_src:
        ref_crs = ref_src.crs  # Coordinate reference system of the reference
        ref_transform = ref_src.transform

    # Open the target image and reproject
    with rasterio.open(src_file) as target_src:

        # Calculate the transform and shape for the reprojected image
        transform, width, height = calculate_default_transform(
            target_src.crs, ref_crs, target_src.width, target_src.height, *target_src.bounds
        )

        # Define metadata for the output file
        kwargs = target_src.meta.copy()
        kwargs.update({
            'crs': ref_crs,
            'transform': transform
        })

        # Reproject and save the reprojected image
        with rasterio.open(dst_file, 'w', **kwargs) as dst:
            for i in range(1, target_src.count + 1):  # Loop over each band
                reproject(
                    source=rasterio.band(target_src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=target_src.transform,
                    src_crs=target_src.crs,
                    dst_transform=transform,
                    dst_crs=ref_crs,
                    resampling=Resampling.cubic  # Use cubic resampling
                )

def xml_to_dict(element):

    if len(element) == 0:
        return element.text

    result = {}

    for child in element:

        child_dict = xml_to_dict(child)

        if child.tag in result:
            if isinstance(result[child.tag], list):
                result[child.tag].append(child_dict)
            else:
                result[child.tag] = [result[child.tag], child_dict]
        else:
            result[child.tag] = child_dict

    return result

def _setup_paths(input_path: str, output_path: str) -> None:
    """Initialize paths for processing."""
    image_path = input_path
    output_path = output_path
    os.makedirs(output_path, exist_ok=True)

def _copy_metadata_file(image_path, output_path) -> None:
    """Copy metadata file to output directory."""
    shutil.copy(os.path.join(image_path, "MTD.xml"),
                os.path.join(output_path, "MTD.xml"))

def _prepare_image_data(src: rasterio.DatasetReader) -> np.ndarray:
    """Prepare image data by handling nodata values."""
    image_data = src.read(1)
    return np.where(np.isin(image_data, [-9999, 0]), np.nan, image_data)

def _save_corrected_image(data: np.ndarray, crs: Any, transform: Any, nodata: Any, path: str) -> None:
    """Save corrected image to GeoTIFF."""
    with rasterio.open(
            path,
            'w',
            driver='GTiff',
            count=1,
            dtype=data.dtype,
            width=data.shape[1],
            height=data.shape[0],
            crs=crs,
            transform=transform,
            compress='lzw',
            nodata=nodata
    ) as dst:
        dst.write(data, 1)

def reproject(src_file, dst_file, epsg_num):

    """
    Reproject a raster file to a new coordinate reference system.

    Args:
        src_file (str): Path to the source raster file.
        dst_file (str): Path to the destination raster file.
        epsg_num (int): EPSG code of the target CRS.
    """

    #data = xr.open_dataarray(src_file)
    data = rioxarray.open_rasterio(src_file)
    data = data.rio.write_crs(data.rio.crs)  # Ensure CRS is written if not present

    # Reproject to the target CRS
    reprojected_data = data.rio.reproject(f"EPSG:{epsg_num}", resampling=Resampling.cubic)
    # try:
    #     reprojected_data.rio.to_raster(dst_file)
    # except:
    #     os.remove(dst_file)
    #     reprojected_data.rio.to_raster(dst_file)

    return reprojected_data

def cut_images_res(data, shapefile, path_output, spatialres):

    # Load the raster as an xarray DataArray
    #data = rioxarray.open_rasterio(path_original)
    data = data.rio.write_crs(data.rio.crs)

    # Clip the raster using the shapefile geometry
    geometries = [mapping(geom) for geom in shapefile.geometry]
    clipped_data = data.rio.clip(geometries, shapefile.crs, drop=True)

    # Get the bounds of the Sentinel tile
    minx, miny, maxx, maxy = shapefile.total_bounds

    # Create a target transform
    transform = rasterio.transform.from_origin(minx, maxy, spatialres, spatialres)
    out_shape = (int((maxy - miny) / spatialres), int((maxx - minx) / spatialres))

    # Reproject to match the tile extent
    resampled_data = clipped_data.rio.reproject(
        dst_crs=shapefile.crs,
        shape=out_shape,
        transform=transform,
        resampling=Resampling.bilinear
    )

    # Create an empty array with the Sentinel tile's dimensions (with NoData values)
    empty_array = np.full((1, *out_shape), -9999, dtype=resampled_data.dtype)

    # Compute valid indices for insertion
    affine_transform = resampled_data.rio.transform()
    clipped_x_start = max(0, int((resampled_data.x[0] - minx) / spatialres))
    clipped_y_start = max(0, int((maxy - resampled_data.y[0]) / spatialres))

    # Determine end indices
    clipped_x_end = min(out_shape[1], clipped_x_start + resampled_data.shape[2])
    clipped_y_end = min(out_shape[0], clipped_y_start + resampled_data.shape[1])

    # Insert data into the empty array
    empty_array[:, clipped_y_start:clipped_y_end, clipped_x_start:clipped_x_end] = \
        resampled_data.values[:, :clipped_y_end - clipped_y_start, :clipped_x_end - clipped_x_start]

    empty_array[empty_array <= 0] = -9999

    # Convert back to DataArray
    result = xr.DataArray(
        empty_array,
        dims=resampled_data.dims,
        coords=resampled_data.coords,
        attrs=resampled_data.attrs
    )

    result.rio.write_crs(shapefile.crs, inplace=True)

    # Save the final image
    result.rio.to_raster(path_output)

def create_dir(path):

    if not os.path.exists(path):

        os.makedirs(path)

def intersection(sen_tile_target, params):

    # Read in the shapefiles for Sentinel and Landsat
    df_sentinel = gpd.read_file(params['sentinel']['tiles_shp'], encoding="UTF-8")
    df_landsat = gpd.read_file(params['landsat']['tiles_shp'], encoding="UTF-8")

     # Set the index of the Sentinel dataframe to 'Name'
    gdf_sentinel = df_sentinel.set_index('Name')

     # Get the geometry of the target Sentinel tile
    target_tile = gdf_sentinel.loc[sen_tile_target]['geometry']

     # Find the intersecting Landsat tiles
    res_intersection = df_landsat[df_landsat.intersects(target_tile)]

     # Get the pathrow numbers of the intersecting Landsat tiles
    pathrows = [f"{i[0]:03d}_{i[1]:03d}" for i in zip(res_intersection['PATH'].to_list(), res_intersection['ROW'].to_list())]

    return pathrows

def intersection_landsat(params):

    # Read in the shapefiles for Sentinel and Landsat
    df_sentinel = gpd.read_file(params['sentinel']['tiles_shp'], encoding="UTF-8")
    df_landsat = gpd.read_file(params['landsat']['tiles_shp'], encoding="UTF-8")

    cleaned_data = [f"{int(part.split('_')[0]):d}_{int(part.split('_')[1]):d}" for part in params['landsat']['tiles']]
    target_landsat_tile = df_landsat[df_landsat['Name'].isin(cleaned_data)]

    intersection_aux = gpd.overlay(df_sentinel, target_landsat_tile, how='intersection').rename(columns={'Name_1': 'Name'})

    return intersection_aux

def get_tile_shp(sen_tile_target, params, epsg_code):

    df_sentinel = gpd.read_file(params['sentinel']['tiles_shp'], encoding="UTF-8")
    df_select = df_sentinel[df_sentinel['Name'] == sen_tile_target]

    # Specify the desired output projection
    output_crs = f'EPSG:{epsg_code}'  # For example, this is the Web Mercator projection
    gdf_reprojected = df_select.to_crs(output_crs)

    #save shapefile
    shpfile = fr"{params['output_dir']}\temp\temp_{sen_tile_target}.shp"
    gdf_reprojected.to_file(shpfile)

    return shpfile

def get_band_name(band_nm, sat):

    '''
    Band name from HLS
    '''

    band_nm=int(band_nm)

    if sat=='landsat':

        if band_nm >= 430 and band_nm <= 450:
            band_nmi = 'B01'  # coastal
        elif band_nm >= 450 and band_nm <= 510:
            band_nmi = 'B02'   # blue
        elif band_nm >= 530 and band_nm <= 590:
            band_nmi = 'B03'     # green
        elif band_nm >= 640 and band_nm <= 670:
            band_nmi = 'B04'     # red
        elif band_nm >= 850 and band_nm <= 880:
            band_nmi = 'B05'     # nir narrow
        elif band_nm >= 1570 and band_nm <= 1650:
            band_nmi = 'B06'     # SWIR1
        elif band_nm >= 2110 and band_nm <= 2290:
            band_nmi = 'B07'  # swir2
        elif band_nm >= 591 and band_nm <= 595:
            band_nmi = 'PAN'  # swir2
        else:
            raise Exception("Band name not found")

    else: #sat == 'sentinel':

        if band_nm >= 430 and band_nm <= 450:
            band_nmi = 'B01' # coastal
        elif band_nm >= 450 and band_nm <= 510:
            band_nmi = 'B02' # blue
        elif band_nm >= 530 and band_nm <= 590:
            band_nmi = 'B03' # green
        elif band_nm >= 640 and band_nm <= 670:
            band_nmi = 'B04' # red
        elif band_nm >= 690 and band_nm <= 710:
            band_nmi = 'B05' #rededge1
        elif band_nm >= 730 and band_nm <= 750:
            band_nmi = 'B06' #rededge2
        elif band_nm >= 770 and band_nm <= 790:
            band_nmi = 'B07' #rededge3
        elif band_nm >= 820 and band_nm <= 845:
            band_nmi = 'B08' #nir broad
        elif band_nm >= 850 and band_nm <= 880:
            band_nmi = 'B8A' #nir narrow
        elif band_nm >= 1570 and band_nm <= 1650:
            band_nmi = 'B11' # swir1
        elif band_nm >= 2110 and band_nm <= 2290:
            band_nmi = 'B12' # swir1
        else:
            raise Exception("Band name not found")

    return band_nmi

def fix_band_name(band_nm, sat):

    '''
    Band name from HLS
    '''

    if sat=='landsat':

        if band_nm == 'B2': band_nmi = 'B02'   # blue
        elif band_nm == 'B3': band_nmi = 'B03' # green
        elif band_nm == 'B4': band_nmi = 'B04' # red
        elif band_nm == 'B5': band_nmi = 'B05' # nir narrow
        elif band_nm == 'B6': band_nmi = 'B06' # SWIR1
        elif band_nm == 'B7': band_nmi = 'B07' # swir2
        else:
            band_nmi = band_nm
    else:
        band_nmi = band_nm

    return band_nmi