import os
import os.path
import rasterio
import rioxarray
import numpy as np
import xarray as xr
import rioxarray as rxr
import geopandas as gpd
from shapely.geometry import mapping
from rasterio.warp import calculate_default_transform, reproject, Resampling

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
    try:
        reprojected_data.rio.to_raster(dst_file)
    except:
        os.remove(dst_file)
        reprojected_data.rio.to_raster(dst_file)

def cut_images_res(path_original, shapefile_tile, path_output, spatialres):

    try:
        shapefile = gpd.read_file(shapefile_tile)
    except:
        shapefile = shapefile_tile

    # Load the raster as an xarray DataArray
    data = rioxarray.open_rasterio(path_original)
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