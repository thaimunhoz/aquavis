# Libraries
import ee
import os
import glob
import geemap
import rasterio
import numpy as np
import rioxarray as rxr
import geopandas as gpd
from rasterio.mask import mask
from rasterio.merge import merge
from shapely.ops import transform
from shapely.geometry import shape
from rasterio.features import shapes
from rasterio.warp import calculate_default_transform, reproject, Resampling

ee.Initialize(project = "ee-thaimunhoz98")

'''
Functions to generate water mask from Landsat and Sentinel images
'''

def raster2shp(raster_numpy, raster_transform, raster_crs):

    results = (
        {'properties': {'raster_val': v}, 'geometry': s}
        for s, v in shapes(raster_numpy.astype(np.int16), transform=raster_transform)
        if v  # Only take shapes with raster_val = True (i.e., v=1)
    )

    geometries = [shape(feature['geometry']) for feature in results]

    gdf = gpd.GeoDataFrame(geometry=geometries, crs=raster_crs)

    return gdf

def get_wm_jrc(gdf_input, percentage_occurrence, output_crs):

    '''
    Generates the water mask using the JRC Global Surface Water dataset
    Input:
        gdf_input (GeoDataFrame): input shapefile
        percentage_occurrence (int): percentage of water occurrence
        output_path (str): path to save the output shapefile
        output_crs (int): output crs
    Return:
        str: path to the water mask shapefile
    '''

    def to_2d(geom):
        return transform(lambda x, y, _: (x, y), geom)  # Drops the Z value

    gdf_input.loc[:, "geometry"] = gdf_input["geometry"].apply(to_2d)

    feature_collection = geemap.geopandas_to_ee(gdf_input)

    gsw = ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
    occurrence = gsw.select('occurrence')
    water_mask = occurrence.gt(percentage_occurrence)

    water_mask_bbox = water_mask.clip(feature_collection).reproject('EPSG:' + str(output_crs), None, 30)

    # Save the JRC water mask as .tiff file
    geemap.ee_export_image(water_mask_bbox, filename = 'water_mask_jrc.tif', region=feature_collection.geometry().bounds())

    try:
        with rasterio.open('water_mask_jrc.tif') as src:
            wm_jrc = src.read(1)
    except:
        jrc_img_path = get_jrc_gcer(gdf_input, percentage_occurrence)

        with rasterio.open(jrc_img_path) as src:
            wm_jrc = src.read(1)

    gdf = raster2shp(wm_jrc, src.transform, src.crs)

    # list all files in the directory with "water_mask" in the name and remove them all
    files = glob.glob('*water_mask*')
    for f in files: os.remove(f)

    return gdf

def create_water_mask(image_path, bbox_lake):

    '''
    Generate the water mask using the Modified NDWI and the JRC Water Mask
    Input:
        image_path (str): path to the image folder
        bbox_lake (GeoDataFrame): bounding box of the lake
        output_path (str): path to save the output shapefile
    Return:
        str: path to the modified water mask shapefile
        int: number of pixels
    '''

    # JRC Water Mask
    # Read the NIR band to get the image crs
    with rasterio.open((glob.glob(os.path.join(image_path, '*B04*')) + glob.glob(os.path.join(image_path, '*B4*')))[0]) as _band: nir_transform = _band.transform

    gdf_jrc = get_wm_jrc(bbox_lake, 75, _band.crs.to_epsg())

    return gdf_jrc

def get_jrc_gcer(bbox_geom, percentage_occurrence):

    image_path = r'Z:\guser\tml\mypapers\hls_synthetic\JRC_DATA'

    files = [f for f in os.listdir(image_path) if os.path.isfile(os.path.join(image_path, f)) and f.endswith('.tif')]

    count = 0 # in case of more than one image intersect the lake
    images_list = []

    for file in files:

        with rasterio.open(image_path + '/' + file) as src:

            raster_bounds = src.bounds

            geom = bbox_geom.geometry

            bounds_aux = bbox_geom.bounds

            if (bounds_aux[['minx']].values > raster_bounds[2] or bounds_aux[['maxx']].values < raster_bounds[0] or bounds_aux[['miny']].values > raster_bounds[3] or bounds_aux[['maxx']].values < raster_bounds[1]):

                continue

            else:

                # Read the subset of the raster data within the polygon area
                try:
                    subset, out_transform = mask(src, geom, nodata=src.nodata, crop=True, indexes=1)
                except:
                    continue

                subset_clip = subset > percentage_occurrence

                out_meta = src.meta.copy()

                out_meta.update({"driver": "GTiff",
                                 "height": subset_clip.shape[0],
                                 "width": subset_clip.shape[1],
                                 "transform": out_transform,
                                 "nodata": 0})

                with rasterio.open('water_mask_jrc_' + str(count) + '.tif', "w", **out_meta) as dest:
                    dest.write(subset_clip.astype(int), 1)
                    dest.nodata = 0

                images_list.append('water_mask_jrc_' + str(count) + '.tif')
                count += 1

    if count > 1:

        src_files_to_mosaic = []

        for fp in images_list:
            src = rasterio.open(fp)
            src_files_to_mosaic.append(src)

        # Merge the tiffs into one mosaic
        mosaic, out_trans = merge(src_files_to_mosaic)

        out_meta = src_files_to_mosaic[0].meta.copy()

        out_meta.update({
                    "driver": "GTiff",
                    "height": mosaic.shape[1],
                    "width": mosaic.shape[2],
                    "transform": out_trans,
        })

        with rasterio.open(output_path + '/' + 'water_mask_jrc.tif', "w", **out_meta) as dest:
            dest.write(mosaic)

        return output_path + '/' + 'water_mask_jrc.tif'

    else:

        return output_path + '/' + 'water_mask_jrc_0.tif'