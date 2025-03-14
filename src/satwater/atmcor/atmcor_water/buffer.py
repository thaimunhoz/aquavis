# -*- mode: python -*

import os
import glob
import rasterio
import numpy as np
import geopandas as gpd
from osgeo import osr, gdal, ogr
gdal.UseExceptions()
import rasterio
import toolbox as tool

from shapely.geometry import shape
from rasterio.features import shapes

os.environ['GDAL_DRIVER_PATH'] = r'C:\Users\tml411\AppData\Local\anaconda3\pkgs\libgdal-jp2openjpeg-3.9.2-hed4c6cb_1\Library\lib\gdalplugins'

class Buffer_OLCI_S3:

    def __init__(self,
                 path_roi,
                 dest):

        self.path_roi = path_roi
        self.dest = dest
        self.WMASKID = 'wmask'
        self.BUFFER = str("nan")

    def run(self):
        """
        Returns a buffer (10-km) from water boundary with only land cover.
        """
        # Water mask -- Oa06 and 0a21:
        arr = {}
        band = [i for i in glob.glob(os.path.join(self.dest, '*.tif')) if 'Oa06' in i or 'Oa21' in i]
        try:
            for num, i in enumerate(band):
                arr[num] = tool.loadarray(i)
            wmask = self.loadwmask(arr, 0.05) # default --threshold equal to 0.05.
            tool.export(wmask, self.WMASKID + 'xxxx', band[0], self.dest)
        except:
            pass
        # Buffer -- 50km -- from water boundary:
        buffer10km = self.buffer(self.path_roi, self.dest)
        # Crops the water mask to ROI:
        cut = self.cutbands(buffer10km, self.dest + '/' + self.WMASKID + '.tif', self.WMASKID + 'xroi', self.dest)
        # Water mask from .tifF to .SHP:
        wmask_shp = self.rasterToshapefile(cut, self.dest)
        # Selects the land cover and export the NEW BUFFER:
        wmask_shp = gpd.read_file(wmask_shp)
        get_land = wmask_shp.loc[wmask_shp['data'] == 0]
        merge_ = get_land.dissolve(by='data')
        self.BUFFER = self.dest + '/buffer.shp'
        merge_.to_file(self.BUFFER)

    def loadwmask(self, array: dict, threshold: float) -> dict:
            """
            Returns the water mask based on spectral index MNDWI.
            """
            mndwi = (array[0] - array[1]) / (array[0] + array[1])
            return np.where(mndwi >= threshold, 1, 0) # default --threshold equal to 0.05.

    def rasterToshapefile(self, raster, dest: str) -> str:
            """
            Converts raster (.tifF) to shapefile (.shp).
            """
            src_ds = gdal.Open(raster)
            srcband = src_ds.GetRasterBand(1)
            dst_layername = 'atm/data'
            drv = ogr.GetDriverByName("ESRI Shapefile")
            dst_ds = drv.CreateDataSource(dest + '/wmask.shp')
            sp_ref = osr.SpatialReference()
            sp_ref.SetFromUserInput('EPSG: 4326')
            dst_layer = dst_ds.CreateLayer(dst_layername, srs=sp_ref)
            fld = ogr.FieldDefn("data", ogr.OFTInteger)
            dst_layer.CreateField(fld)
            dst_field = dst_layer.GetLayerDefn().GetFieldIndex("data")
            gdal.Polygonize(srcband, None, dst_layer, dst_field, [], callback=None)
            return dest + '/wmask.shp'

    def buffer(self, shapefile: str, dest: str) -> str:
            """
            Generates a buffer of 10km from waterbody edge.
            """
            geometry = gpd.read_file(shapefile)
            buffer_distance = 10 / 111  # degree # default
            buffered_geometry = geometry.copy()
            buffered_geometry['geometry'] = geometry.buffer(buffer_distance)
            buffered_geometry.to_file(dest + '/roi_buffer.shp')
            return dest + '/roi_buffer.shp'

    def cutbands(self, path_shapefile: str, path_image: str, index: str, dest: str) -> str:
        """
        Crops the images.
        """
        kwargs = {'cutlineDSName': 'True', 'dstNodata': 'np.nan', '-to': 'Float32'}
        gdal.Warp(dest + '/' + index + '.tif', path_image,
                  cutlineDSName=path_shapefile,
                  cropToCutline=True,
                  dstNodata=-9999)
        return dest + '/' + index + '.tif'


class Buffer_MSI_S2:

    def __init__(self,
                 path_roi,
                 dest):

        self.path_roi = path_roi
        self.dest = dest
        self.WMASKID = 'wmask'
        self.BUFFER = str("nan")

    def run(self):

        """
        Returns a buffer (50-km) from water boundary with only land cover.
        """

        # Water mask -- B3 and B12:
        arr = {}

        self.resample([i for i in glob.glob(os.path.join(self.dest, '*.tif')) if 'B12' in i][0], 10)

        band = [i for i in glob.glob(os.path.join(self.dest, '*.tif')) if 'B03' in i or 'B12_resamp' in i]

        for num, i in enumerate(band):
            with rasterio.open(i) as src:
                out_meta = src.meta.copy()
                out_transform = src.transform
                arr[num] = src.read(1)

        wmask = self.loadwmask(arr, 0.05) # default --threshold equal to 0.05.

        out_meta.update({"driver": "GTiff",
                         "height": wmask.shape[0],
                         "width": wmask.shape[1],
                         "transform": out_transform,
                         "dtype": 'int32',
                         })

        index = self.WMASKID + 'xxxx'

        with rasterio.open(self.dest + '/' + index[0:-4] + '.tif', "w", **out_meta) as dest:
            dest.write(wmask, 1)

        #tool.export(wmask, self.WMASKID + 'xxxx', band[0], self.dest)

        # Buffer -- 50km -- from water boundary:
        src = rasterio.open(band[0])
        buffer10km = self.buffer(self.path_roi, self.dest, src.crs)

        # Crops the water mask to ROI:
        cut = self.cutbands(buffer10km, self.dest + '/' + self.WMASKID + '.tif', self.WMASKID + 'xroi', self.dest)
        self.resample(cut, 300)

        # Water mask from .tifF to .SHP:
        wmask_shp = self.rasterToshapefile(cut[:-4] + '_resamp.tif', self.dest, src.crs['init'])

        # Selects the land cover and export the NEW BUFFER:
        wmask_shp = gpd.read_file(wmask_shp)
        wmask_shp = wmask_shp.to_crs(epsg=4326)
        get_land = wmask_shp.loc[wmask_shp['data'] == 0]
        merge_ = get_land.dissolve(by='data')
        if len(merge_) == 0:
            get_land = wmask_shp.loc[wmask_shp['data'] == 1]
            merge_ = get_land.dissolve(by='data')
        else:
            pass
        self.BUFFER = self.dest + '/buffer.shp'
        merge_.to_file(self.BUFFER)

    def loadwmask(self, array: dict, threshold: float) -> dict:
            """
            Returns the water mask based on spectral index MNDWI.
            """
            mndwi = (array[0] - array[1]) / (array[0] + array[1])
            return np.where(mndwi >= threshold, 1, 0) # default --threshold equal to 0.05.

    def rasterToshapefile(self, raster, dest: str, crs) -> str:
            """
            Converts raster (.tifF) to shapefile (.shp).
            """
            src_ds = gdal.Open(raster)
            srcband = src_ds.GetRasterBand(1)
            dst_layername = 'atm/data'
            drv = ogr.GetDriverByName("ESRI Shapefile")
            dst_ds = drv.CreateDataSource(dest + '/wmask.shp')
            sp_ref = osr.SpatialReference()
            # sp_ref.SetFromUserInput('EPSG: 4326')
            sp_ref.SetFromUserInput(crs)
            dst_layer = dst_ds.CreateLayer(dst_layername, srs=sp_ref)
            fld = ogr.FieldDefn("data", ogr.OFTInteger)
            dst_layer.CreateField(fld)
            dst_field = dst_layer.GetLayerDefn().GetFieldIndex("data")
            gdal.Polygonize(srcband, None, dst_layer, dst_field, [], callback=None)
            return dest + '/wmask.shp'

    def buffer (self, shapefile: str, dest: str, crs) -> str:
            """
            Generates a buffer of 50km from waterbody edge.
            """
            geometry = gpd.read_file(shapefile)
            geometry = geometry.to_crs(crs=crs) # geotransform.
            buffer_distance = 10000  # meter # default
            buffered_geometry = geometry.copy()
            buffered_geometry['geometry'] = geometry.buffer(buffer_distance)
            buffered_geometry.to_file(dest + '/roi_buffer.shp')
            return dest + '/roi_buffer.shp'

    def cutbands(self, path_shapefile: str, path_image: str, index: str, dest: str) -> str:
        """
        Crops the images.
        """
        kwargs = {'cutlineDSName': 'True', 'dstNodata': 'np.nan', '-to': 'Float32'}
        gdal.Warp(dest + '/' + index + '.tif', path_image,
                  cutlineDSName=path_shapefile,
                  cropToCutline=True,
                  dstNodata=-9999)
        return dest + '/' + index + '.tif'

    def resample(self, filename: str, res: int) -> None:
        """
        Re-samples the pixel size to 10 m.
        """
        gdal.Warp(filename[:-4] + '_resamp.tif', filename, xRes=res, yRes=res, resampleAlg='bilinear')
        return None


class Buffer_OLI_L89:

    def __init__(self,
                 path_main,
                 path_roi,
                 dest):

        self.path_main = path_main
        self.path_roi = path_roi
        self.dest = dest
        self.WMASKID = 'wmask'
        self.BUFFER = str("nan")

    def run(self):
        """
        Returns a buffer (50-km) from water boundary with only land cover.
        """
        # Water mask -- B3 and B7:
        arr = {}
        band = [i for i in glob.glob(os.path.join(self.path_main, '*.tif')) if 'B3' in i or 'B7' in i]
        try:
            for num, i in enumerate(band):
                #arr[num] = tool.loadarray(i)
                with rasterio.open(i) as src:
                    out_meta = src.meta.copy()
                    out_transform = src.transform
                    arr[num] = src.read(1)

            wmask = self.loadwmask(arr, 0.05) # default --threshold equal to 0.05.

            out_meta.update({"driver": "GTiff",
                             "height": wmask.shape[0],
                             "width": wmask.shape[1],
                             "transform": out_transform,
                             "dtype": 'float32',
                             })

            #wmask_uint16 = wmask.astype('uint16')

            index = self.WMASKID + 'xxxx'
            with rasterio.open(self.dest + '/' + index[0:-4] + '.tif', "w", **out_meta) as dest:
                dest.write(wmask, 1)
            #tool.export(wmask, self.WMASKID + 'xxxx', band[0], self.dest)

        except:
            pass
        # Buffer -- 10km -- from water boundary:
        src = rasterio.open(band[0])
        buffer10km = self.buffer(self.path_roi, self.dest, src.crs)
        # Crops the water mask to ROI:
        cut = self.cutbands(buffer10km, self.dest + '/' + self.WMASKID + '.tif', self.WMASKID + 'xroi', self.dest)
        # Water mask from .tifF to .SHP:
        wmask_shp = self.rasterToshapefile(cut, self.dest, src.crs['init'])
        # Selects the land cover and export the NEW BUFFER:
        #wmask_shp = gpd.read_file(wmask_shp)
        wmask_shp = wmask_shp.to_crs(epsg=4326)
        get_land = wmask_shp
        #get_land = wmask_shp.loc[wmask_shp['data'] == 0]
        merge_ = get_land.dissolve()#.drop(columns=['FID'])
        self.BUFFER = self.dest + '/buffer.shp'
        merge_.to_file(self.BUFFER)

    def loadwmask(self, array: dict, threshold: float) -> dict:
        """
        Returns the water mask based on spectral index MNDWI.
        """
        mndwi = (array[0] - array[1]) / (array[0] + array[1])
        return np.where(mndwi >= threshold, 1, 0) # default --threshold equal to 0.05.

    def rasterToshapefile(self, raster, dest: str, crs) -> str:
            """
            Converts raster (.tifF) to shapefile (.shp).
            """
            with rasterio.open(raster) as src:
                image_arr = src.read(1)
                ref = src.crs
                src_transform = src.transform

            #mask = image_arr == 0

            # Convert mask to polygon shapes
            results = (
                {'properties': {'raster_val': v}, 'geometry': s}
                for s, v in shapes(image_arr.astype(np.int16), transform=src_transform)
                if v == 0  # Only take shapes with raster_val = True (i.e., v=1)
            )

            geometries = [shape(feature['geometry']) for feature in results]

            gdf = gpd.GeoDataFrame(geometry=geometries, crs=ref)

            gdf_reprojected = gpd.GeoDataFrame(geometry=gdf.geometry.set_crs(ref, allow_override=True), crs=ref)

            gdf_reprojected.to_file(dest + '/wmask.shp')

            #src_ds = gdal.Open(raster)
            #srcband = src_ds.GetRasterBand(1)
            #dst_layername = 'data'
            #drv = ogr.GetDriverByName("ESRI Shapefile")
            #dst_ds = drv.CreateDataSource(dest + '/wmask.shp')
            #sp_ref = osr.SpatialReference()

            ## sp_ref.SetFromUserInput('EPSG: 4326')
            #sp_ref.SetFromUserInput(crs)
            #dst_layer = dst_ds.CreateLayer(dst_layername, srs=sp_ref)
            #fld = ogr.FieldDefn("data", ogr.OFTInteger)
            #dst_layer.CreateField(fld)
            #dst_field = dst_layer.GetLayerDefn().GetFieldIndex("data")
            #gdal.Polygonize(srcband, None, dst_layer, dst_field, [], callback=None)

            return gdf_reprojected

    def buffer (self, shapefile: str, dest: str, crs) -> str:
            """
            Generates a buffer of 10km from waterbody edge.
            """
            geometry = gpd.read_file(shapefile)
            geometry = geometry.to_crs(crs=crs) # geotransform.
            buffer_distance = 10000  # meter # default
            buffered_geometry = geometry.copy()
            buffered_geometry['geometry'] = geometry.buffer(buffer_distance)
            buffered_geometry.to_file(dest + '/roi_buffer.shp')
            return dest + '/roi_buffer.shp'

    def cutbands(self, path_shapefile: str, path_image: str, index: str, dest: str) -> str:
        """
        Crops the images.
        """
        kwargs = {'cutlineDSName': 'True', 'dstNodata': 'np.nan', '-to': 'Float32'}
        gdal.Warp(dest + '/' + index + '.tif', path_image,
                  cutlineDSName=path_shapefile,
                  cropToCutline=True,
                  dstNodata=-9999)
        return dest + '/' + index + '.tif'
