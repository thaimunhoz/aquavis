import os
import glob
import calendar
import numpy as np
import pandas as pd
import rioxarray as rxr
from src.satwater.atmcor.atmcor_water import toolbox as tool
from datetime import datetime, timedelta
from src.satwater.atmcor.atmcor_water.atm.coefficient import MCDExtractWindow

class Metadata_MSI_S2:

    def __init__(self,
                 path_main: str,
                 path_dest: str,
                 networkdrive_letter: str,
                 satellite: str,
                 aero_type: str,
                 mode=None,
                 msi_tile=None):

        self.MTD = '/MTD_TL.xml'
        self.BAND_ID = '_B'
        self.GRANULE = '/GRANULE'
        self.IMG_DATA = '/IMG_DATA'
        self.MTD_MSIL1C = '/MTD_MSIL1C.xml'
        self.MOD08_D3 = ':/dbcenter/products/atm/modis/C61/MOD08_D3'
        self.MDE = ':/dbcenter/products/land/dem30m'
        self.TEMP_COEF = ':/public/temp_dir'
        self.TEMP = path_dest + '/' + path_main[-65:] + '/tempdir'

        self.path_main = path_main
        self.path_dest = path_dest
        self.networkdrive_letter = networkdrive_letter
        self.satellite = satellite
        self.aero_type = aero_type
        self.mode = mode
        self.msi_tile = msi_tile

        self.type = str("nan")
        self.bandname = list("nan")
        self.s2path = str("nan")
        self.aod = float("nan")
        self.water_vapour = float("nan")
        self.ozone = float("nan")
        self.altitude = float("nan")
        self.geometry = {}
        self.datetime = {}
        self.rescale = {}

    def run(self):

        """
        Scans the metadata.
        """

        path = [i for i in glob.glob(os.path.join(self.path_main + self.GRANULE, '*')) if 'L1C_' in i]

        self.rescale_factor()

        self.s2path = path[0] + self.IMG_DATA

        # Return the bouding box of the image:
        self.roi = tool.return_water(self.s2path, self.rescale, self.msi_tile)

        if len(self.roi) == 0:
            self.roi = tool.return_bbox(os.path.join(self.s2path, os.listdir(self.s2path)[3]))

        # Defining the preferred order of spectral bands, including multiple naming conventions for the same band.
        bands_order = ['B01', 'B02', 'B03', 'B04','B05', 'B06','B07','B08','B8A', 'B09', 'B10', 'B11', 'B12']

        aux_bandnames = [i for i in os.listdir(self.s2path) if self.BAND_ID in i]
        aux_bandnames.insert(8, aux_bandnames.pop(-1))  # band name

        # ensure selected bands are sorted in the order defined in bands_order
        self.bandname = sorted(aux_bandnames, key=lambda x: next((i for i, band in enumerate(bands_order) if band in x), float('inf')))

        self.dict_metadata = tool.xml_to_json(str(path[0]) + self.MTD) # metadata from sensor
        self.type = str(self.dict_metadata['n1:Level-1C_Tile_ID']['n1:General_Info']['TILE_ID']['#text'][0:3]) # safe number A or B
        self.date_and_time()
        self.geo()
        self.read_coefficient()

        df = pd.DataFrame({'img': [self.path_main], 'aod': [self.aod], 'wv': [self.water_vapour], 'oz': [self.ozone], 'alt': [self.altitude]})
        df.to_csv(self.path_dest + '/' + 'atm_parameters.csv')

    def date_and_time(self):

        """
        Returns the date and the time.
        """

        # Verifies the metadata:
        date_and_time = self.dict_metadata["n1:Level-1C_Tile_ID"]['n1:General_Info']['SENSING_TIME']['#text']
        date_acquired = date_and_time[0:10]
        date = datetime.strptime(date_acquired, '%Y-%m-%d').timetuple()
        scene_center_time = date_and_time[11:-1].split(':')
        time_hh = int(scene_center_time[0]) + (float(scene_center_time[1]) / 60) + (float(scene_center_time[2]) / 3600)

        # Export date and time in a metadata structure:
        value = DateTime()
        value.day = date.tm_mday
        value.month = date.tm_mon
        value.year = date.tm_year
        value.time_hh = time_hh

        self.datetime = value

    def geo(self):

        """
        Returns the geometry of observation and illumination.
        """

        # Sun azimuth [az] and zenith [zn] angle:
        # View azimuth [az] and zenith [zn] angles:
        # It considers the angle average of the scene.

        for i in range(0, 13):

            output = {}
            output['solar_az'] = float(self.dict_metadata['n1:Level-1C_Tile_ID']['n1:Geometric_Info']['Tile_Angles']['Mean_Sun_Angle']['AZIMUTH_ANGLE']['#text'])
            output['solar_zn'] = float(self.dict_metadata['n1:Level-1C_Tile_ID']['n1:Geometric_Info']['Tile_Angles']['Mean_Sun_Angle']['ZENITH_ANGLE']['#text'])
            output['view_az'] = np.array(float(self.dict_metadata['n1:Level-1C_Tile_ID']['n1:Geometric_Info']['Tile_Angles']["Mean_Viewing_Incidence_Angle_List"]['Mean_Viewing_Incidence_Angle'][i]['AZIMUTH_ANGLE']['#text'])).mean()
            output['view_zn'] = np.array(float(self.dict_metadata['n1:Level-1C_Tile_ID']['n1:Geometric_Info']['Tile_Angles']["Mean_Viewing_Incidence_Angle_List"]['Mean_Viewing_Incidence_Angle'][i]['ZENITH_ANGLE']['#text'])).mean()

            self.geometry[i] = output

    def rescale_factor(self):

        """
        Returns the factor values to convert from DN_TOA to REFLECTANCE_TOA
        """

        MSIL1C = tool.xml_to_json(self.path_main + self.MTD_MSIL1C)
        QUANTIFICATION_VALUE = float(MSIL1C['n1:Level-1C_User_Product']['n1:General_Info']['Product_Image_Characteristics']['QUANTIFICATION_VALUE']['#text'])

        for i in range(0, 13):

            try:
                RADIO_ADD_OFFSET = float(MSIL1C['n1:Level-1C_User_Product']['n1:General_Info']['Product_Image_Characteristics']['Radiometric_Offset_List']['RADIO_ADD_OFFSET'][i]['#text'])
            except:
                RADIO_ADD_OFFSET = 0

            self.rescale[i] = {'qvalue': QUANTIFICATION_VALUE, 'offset': RADIO_ADD_OFFSET}

    def read_coefficient(self):

        """
        Recovers the atmospheric coefficients from MODIS (MOD08_D3) and SRTM.
        """

        # Average values are obtained based on a time window (in days). The maximum window is 7-days.
        # The values are returned from the time window size closer the target date.

        start_date = str(self.datetime.year) + '-' + str(self.datetime.month) + '-' + str(self.datetime.day)
        end_date = str(self.datetime.year) + '-' + str(self.datetime.month) + '-' + str(self.datetime.day)

        mcd_scanner = MCDExtractWindow(dir_mod08=str(self.networkdrive_letter) + self.MOD08_D3,
                                        dir_mde=str(self.networkdrive_letter) + self.MDE,
                                        dir_temp=str(self.networkdrive_letter) + self.TEMP_COEF,
                                        ini_date=start_date,
                                        end_date=end_date,
                                        bounding_shp=self.roi)

        dataset_info_mod08 = mcd_scanner.run_extraction_mod08d3()
        dataset_info_mde = mcd_scanner.run_extract_mde()

        self.aod = dataset_info_mod08['AOD_mean'].mean()
        self.water_vapour = dataset_info_mod08['WV_mean'].mean()
        self.ozone = dataset_info_mod08['OZ_mean'].mean() / 1000 # in cm_atm
        self.altitude = dataset_info_mde['MDE_mean'].mean() / 1000 # in km

        # Weekly mean values:
        if (np.isnan(self.aod) or self.aod == 0.0) or (np.isnan(self.water_vapour) or self.water_vapour == 0.0) or (np.isnan(self.ozone) or self.ozone == 0.0):

            date = datetime.strptime(start_date, '%Y-%m-%d')

            start_of_week = (date - timedelta(days=date.weekday()))
            end_of_week = (start_of_week + timedelta(days=6))

            start_of_week_str = start_of_week.strftime('%Y-%m-%d')
            end_of_week_str = end_of_week.strftime('%Y-%m-%d')

            mcd_scanner = MCDExtractWindow(dir_mod08=str(self.networkdrive_letter) + self.MOD08_D3,
                                            dir_mde=str(self.networkdrive_letter) + self.MDE,
                                            dir_temp=str(self.networkdrive_letter) + self.TEMP_COEF,
                                            ini_date=start_of_week_str,
                                            end_date=end_of_week_str,
                                            bounding_shp=self.roi)

            dataset_info_mod08 = mcd_scanner.run_extraction_mod08d3()

            self.aod = dataset_info_mod08['AOD_mean'].mean()
            self.water_vapour = dataset_info_mod08['WV_mean'].mean()
            self.ozone = dataset_info_mod08['OZ_mean'].mean() / 1000  # in cm_atm

        # Montly mean values:
        if (np.isnan(self.aod) or self.aod == 0.0) or (np.isnan(self.water_vapour) or self.water_vapour == 0.0) or (np.isnan(self.ozone) or self.ozone == 0.0):

            month_name = calendar.month_name[self.datetime.month]

            gdf = self.roi

            self.aod = mcd_scanner.get_modis_monthly_mean(month_name, 'AOD', gdf)
            self.water_vapour = mcd_scanner.get_modis_monthly_mean(month_name, 'Water_Vapor', gdf)
            self.ozone = mcd_scanner.get_modis_monthly_mean(month_name, 'Total_Ozone', gdf) / 1000

class Metadata_OLI_L89:

    def __init__(self,
                 path_main: str,
                 path_dest: str,
                 networkdrive_letter: str,
                 satellite: str,
                 aero_type: str,
                 mode=None,
                 msi_tile=None):

        self.MTD = '/MTD_TL.xml'
        self.BAND_ID = '_B'
        self.MTD_ID = 'MTL'
        self.ANG_ID = 'ANG'
        self.MOD08_D3 = ':/dbcenter/products/atm/modis/C61/MOD08_D3'
        self.MDE = ':/dbcenter/products/land/dem30m'
        self.TEMP_COEF = ':/public/temp_dir'
        self.TEMP = path_dest + '/' + path_main[-40:] + '/tempdir'

        self.path_main = path_main
        self.path_dest = path_dest
        self.networkdrive_letter = networkdrive_letter
        self.satellite = satellite
        self.aero_type = aero_type
        self.mode = mode
        self.msi_tile = msi_tile

        self.type = str("nan")
        self.bandname = list("nan")
        self.s2path = str("nan")
        self.aod = float("nan")
        self.water_vapour = float("nan")
        self.ozone = float("nan")
        self.altitude = float("nan")
        self.geometry = {}
        self.datetime = {}
        self.rescale = {}

    def run(self):

        """
        Scans the metadata.
        """

        # Defining the preferred order of spectral bands, including multiple naming conventions for the same band.
        bands_order = ['B1', 'B2','B3', 'B4','B5','B6', 'B7','B8','B9','B10','B11','B12']

        aux_bandnames = [i for i in os.listdir(self.path_main) if
                         self.BAND_ID in i and 'B9' not in i and 'B10' not in i and 'B11' not in i and 'B12' not in i]

        # ensure selected bands are sorted in the order defined in bands_order
        self.bandname = sorted(aux_bandnames, key=lambda x: next((i for i, band in enumerate(bands_order) if band in x), float('inf')))

        path = [i for i in glob.glob(os.path.join(self.path_main, '*.xml')) if self.MTD_ID in i]

        self.dict_metadata = tool.xml_to_json(str(path[0]))  # metadata from sensor
        self.rescale_factor()

        # Return the bouding box of the image:
        self.roi = tool.return_water(self.path_main, self.rescale, self.msi_tile)

        if len(self.roi) == 0:
            self.roi = tool.return_bbox(glob.glob(os.path.join(self.path_main, self.bandname[0]))[0])

        #self.roi = tool.return_bbox(glob.glob(os.path.join(self.path_main, self.bandname[0]))[0])
        #self.roi = tool.return_bbox(self.path_main)

        self.type = str(self.dict_metadata['LANDSAT_METADATA_FILE']['PRODUCT_CONTENTS']['LANDSAT_PRODUCT_ID'][0:4]) # safe number L8 or L9
        self.date_and_time()
        self.geo()
        self.read_coefficient()

        os.makedirs(self.path_dest, exist_ok = True)
        df = pd.DataFrame({'img': [self.path_main], 'aod': [self.aod], 'wv': [self.water_vapour], 'oz': [self.ozone], 'alt': [self.altitude]})
        df.to_csv(self.path_dest + '/' + 'atm_parameters.csv')

    def date_and_time(self):

        """
        Returns the date and the time.
        """

        # Verifies the metadata:
        date_acquired = self.dict_metadata['LANDSAT_METADATA_FILE']['IMAGE_ATTRIBUTES']['DATE_ACQUIRED']
        date = datetime.strptime(date_acquired, '%Y-%m-%d').timetuple()
        scene_center_time = self.dict_metadata['LANDSAT_METADATA_FILE']['IMAGE_ATTRIBUTES']['SCENE_CENTER_TIME'][0:16].split(':')
        time_hh = int(scene_center_time[0]) + (float(scene_center_time[1]) / 60) + (float(scene_center_time[2]) / 3600)

        # Export date and time in a metadata structure:
        value = DateTime()
        value.day = date.tm_mday
        value.month = date.tm_mon
        value.year = date.tm_year
        value.time_hh = time_hh
        self.datetime = value

    def geo(self):

        """
        Returns the geometry of observation and illumination.
        """

        # Sun azimuth [az] and zenith [zn] angle:
        # View azimuth [az] and zenith [zn] angles:
        # It considers the angle averages of the scene.
        # The band 4 (red) is used as reference because it is near the center of the OLI/Landsat-8/9 focal plane.

        solar_az = rxr.open_rasterio(glob.glob(os.path.join(self.path_main, '*SAA.tif'))[0]).values.astype(float)
        solar_zen = rxr.open_rasterio(glob.glob(os.path.join(self.path_main, '*SZA.tif'))[0]).values.astype(float)
        view_az = rxr.open_rasterio(glob.glob(os.path.join(self.path_main, '*VAA.tif'))[0]).values.astype(float)
        view_zen = rxr.open_rasterio(glob.glob(os.path.join(self.path_main, '*VZA.tif'))[0]).values.astype(float)

        solar_az[(solar_az == 0) | (solar_az == -9999)] = np.nan
        solar_zen[(solar_zen == 0) | (solar_zen == -9999)] = np.nan
        view_az[(view_az == 0) | (view_az == -9999)] = np.nan
        view_zen[(view_zen == 0) | (view_zen == -9999)] = np.nan

        for i in range(0, 8):

            output = {}

            output['solar_az'] = np.mean(solar_az[~np.isnan(solar_az)]) / 100
            output['solar_zn'] = np.mean(solar_zen[~np.isnan(solar_zen)]) / 100
            output['view_az'] = np.mean(view_az[~np.isnan(view_az)]) / 100
            output['view_zn'] = np.mean(view_zen[~np.isnan(view_zen)]) / 100

            self.geometry[i] = output

    def rescale_factor(self):

        """
        Returns the factor values to convert from DN_TOA to REFLECTANCE_TOA
        """

        # It obtains the factor values to convert from DN_TOA to REFLECTANCE_TOA:
        for i in range(1, 9):
            ADD_BAND = float(self.dict_metadata['LANDSAT_METADATA_FILE']['LEVEL1_RADIOMETRIC_RESCALING']['REFLECTANCE_ADD_BAND_' + str(i)])
            MULT_BAND = float(self.dict_metadata['LANDSAT_METADATA_FILE']['LEVEL1_RADIOMETRIC_RESCALING']['REFLECTANCE_MULT_BAND_' + str(i)])
            self.rescale[i - 1] = {'add': ADD_BAND, 'mult': MULT_BAND}

    def read_coefficient(self):

        """
        Recovers the atmospheric coefficients from MODIS (MOD08_D3) and SRTM.
        """

        # Atmospheric parameters:
        # Average values are obtained based on a time window (in days). The maximum window is 7-days.
        # The values are returned from the time window size closer the target date.

        start_date = str(self.datetime.year) + '-' + str(self.datetime.month) + '-' + str(self.datetime.day)
        end_date = str(self.datetime.year) + '-' + str(self.datetime.month) + '-' + str(self.datetime.day)

        mcd_scanner = MCDExtractWindow(dir_mod08=str(self.networkdrive_letter) + self.MOD08_D3,
                                           dir_mde=str(self.networkdrive_letter) + self.MDE,
                                           dir_temp=str(self.networkdrive_letter) + self.TEMP_COEF,
                                           ini_date=start_date,
                                           end_date=end_date,
                                           bounding_shp=self.roi)

        # dataset_info_mcd19a2 = mcd_scanner.run_extraction_mcd19a2()
        dataset_info_mod08 = mcd_scanner.run_extraction_mod08d3()
        dataset_info_mde = mcd_scanner.run_extract_mde()

        self.aod = float(dataset_info_mod08['AOD_mean'].mean())
        self.water_vapour = dataset_info_mod08['WV_mean'].mean()
        self.ozone = dataset_info_mod08['OZ_mean'].mean() / 1000  # in cm_atm
        self.altitude = dataset_info_mde['MDE_mean'].mean() / 1000  # in km

        # Weekly mean values:
        if (np.isnan(self.aod) or self.aod == 0.0) or (np.isnan(self.water_vapour) or self.water_vapour == 0.0) or (np.isnan(self.ozone) or self.ozone == 0.0):
            date = datetime.strptime(start_date, '%Y-%m-%d')

            start_of_week = (date - timedelta(days=date.weekday()))
            end_of_week = (start_of_week + timedelta(days=6))

            start_of_week_str = start_of_week.strftime('%Y-%m-%d')
            end_of_week_str = end_of_week.strftime('%Y-%m-%d')

            mcd_scanner = MCDExtractWindow(dir_mod08=str(self.networkdrive_letter) + self.MOD08_D3,
                                               dir_mde=str(self.networkdrive_letter) + self.MDE,
                                               dir_temp=str(self.networkdrive_letter) + self.TEMP_COEF,
                                               ini_date=start_of_week_str,
                                               end_date=end_of_week_str,
                                               bounding_shp=self.roi)

            dataset_info_mod08 = mcd_scanner.run_extraction_mod08d3()

            self.aod = dataset_info_mod08['AOD_mean'].mean()
            self.water_vapour = dataset_info_mod08['WV_mean'].mean()
            self.ozone = dataset_info_mod08['OZ_mean'].mean() / 1000  # in cm_atm

        # Montly mean values:
        if (np.isnan(self.aod) or self.aod == 0.0) or (np.isnan(self.water_vapour) or self.water_vapour == 0.0) or (np.isnan(self.ozone) or self.ozone == 0.0):
            month_name = calendar.month_name[self.datetime.month]

            gdf = self.roi

            self.aod = mcd_scanner.get_modis_monthly_mean(month_name, 'AOD', gdf)
            self.water_vapour = mcd_scanner.get_modis_monthly_mean(month_name, 'Water_Vapor', gdf)
            self.ozone = mcd_scanner.get_modis_monthly_mean(month_name, 'Total_Ozone', gdf) / 1000

class DateTime(object):

    """
    Stores date and time values.
    """

    day = float("nan")
    month = float("nan")
    year = float("nan")
    time_hh = float("nan")

    def __str__(self):
        return 'day: %f, month: %f, year: %f, timehh: %f' % (self.day, self.month, self.year, self.time_hh)