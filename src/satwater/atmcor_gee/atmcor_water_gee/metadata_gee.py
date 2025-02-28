import os
import glob
import calendar
import numpy as np
import pandas as pd
import rioxarray as rio
import geopandas as gpd
from shapely.geometry import box
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
                 xda_vnir,
                 metadata_xda,
                 mode=None):

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
        self.xda_vnir = xda_vnir
        self.metadata_xda = metadata_xda
        self.mode = mode

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

        # path = [i for i in glob.glob(os.path.join(self.path_main + self.GRANULE, '*')) if 'L1C_' in i]
        #
        # self.s2path = path[0] + self.IMG_DATA
        #
        # # Return the bouding box of the image:
        # self.roi = tool.return_bbox(os.path.join(self.s2path, os.listdir(self.s2path)[3]))

        # Return the bouding box of the image:
        self.xda_vnir = self.xda_vnir.rename({"X": "x", "Y": "y"})
        self.xda_vnir.rio.set_spatial_dims("x", "y")
        bounds = self.xda_vnir.rio.bounds()
        bbox_polygon = box(*bounds)
        self.roi = gpd.GeoDataFrame({'geometry': [bbox_polygon]}, crs=self.xda_vnir.rio.crs)

        self.bandname = ["B2", "B3", "B4", "B8A", "B11"]
        self.type = self.path_main.split('_')[0] # safe number A or B

        self.date_and_time()
        self.geo()
        self.read_coefficient()
        #self.rescale_factor()

        df = pd.DataFrame({'img': [self.path_main], 'aod': [self.aod], 'wv': [self.water_vapour], 'oz': [self.ozone], 'alt': [self.altitude]})

        if not os.path.exists(self.path_dest): os.makedirs(self.path_dest)
        df.to_csv(self.path_dest + '/' + 'atm_parameters.csv')

    def date_and_time(self):

        """
        Returns the date and the time.
        """

        date = self.xda_vnir["time"].values[0].astype('datetime64[ms]').item()

        # Export date and time in a metadata structure:
        value = DateTime()
        value.day = date.day
        value.month = date.month
        value.year = date.year
        value.time_hh = date.hour + (date.minute / 60) + (date.second / 3600)

        self.datetime = value

    def geo(self):

        """
        Returns the geometry of observation and illumination.
        """

        # Sun azimuth [az] and zenith [zn] angle:
        # View azimuth [az] and zenith [zn] angles:
        # It considers the angle average of the scene.
        count = 0
        bands = ["B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B8A", "B9", "B10", "B11", "B12"]

        for i in bands:

            output = {}
            output['solar_az'] = float(self.metadata_xda["features"][0]["properties"]["MEAN_SOLAR_AZIMUTH_ANGLE"]) # MEAN_SOLAR_AZIMUTH_ANGLE
            output['solar_zn'] = float(self.metadata_xda["features"][0]["properties"]["MEAN_SOLAR_ZENITH_ANGLE"]) #MEAN_SOLAR_ZENITH_ANGLE

            band_name_view_az = "MEAN_INCIDENCE_AZIMUTH_ANGLE_" + str(i)
            output['view_az'] = float(self.metadata_xda["features"][0]["properties"][band_name_view_az]) # MEAN_INCIDENCE_AZIMUTH_ANGLE

            band_name_view_zn = "MEAN_INCIDENCE_ZENITH_ANGLE_" + str(i)
            output['view_zn'] = float(self.metadata_xda["features"][0]["properties"][band_name_view_zn]) # MEAN_INCIDENCE_ZENITH_ANGLE

            self.geometry[count] = output
            count += 1

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
                 xda_vnir,
                 xda_angles,
                 metadata_xda,
                 mode=None):

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
        self.xda_vnir = xda_vnir
        self.xda_angles = xda_angles
        self.metadata_xda = metadata_xda

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

        self.bandname = ["B2", "B3", "B4", "B5", "B6"] # band name
        #path = [i for i in glob.glob(os.path.join(self.path_main, '*.xml')) if self.MTD_ID in i]

        # Return the bouding box of the image:
        self.xda_vnir = self.xda_vnir.rename({"X": "x", "Y": "y"})
        self.xda_vnir.rio.set_spatial_dims("x", "y")
        bounds = self.xda_vnir.rio.bounds()
        bbox_polygon = box(*bounds)
        self.roi = gpd.GeoDataFrame({'geometry': [bbox_polygon]}, crs=self.xda_vnir.rio.crs)

        #self.dict_metadata = tool.xml_to_json(str(path[0])) # metadata from sensor
        self.type = self.path_main.split('_')[0] # safe number L8 or L9
        self.date_and_time()
        self.geo()
        self.read_coefficient()
        self.rescale_factor()

        df = pd.DataFrame({'img': [self.path_main], 'aod': [self.aod], 'wv': [self.water_vapour], 'oz': [self.ozone], 'alt': [self.altitude]})

        if not os.path.exists(self.path_dest): os.makedirs(self.path_dest)

        df.to_csv(self.path_dest + '/' + 'atm_parameters.csv')

    def date_and_time(self):

        """
        Returns the date and the time.
        """

        date = self.xda_vnir["time"].values[0].astype('datetime64[ms]').item()

        # Export date and time in a metadata structure:
        value = DateTime()
        value.day = date.day
        value.month = date.month
        value.year = date.year
        value.time_hh = date.hour + (date.minute / 60) + (date.second / 3600)

        self.datetime = value

    def geo(self):

        """
        Returns the geometry of observation and illumination.
        """

        # Sun azimuth [az] and zenith [zn] angle:
        # View azimuth [az] and zenith [zn] angles:
        # It considers the angle averages of the scene.
        # The band 4 (red) is used as reference because it is near the center of the OLI/Landsat-8/9 focal plane.

        for i in range(0, 8):

            output = {}

            output['solar_az'] = np.nanmean(self.xda_angles["SAA"].values[0,:,:])
            output['solar_zn'] = np.nanmean(self.xda_angles["SZA"].values[0,:,:])
            output['view_az'] = np.nanmean(self.xda_angles["VAA"].values[0,:,:])
            output['view_zn'] = np.nanmean(self.xda_angles["VZA"].values[0,:,:])

            #print(f"ANGLE OUTPUT FOR BAND {i + 1}: {output}")

            self.geometry[i] = output

    def rescale_factor(self):

        """
        Returns the factor values to convert from DN_TOA to REFLECTANCE_TOA
        """

        # It obtains the factor values to convert from DN_TOA to REFLECTANCE_TOA:
        for i in range(1, 9):
            ADD_BAND = float(self.metadata_xda["features"][0]["properties"]['REFLECTANCE_ADD_BAND_' + str(i)])
            MULT_BAND = float(self.metadata_xda["features"][0]["properties"]['REFLECTANCE_MULT_BAND_' + str(i)])
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