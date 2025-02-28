# -*- mode: python -*-

import Py6S as py6s
import numpy as np
import pandas as pd
import warnings

class Atmosphere:

    def __init__(self,parameters):

        self.parameters = parameters
        self.values = {}

    def run(self):
        """
        Returns the atmospheric parameters by 6S-code.
        """
        import os

        # Instantiate the 6S model-------------------------------------------------------------------------------------:
        s = py6s.SixS()  # Main class. It has attributes that allow you to set parameters, run 6S and access the outputs.
        # {Attributes - Input}:
        # Atmospheric Profiles-----------------------------------------------------------------------------------------:
        s.atmos_profile = py6s.AtmosProfile.UserWaterAndOzone(self.parameters.water_vapour, self.parameters.ozone) # Data retrieved of the MODIS products (MCD19A2.006 and MOD08_D3 V6).
        # Defines AOD (Aerosol Optical Depth)--------------------------------------------------------------------------:
        s.aot550 = self.parameters.aod
        # Aerosol Profiles---------------------------------------------------------------------------------------------:
        if self.parameters.aero_type == 'Continental':
            s.aero_profile = py6s.AeroProfile.PredefinedType(py6s.AeroProfile.Continental)
        else:
            s.aero_profile = py6s.AeroProfile.PredefinedType(py6s.AeroProfile.Maritime)  # Default
        # Altitudes---------------------------------------------------------------------------------------------------:
        s.altitudes = py6s.Altitudes()
        s.altitudes.set_sensor_satellite_level()  # Set the sensor altitude to be satellite level.
        s.altitudes.set_target_custom_altitude(self.parameters.altitude)  # The altitude of the target (in km).
        # Satellite selection:
        if self.parameters.satellite == 'MSI_S2':
            if self.parameters.type == 'S2A':
                rsr = pd.read_csv(r'C:\Users\tml411\Documents\Python Scripts\hls_water\src\satwater\atmcor\gceratmos_hls\atm\rsr_S2A_MSI.txt', sep='\t', skiprows=5, names=['wavelength', 'rsr', 'band', 'id'])
                bands = self.rsr_interp(rsr, 13, [0.412, 2.321])
            else:
                rsr = pd.read_csv(r'C:\Users\tml411\Documents\Python Scripts\hls_water\src\satwater\atmcor\gceratmos_hls\atm\rsr_S2B_MSI.txt', sep='\t', skiprows=5, names=['wavelength', 'rsr', 'band', 'id'])
                bands = self.rsr_interp(rsr, 13, [0.412, 2.321])
        elif self.parameters.satellite == 'OLI_L8/9':
            if self.parameters.type == 'LC08':
                rsr = pd.read_csv(r'C:\Users\tml411\Documents\Python Scripts\hls_water\src\satwater\atmcor\gceratmos_hls\atm\rsr_L8_OLI.txt', sep='\t', skiprows=5, names=['wavelength', 'rsr', 'band', 'id'])
                bands = self.rsr_interp(rsr, 8, [0.427, 2.356])
            else:
                rsr = pd.read_csv(r'C:\Users\tml411\Documents\Python Scripts\hls_water\src\satwater\atmcor\gceratmos_hls\atm\rsr_L9_OLI2.txt', sep='\t', skiprows=5, names=['wavelength', 'rsr', 'band', 'id'])
                bands = self.rsr_interp(rsr, 8, [0.427, 2.356])
        else:
            warnings.warn("The sensor type was not identified in ATMOSPHERE.", UserWarning)
        # Simulation:
        x = 1
        for i in range(x, len(bands) + x):
            # Geometries of view and illumination:
            s.geometry = py6s.Geometry.User()
            s.geometry.day = self.parameters.datetime.day
            s.geometry.month = self.parameters.datetime.month
            s.geometry.solar_a = self.parameters.geometry[i - x]['solar_az']
            s.geometry.solar_z = self.parameters.geometry[i - x]['solar_zn']
            s.geometry.view_a = self.parameters.geometry[i - x]['view_az']
            s.geometry.view_z = self.parameters.geometry[i - x]['view_zn']
            s.wavelength = py6s.Wavelength(np.min(bands[i][1]), np.max(bands[i][1]), bands[i][0])
            s.run()
            # Output:
            self.values[i - x] = {'view_zn': self.parameters.geometry[i - x]['view_zn'],
                                  'view_az': self.parameters.geometry[i - x]['view_az'],
                                  'solar_zn': self.parameters.geometry[i - x]['solar_zn'],
                                  'solar_az': self.parameters.geometry[i - x]['solar_az'],
                                  'Tg_O3': float(s.outputs.transmittance_ozone.total),
                                  'p_atm': float(s.outputs.atmospheric_intrinsic_reflectance),
                                  'optical_depth_rayleigh': float(s.outputs.optical_depth_total.rayleigh)}

    def rsr_interp(self, rsr, band_number: int, range_w: list) -> dict:
        """
        It interpolates the RSR in 2.5 nm.
        """
        output = {}
        for i in range(1, band_number + 1):
            filter_ = rsr.loc[rsr['id'] == i]
            w_min_ = np.min(filter_['wavelength'])
            w_max_ = np.max(filter_['wavelength'])
            # Range required --it interpolates at 2.5 nm:
            wavelength_ = [np.around(w, 3) for w in np.arange(w_min_, w_max_, 0.001)]
            wavelength_int = [np.around(w, 3) for w in np.arange(range_w[0], range_w[1], 0.001)]
            wavelength_sub = list(set(wavelength_int) - set(wavelength_))
            df_ = pd.DataFrame({'wavelength': wavelength_sub, 'rsr': int(0), 'band': str(i)})
            df_ = pd.concat([filter_, df_]).sort_values(by=['wavelength'], ascending=True)
            d_select = df_['rsr'].to_list()
            wavelength_interp = [np.around(w, 4) for w in np.arange(range_w[0], range_w[1], 0.0025)]
            norm = np.interp(wavelength_interp, df_['wavelength'].to_list(), d_select).tolist()
            out_ = [norm, [i for i in wavelength_interp]]
            output[i] = out_
        return output

