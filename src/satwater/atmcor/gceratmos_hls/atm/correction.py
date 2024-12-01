# -*- mode: python -*

import warnings
import numpy as np


class Correction:

    def __init__(self,
                 parameters,
                 atmosphere,
                 arr,
                 index):

        self.parameters = parameters
        self.atmosphere = atmosphere
        self.arr = arr
        self.index = index
        self.arr_r = np.inf

    def run(self):
        """
        It runs the estimation of the surface reflectance for each satellite.
        """
        # Atmospheric Correction -> Adapted equation from Vermote et al. (2016) -> Atmosphere spherical albedo is removed:
        # Description:
        # Tg_OG: Total transmission of Other Gases;
        # Tg_Tg_O3: Total transmission of the Ozone;
        # Tg_H20: Total transmission of the Water Vapor;
        # T_upward: Total transmittance upward (Rayleigh + Aerosol);
        # T_downward: Total transmittance downward (Rayleigh + Aerosol);
        # r_atm: Atmosphere intrinsic reflectance;
        # arr_TOA: Top Of Atmosphere reflectance;
        # arr_sr: Surface reflectance.
        # Gaseous total transmission:
        tg_OG_co = float(self.atmosphere.values[self.index]['tg_OG_co'])
        tg_OG_c02 = float(self.atmosphere.values[self.index]['tg_OG_c02'])
        tg_OG_o2 = float(self.atmosphere.values[self.index]['tg_OG_o2'])
        tg_OG_no2 = float(self.atmosphere.values[self.index]['tg_OG_no2'])
        tg_OG_ch4 = float(self.atmosphere.values[self.index]['tg_OG_ch4'])
        # Total transmission of Other Gases -> Tg_OG:
        Tg_OG = float(tg_OG_co * tg_OG_c02 * tg_OG_o2 * tg_OG_no2 * tg_OG_ch4)
        # Total transmission of the Ozone -> Tg_O3:
        Tg_O3 = float(self.atmosphere.values[self.index]['Tg_O3'])
        # Total transmission of the Water Vapor -> Tg_H2O:
        Tg_H20 = float(self.atmosphere.values[self.index]['Tg_H20'])
        # Total transmittance upward (Rayleigh + Aerosol) -> T_upward:
        T_upward = float(self.atmosphere.values[self.index]['T_upward'])
        # Total transmittance downward (Rayleigh + Aerosol) -> T_downward:
        T_downward = float(self.atmosphere.values[self.index]['T_downward'])
        # Total transmission of the atmosphere -> T_atm:
        T_atm = float(T_upward * T_downward)
        # Atmosphere intrinsic reflectance -> p_atm:
        r_atm = float(self.atmosphere.values[self.index]['p_atm'])
        # From DN_TOA to REFLECTANCE_TOA:
        if self.parameters.satellite == 'OLCI_S3':
            arr_toa = (self.arr * np.pi) / (float(self.parameters.rescale[self.index]) * (np.cos((float(self.parameters.geometry[self.index]['solar_zn']) * np.pi) / 180)))
        elif self.parameters.satellite == 'MSI_S2':
            nodata = np.where(self.arr == 0, 0, 1) # 0 refers to NaN values.
            arr_toa = (self.arr + float(self.parameters.rescale[self.index]['offset'])) / float(self.parameters.rescale[self.index]['qvalue'])
            arr_toa = np.where((arr_toa * nodata) == 0, -9999, (arr_toa * nodata))
        elif self.parameters.satellite == 'OLI_L8/9':
            nodata = np.where(self.arr == 0, 0, 1) # 0 refers to NaN values.
            arr_toa = ((self.arr * float(self.parameters.rescale[self.index]['mult'])) + float(self.parameters.rescale[self.index]['add'])) / (np.cos((float(self.parameters.geometry[self.index]['solar_zn']) * np.pi) / 180))
            arr_toa = np.where(abs(arr_toa * nodata) == 0, 0, (arr_toa * nodata))
        elif self.parameters.satellite == 'OCI_PACE':
            arr_toa = self.arr
        else:
            warnings.warn("The sensor type was not identified in CORRECTION >> run", UserWarning)
        # Correction --SURFACE REFLECTANCE:
        self.arr_sr = ((arr_toa / (Tg_OG*Tg_O3*Tg_H20)) - r_atm) / (T_atm)
