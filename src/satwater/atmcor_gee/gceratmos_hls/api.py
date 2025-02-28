# -*- mode: python -*-

from gceratmos_gee import Gceratmos
from cloud.SatClouds.satclouds import Satcloud

"""
GCERatmos (version 1) is a proto-algorithm designed for atmospheric correction in coastal and inland waters based on MODIS data and 6SV radiative transfer code.
In this version you will find a processor for Sentinel-2A/B | Landsat-8/9 | Sentinel-3A/B | PACE/OCI.
<GeneralInfo>
    <version>2<version>
    <date>2024-07-19<date>
    <team>Rejane_Paulino<team>
"""

# {Input}
path_main = r"Z:\dbcenter\images\sentinel\scenes\level_toa\39TXM\20230427\S2B_MSIL1C_20230427T071619_N0509_R006_T39TXM_20230427T075636.SAFE" # str
path_roi = r"C:\Users\tml411\Documents\papers\hls_synthetic\database\1\S2\S2B_MSIL1C_20230427T071619_N0509_R006_T39TXM_20230427T075636.SAFE\*.shp" # str + .shp
path_dest = r"C:\Users\tml411\Documents\papers\hls_synthetic\database\1\S2\S2B_MSIL1C_20230427T071619_N0509_R006_T39TXM_20230427T075636.SAFE\toa_corrected" # str
path_buffer = r"C:.../*.shp" # str + .shp # requested to PACE/OCI
networkdrive_letter = 'Z' # str #networkdrive associated to GCER data server
sat = 'MSI_S2'  # str # e.g., ['OLCI_S3', 'MSI_S2', 'OLI_L8/9', 'OCI_PACE']
mode = 'gee' # e.g., ['gee' or None]

# {Usage}
# Sentinel-2A/B:
gceratmos_r = Gceratmos(path_main, path_roi, path_dest, networkdrive_letter, sat, mode)
gceratmos_r.run()

cloud_mask = Satcloud(path_main, sat, path_dest + '/' + path_main[-65:])
cloud_mask.run()

# Landsat-8/9:
#gceratmos_r = Gceratmos(path_main, path_roi, path_dest, networkdrive_letter, sat, mode)
#gceratmos_r.run()

#cloud_mask = Satcloud(path_main, sat, path_dest)
#cloud_mask.run()

# Sentinel-3A/B:
#gceratmos_r = Gceratmos(path_main, path_roi, path_dest, networkdrive_letter, sat, mode)
#gceratmos_r.run()

#cloud_mask = Satcloud(path_main, sat, path_dest + '/' + path_main[-99:])
#cloud_mask.run()

# OCI-PACE:
"""Note: The correction is based on Level-1B (TOA radiance).
            */ The path main refers to path r"C:.../PACE_OCI_*_LIC_V2.nc;
            */ The path_buffer is requested here;
            */ The output array are saved in a NetCDF structure."""

#gceratmos_r = Gceratmos(path_main, path_roi, path_dest, networkdrive_letter, sat, mode, path_buffer)
#gceratmos_r.run()


#########################################
from gceratmos_gee import Gceratmos
from cloud.SatClouds.satclouds import Satcloud

import pandas as pd

datax = pd.read_csv(r"R:\guser\rsp\PhD\papers\hls_synthetic\dataset\images\S2X_MSIL1C__\sentinel2_catalog.csv").reset_index()

networkdrive_letter = 'R' # str
sat = 'MSI_S2'  # e.g., ['OLCI_S3', 'MSI_S2', 'OLI_L8/9']
mode = 'gee'

z = []
for i in datax.index:
    try:
        data = datax.loc[datax.index == i]
        # # Atmospheric Correction:
        gceratmos_r = Gceratmos(data['paths_S2'].values[0], data['ROIbbox'].values[0], data['dest'].values[0], networkdrive_letter, sat, mode)
        gceratmos_r.run()
        # Cloud masks:
        cloud_mask = Satcloud(data['paths_S2'].values[0], sat, data['dest'].values[0] + '/' + data['paths_S2'].values[0][-65:])
        cloud_mask.run()
    except:
        z.append(data)
        pass
pd.concat(z).to_csv(r'R:\guser\rsp\PhD\papers\hls_synthetic\dataset\images\S2X_MSIL1C__\error_s2.csv')

import os
import glob
from gceratmos_gee import Gceratmos
import pandas as pd

# {Input}
path_main = r'R:\guser\rsp\PhD\papers\pace_sentinel3\AERONET-OCLWN15\csv_paths' # str
path_buffer = r'x__x'
destx = r'R:\guser\rsp\PhD\papers\pace_sentinel3\dataset\SR\AERONET-OCLWN15\PACE_OCI_L1B_'
networkdrive_letter = 'R' # str #networkdrive associated to GCER data server
sat = 'OCI_PACE'  # str # e.g., ['OLCI_S3', 'MSI_S2', 'OLI_L8/9', 'OCI_PACE']
mode = None # e.g., ['gee' or None]

for i in os.listdir(path_main):
    data = pd.read_csv(path_main + '/' + i)
    os.makedirs(destx + '/' + i[:-4], exist_ok=True)
    for index, row in data.iterrows():
        aod = float(row['AOD'])
        wv = float(row['WV(cm)'])
        oz = float(row['OZ'])
        alt = float(row['ALT'])
        img_path = str(row['PACE_OCI_L1B_'])
        path_roi = str(row['roi'])
        path_dest = str(row['dest'])
        gceratmos_r = Gceratmos(img_path, path_roi, path_dest, networkdrive_letter, sat, aod, wv, oz, alt, mode, path_buffer)
        gceratmos_r.run()







#
# paths = [i for i in glob.glob(os.path.join(path_main , '*.nc'))]
#
# for i in paths:
#
#     gceratmos_r = Gceratmos(i, path_roi, path_dest, networkdrive_letter, sat, mode, path_buffer)
#     gceratmos_r.run()

