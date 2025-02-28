from src.satwater.atmcor.atmcor_water.atmcor_water import Gceratmos
from src.satwater.water_mask.WaterMaskClass import WaterMaskClass
from src.satwater.atmcor.atmcor_water.cloud.SatClouds.satclouds import Satcloud

def run_gceratmos(path_main, path_dest, satellite, aero_type="Maritime", path_roi=None):

    """
        Applies GCER Atmos correction to Sentinel-2 images.
    """

    networkdrive_letter = 'Z'

    if satellite == 'sentinel':
        sat = 'MSI_S2'
    else:
        sat = 'OLI_L8/9'

    mode = None
    fmask_env = r'C:\Users\tml411\AppData\Local\anaconda3\envs\fmask_env\python'
    aero_type = "Maritime"
    gceratmos_r = Gceratmos(path_main, path_dest, networkdrive_letter, sat, aero_type, mode)
    gceratmos_r.run()
    print(f'GCER Atmos correction applied.')

    cloud_mask = Satcloud(path_main, sat, path_dest, fmask_env)
    cloud_mask.run()
    print(f'Cloud mask applied.')

    # Apply water mask
    # wm = WaterMaskClass()
    # dest_path_water = path_dest + '\\water_mask'
    # wm.clip_water(path_dest, dest_path_water)
    # print(f'Water mask and glint correction applied.')