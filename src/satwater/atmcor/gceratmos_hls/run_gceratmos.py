from gceratmos import Gceratmos

def run_gceratmos(path_main, path_dest, satellite, path_roi=None):

    """
        Applies GCER Atmos correction to Sentinel-2 images.
    """

    networkdrive_letter = 'Z'

    if satellite == 'sentinel':
        sat = 'MSI_S2'
    else:
        sat = 'OLI_L8/9'

    mode = None

    gceratmos_r = Gceratmos(path_main, path_dest, networkdrive_letter, sat, mode)
    gceratmos_r.run()

    print(f'GCER Atmos correction applied.')