import logging
from src.satwater import build_aquavis

def run_satwater(select_sat: str, tile: str, period_ini: str, period_end: str, output_dir: str, keep_atmcor: str, input_toa: str, brdf_corr: str, output_type: str) -> None:

    '''
    Run the AQUAVis processing chain.

    Input:
        select_sat (str): Satellite to process (landsat or sentinel)
        tile (str): Tile name (e.g., '35VLG')
        period_ini (str): Initial date in the format 'YYYYMMDD'
        period_end (str): End date in the format 'YYYYMMDD'
        output_dir (str): Output directory
        keep_atmcor (bool): Keep atmospheric correction files
        input_toa (str): Input TOA source ('GEE' or 'local')
        brdf_corr (bool): Apply BRDF correction
        output_type (str): Output type ('rho' or 'rrs')

    Output:
        None
    '''

    tiles_sentinel = r'C:\Users\tml411\Documents\Python Scripts\hls_water\src\satwater\auxfiles\tiles\MGRS_tiles.shp'
    tiles_landsat = r'C:\Users\tml411\Documents\Python Scripts\hls_water\src\satwater\auxfiles\tiles\ms_landsat.shp'

    # Validate inputs
    if select_sat not in {"landsat", "sentinel"}:
        raise ValueError(f"Invalid satellite selected: {select_sat}. Choose 'landsat' or 'sentinel'.")

    if len(period_ini) != 8 or len(period_end) != 8:
        raise ValueError("Dates must be in the format 'YYYYMMDD'.")

    n_cores=12

    # Parameters setup
    params = {
        "aux_info": {
            "period": (period_ini, period_end),
            "n_cores": n_cores,
            "sat_name": select_sat,
            "keep_atmor": keep_atmcor,
            "brdf_corr": brdf_corr,
            "output_type": output_type,
        },
        "sentinel": {
            "input_dir": r"Z:\guser\tml\mypapers\HLS_package_paper\sentinel_toa",
            "tiles_shp": tiles_sentinel,
        },
        "landsat": {
            "input_dir": r"Z:\guser\tml\mypapers\HLS_package_paper\landsat_toa",
            "generation": "L89",
            "tiles_shp": tiles_landsat,
        },
        "output_dir": output_dir,
    }

    params[select_sat]["tiles"] = tile

    # Initialize SatWater object
    try:
        SatWater_i = build_aquavis.SatWater(select_sat, params)
    except Exception as e:
        logging.error(f"Failed to initialize SatWater: {e}")
        raise

    # Run processing steps
    if input_toa == "GEE":
        SatWater_i.run_atmcor_gee()  # 1. Atmospheric correction using GEE as input TOA source

    else:
        SatWater_i.run_atmcor() # 1. Atmospheric correction using local TOA source

    SatWater_i.run_adjcorr() # 2. Adjacent correction

    SatWater_i.run_glint_corr() # 3. Glint correction

    if select_sat == "landsat":
        SatWater_i.run_tiling() # 4. Tiling for Landsat

    elif select_sat == "sentinel":
        SatWater_i.run_resample() # 4. Resampling and bandpass adjustment for Sentinel-2

    SatWater_i.run_water_mask() # 5. Water mask generation

    SatWater_i.run_hlswater() # 6. HLS water generation

    #SatWater_i.run_plot() # 7. Plotting

    logging.info("Process finished successfully.")