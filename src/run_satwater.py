import logging
from src.satwater import build_satwater

def run_satwater(select_sat: str, tile: str, period_ini: str, period_end: str, output_dir: str,) -> None:

    '''
    Run the SatWater processing chain.
    Input:
        select_sat (str): Satellite to process (landsat or sentinel)
        tile (str): Tile name
        period_ini (str): Initial date in the format 'YYYYMMDD'
        period_end (str): End date in the format 'YYYYMMDD'
        output_dir (str): Output directory
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
        },
        "sentinel": {
            "input_dir": r"Z:\dbcenter\images\sentinel\scenes\level_toa",
            "tiles_shp": tiles_sentinel,
        },
        "landsat": {
            "input_dir": r"Z:\dbcenter\images\landsat\scenes\level_toa",
            "generation": "L89",
            "tiles_shp": tiles_landsat,
        },
        "output_dir": output_dir,
    }

    params[select_sat]["tiles"] = tile

    # Initialize SatWater object
    try:
        SatWater_i = build_satwater.SatWater(select_sat, params)
    except Exception as e:
        logging.error(f"Failed to initialize SatWater: {e}")
        raise

    # Run processing steps
    SatWater_i.run_atmcor() # 1. Atmospheric correction

    if select_sat == "landsat":
        SatWater_i.run_tiling() # 2. Tiling for Landsat

    elif select_sat == "sentinel":
        SatWater_i.run_resample() #3. Resampling and bandpass adjustment for Sentinel-2

    #SatWater_i.run_glint() # 4. Glint correction

    SatWater_i.run_hlswater() # 5. HLS water generation

    SatWater_i.run_plot() # 6. Plotting

    logging.info("Process finished successfully.")