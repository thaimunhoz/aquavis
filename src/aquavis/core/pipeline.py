import logging
from src.aquavis.config.config import SatWaterConfig
from src.aquavis.processors import base
from src.aquavis.processors.data_class import AquaVisDataLoader

def run_aquavis(select_sat: str, tile: str, period_ini: str, period_end: str, output_dir: str, output_type: str) -> None:

    '''
    Run the AQUAVis processing chain.

    Input:
        select_sat (str): Satellite to process (landsat or sentinel)
        tile (str): Tile name (e.g., '35VLG')
        period_ini (str): Initial date in the format 'YYYYMMDD'
        period_end (str): End date in the format 'YYYYMMDD'
        output_dir (str): Output directory
        output_type (str): Output type ('rho' or 'rrs')
    '''

    # Initialize configuration
    config = SatWaterConfig()

    params = config.get_params(
        select_sat=select_sat,
        tile=tile,
        period_ini=period_ini,
        period_end=period_end,
        output_dir=output_dir,
        output_type=output_type
    )

    loader = AquaVisDataLoader()
    shared_aquavis_data = loader.load(params)
    loader.save_aquavis_data(shared_aquavis_data)

    # Initialize AQUAVis object
    try:
        AQUAVis_i = base.ProcessorFunctions()
    except Exception as e:
        logging.error(f"Failed to initialize SatWater: {e}")
        raise

    AQUAVis_i.run()

    logging.info("Process finished successfully.")