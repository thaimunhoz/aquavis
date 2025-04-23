import logging
from typing import Dict
from src.aquavis.config.config import SatWaterConfig
from src.aquavis.processors.ProcessorABC import Processors
from src.aquavis.processors.atmcor.atm_aquavis import AquaVisAtmos
from src.aquavis.processors.atmcor.cloud.SatClouds.satclouds import Satcloud
from src.aquavis.processors.data_class import AquaVisDataLoader

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class AtmosphericCorrection(Processors):

    def __init__(self) -> None:
        self.params = AquaVisDataLoader().load_aquavis_data()

        config = SatWaterConfig()._load_paths()
        self.FMASK_ENV_PATH = config["fmask_env_path"]

    @staticmethod
    def _get_sensor_name(satellite: str) -> str:

        sensor_map = {
            'sentinel': 'MSI_S2',
            'landsat': 'OLI_L8/9'
        }

        return sensor_map[satellite]

    def run(self, path_main: str, path_dest: str) -> None:
        """Run atmospheric correction and cloud masking for a satellite image."""

        sensor = self._get_sensor_name(self.params.select_sat)

        # Run AquaVisAtmos modified
        aquavisatmos = AquaVisAtmos(
            path_main,
            path_dest,
            sensor,
        )
        aquavisatmos.run()
        logger.info("Atmospheric correction completed successfully")

        # Run cloud masking
        logger.info("Starting cloud masking")
        cloud_mask = Satcloud(
            path_main,
            sensor,
            path_dest,
            self.FMASK_ENV_PATH
        )
        cloud_mask.run()
        logger.info("Cloud masking completed successfully")