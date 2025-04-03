import logging
from typing import Optional
from src.satwater.atmcor.atmcor_water.gceratmos import GCERAtmos
from src.satwater.atmcor.atmcor_water.cloud.SatClouds.satclouds import Satcloud

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class AtmosphericCorrector:
    """
    Handles atmospheric correction and cloud masking for satellite imagery.

    Attributes:
        NETWORK_DRIVE (str): Default network drive letter
        FMASK_ENV_PATH (str): Path to Fmask environment
        DEFAULT_AERO_TYPE (str): Default aerosol type
    """

    NETWORK_DRIVE = 'Z'
    FMASK_ENV_PATH = r'C:\Users\tml411\AppData\Local\anaconda3\envs\fmask_env\python'
    DEFAULT_AERO_TYPE = "Maritime"

    @staticmethod
    def _get_sensor_name(satellite: str) -> str:

        sensor_map = {
            'sentinel': 'MSI_S2',
            'landsat': 'OLI_L8/9'
        }

        if satellite not in sensor_map:
            raise ValueError(f"Unsupported satellite type: {satellite}. Must be 'sentinel' or 'landsat'")

        return sensor_map[satellite]

    def run_correction(self,
            path_main: str,
            path_dest: str,
            satellite: str,
            msi_tile: str,
            aero_type: Optional[str] = None,
            mode: Optional[str] = None
    ) -> None:
        """
        Run atmospheric correction and cloud masking for a satellite image.

        Args:
            path_main (str): Path to input image
            path_dest (str): Path for output files
            satellite (str): Satellite type ('sentinel' or 'landsat')
            msi_tile (str): Tile identifier
            aero_type (str, optional): Aerosol type. Defaults to "Maritime"
            mode (str, optional): Processing mode. Defaults to None

        Raises:
            RuntimeError: If processing fails at any stage
        """
        try:
            # Set defaults
            aero_type = aero_type or self.DEFAULT_AERO_TYPE
            sensor = self._get_sensor_name(satellite)

            # Run atmospheric correction
            logger.info(f"Starting atmospheric correction for {path_main}")
            gceratmos = GCERAtmos(
                path_main,
                path_dest,
                self.NETWORK_DRIVE,
                sensor,
                aero_type,
                mode,
                msi_tile
            )
            gceratmos.run()
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

        except Exception as e:
            logger.error(f"Processing failed: {str(e)}")
            raise RuntimeError(f"Failed to process image {path_main}") from e