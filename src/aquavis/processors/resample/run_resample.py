import logging

from src.aquavis.processors.ProcessorABC import Processors
from src.aquavis.processors.resample import resample
from src.aquavis.processors.data_class import AquaVisDataLoader
from src.aquavis.utils import io

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class RunResample(Processors):

    def run(self, path_main: str, path_dest: str) -> None:
        """Run glint correction on the input image."""

        params = AquaVisDataLoader().load_aquavis_data()

        io.validate_file(params, path_main)

        if params.select_sat == 'sentinel':
            resample.gen_resample(path_main, path_dest)

        logger.info("Sentinel resampling completed successfully")