import logging
from src.aquavis.processors.ProcessorABC import Processors
from src.aquavis.processors.tiling import tiles
from src.aquavis.processors.data_class import AquaVisDataLoader

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class RunTiling(Processors):

    def run(self, path_main: str, path_dest: str) -> None:
        """Run tiling on the input image."""

        params = AquaVisDataLoader().load_aquavis_data()

        if params.select_sat == 'landsat':
            tiles.gen_tiles(path_main, path_dest)

        logger.info("Landsat tiling completed successfully")