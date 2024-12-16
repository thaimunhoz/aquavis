import os
import logging
from typing import Optional, Dict, Any

from src.satwater.atmcor import atm6s as atmcor
from src.satwater.tiling import tiles as tiling
from src.satwater.tiling import resample as resample
from .hlswater import generate_hls as hlswater
from .visualization import plot_images as visualization
from .glint_correction import glint_corr as glintcorr

class SatWater(object):

    """
    A class to manage and run the full SatWater processing chain.

    This includes:
        - Atmospheric correction
        - Image tiling (for Landsat images)
        - Image resampling (for Sentinel images)
        - Glint correction
        - Synthetic HLS image generation
        - Visualization and true color composition

    Args:
        select_sat (Optional[str]): Satellite name ("landsat" or "sentinel").
        params (Optional[Dict[str, Any]]): Dictionary of parameters for processing.

    Attributes:
        select_sat (str): Selected satellite.
        params (Dict[str, Any]): Parameters for processing.
    """

    def __init__(self, select_sat: Optional[str] = None, params: Optional[Dict[str, Any]] = None):
        self.select_sat = select_sat
        self.params = params

    def _create_output_dir(self) -> None:
        """
        Create the output directory if it does not already exist.
        """
        output_dir = self.params.get("output_dir")
        if not output_dir:
            raise ValueError("Output directory is not specified in the parameters.")

        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            logging.info(f"Output directory created: {output_dir}")
        else:
            logging.info(f"Output directory already exists: {output_dir}")

    def run_atmcor(self) -> None:
        """
        Perform atmospheric correction using the 6S model.
        """
        try:
            self._create_output_dir()
            atmcor.run(self.select_sat, self.params)
            logging.info("Atmospheric correction completed successfully.")
        except Exception as e:
            logging.error(f"Error during atmospheric correction: {e}")
            raise

    def run_tiling(self) -> None:
        """
        Perform tiling of Landsat images, reprojecting and clipping them to Sentinel MGRS tiles.
        """
        if self.select_sat != "landsat":
            logging.warning("Tiling is only applicable to Landsat images.")
            return

        try:
            tiling.run(self.select_sat, self.params)
            logging.info("Tiling completed successfully.")
        except Exception as e:
            logging.error(f"Error during tiling: {e}")
            raise

    def run_resample(self) -> None:
        """
        Perform resampling of Sentinel images to match Landsat spatial resolution (30m).
        """
        if self.select_sat != "sentinel":
            logging.warning("Resampling is only applicable to Sentinel images.")
            return

        try:
            resample.run(self.params)
            logging.info("Resampling completed successfully.")
        except Exception as e:
            logging.error(f"Error during resampling: {e}")
            raise

    def run_glint(self) -> None:
        """
        Apply glint correction to the images.
        """
        try:
            glintcorr.run(self.params)
            logging.info("Glint correction completed successfully.")
        except Exception as e:
            logging.error(f"Error during glint correction: {e}")
            raise

    def run_hlswater(self) -> None:
        """
        Generate synthetic HLS water images.
        """
        try:
            hlswater.run(self.params)
            logging.info("HLS water image generation completed successfully.")
        except Exception as e:
            logging.error(f"Error during HLS water generation: {e}")
            raise

    def run_plot(self) -> None:
        """
        Generate and save true color compositions of the processed images.
        """
        try:
            visualization.run(self.params)
            logging.info("Visualization and true color composition completed successfully.")
        except Exception as e:
            logging.error(f"Error during visualization: {e}")
            raise