import os
import logging
from typing import Optional, Dict, Any

from src.satwater.atmcor import atm6s as atmcor
from src.satwater.tiling import tiles as tiling
from src.satwater.water_mask.WaterMaskClass import WaterMaskClass
from src.satwater.tiling import resample as resample
from .visualization import plot_images as visualization
from src.satwater.atmcor_gee import atm6s_gee as atmcor_gee
from src.satwater.glint_correction.fresglint import FresGLINT
from src.satwater.adjacent_correction.adj_corr import AdjCorrClass
from .aquavis_product_generator import generate_aquavis as aquavis

class SatWater(object):

    """
    A class to manage and run the full AQUAVis processing chain.

    This includes:
        - Atmospheric and glint correction
        - Image tiling (for Landsat images)
        - Image resampling (for Sentinel images)
        - Water mask generation
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

        self._create_output_dir()
        atmcor.run(self.select_sat, self.params)
        print("Atmospheric correction completed successfully.")

    def run_atmcor_gee(self) -> None:

        """
        Perform atmospheric correction using the 6S model.
        """

        self._create_output_dir()
        atmcor_gee.run(self.select_sat, self.params)
        print("Atmospheric correction completed successfully.")

    def run_adjcorr(self) -> None:

        """
        Perform adjacent correction using the 6S model.
        """

        adjcorr = AdjCorrClass()
        adjcorr.run(self.params)
        print("Adjacent correction completed successfully.")

    def run_glint_corr(self) -> None:

        """
        Perform glint correction using the 6S model.
        """

        glintcorr = FresGLINT()
        glintcorr.run(self.params)
        print("Glint correction completed successfully.")

    def run_tiling(self) -> None:

        """
        Perform tiling of Landsat images, reprojecting and clipping them to Sentinel MGRS tiles.
        """

        tiling.run(self.select_sat, self.params)
        print("Tiling completed successfully.")

    def run_resample(self) -> None:

        """
        Perform resampling of Sentinel images to match Landsat spatial resolution (30m).
        """

        resample.run(self.params)
        print("Resampling completed successfully.")

    def run_water_mask(self) -> None:

        """
        Generate water masks using image-based approach.
        """

        watermask = WaterMaskClass()
        watermask.run(self.params)

        print("Water mask generation completed successfully.")


    def run_hlswater(self) -> None:

        """
        Generate synthetic HLS water images.
        """

        aquavis.run(self.params)
        print("HLS water image generation completed successfully.")


    def run_plot(self) -> None:

        """
        Generate and save true color compositions of the processed images.
        """

        visualization.run(self.params)
        print("Visualization and true color composition completed successfully.")
