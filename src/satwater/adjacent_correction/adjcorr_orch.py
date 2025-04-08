import os
import glob
import shutil
import rasterio
import numpy as np
from typing import Dict, List, Any
import xml.etree.ElementTree as ET
from scipy.signal import fftconvolve

from src.satwater.utils import satwutils
from src.satwater.adjacent_correction import adj_corr as adjcorr


class AdjCorrOrchestrator:

    """Handles adjacent correction for satellite imagery."""

    def __init__(self):

        self.adjcorr = adjcorr.AdjCorrClass()

    def run(self, params: Dict[str, Any]) -> None:

        """
        Run adjacent correction on all scenes in the input directory.

        Args:
            params: Dictionary containing processing parameters
        """

        # Setup input directory
        base_path = os.path.join(params['output_dir'], 'atmcor')
        # Get input scenes based on satellite type
        input_paths = self.adjcorr._get_input_paths(params, base_path)

        # Setup output directory
        params['output_dir_adjcorr'] = os.path.join(params['output_dir'], 'adjcorr')
        satwutils.create_dir(params['output_dir_adjcorr'])

        # Process each scene
        for scene_path in input_paths:

            scene_name = self.adjcorr._get_scene_name(params, scene_path)

            output_path = os.path.join(params['output_dir_adjcorr'], os.path.basename(scene_name))

            if os.path.exists(output_path):
                print(f"Skipping {output_path}, already exists.")
                continue

            satwutils.create_dir(output_path)

            self.adjcorr.apply_adjcorr(scene_name, output_path, params)