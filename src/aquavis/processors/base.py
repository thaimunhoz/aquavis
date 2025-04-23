import os
from dataclasses import dataclass
from typing import Dict, Any, List, Tuple, Optional, Callable
from src.aquavis.utils import toolbox_processors as toolbox

from src.aquavis.processors.data_class import AquaVisDataLoader
from src.aquavis.processors.atmcor.run_atmcor_water import AtmosphericCorrection
from src.aquavis.processors.adjacent_correction.adj_corr import AdjCorrClass
from src.aquavis.processors.glint_correction.fresglint import FresGLINT
from src.aquavis.processors.tiling.run_tiling import RunTiling
from src.aquavis.processors.resample.run_resample import RunResample
from src.aquavis.processors.water_mask.water_mask import RunWaterMask
from src.aquavis.processors.aquavis_product_generator.generate_aquavis import AQUAVisProductGenerator

@dataclass
class ProcessorConfig:
    processor: Any
    input_func: Callable
    output_func: Callable

class ProcessorFunctions():

    ''' Pre-processing funtions'''

    def __init__(self):
        self.params = AquaVisDataLoader().load_aquavis_data()

    def _create_processors(self):
        """Create processors for each step."""

        return [
            ProcessorConfig(
                processor=AtmosphericCorrection(),
                input_func=lambda: toolbox.atmcor_input_files(self.params.select_sat, self.params.input_dir, self.params.tile, self.params.period),
                output_func=lambda: toolbox.atmcor_output_files(self.params.select_sat, self.params.input_dir, self.params.tile, self.params.period, self.params.output_dir_atmcor),
            ),
            ProcessorConfig(
                processor=AdjCorrClass(),
                input_func=lambda: toolbox.input_files(self.params.output_dir_atmcor),
                output_func=lambda: toolbox.output_files(self.params.output_dir_adjcorr, toolbox.input_files(self.params.output_dir_atmcor)),
            ),
            ProcessorConfig(
                processor=FresGLINT(),
                input_func=lambda: toolbox.input_files(self.params.output_dir_adjcorr),
                output_func=lambda: toolbox.output_files(self.params.output_dir_glintcorr, toolbox.input_files(self.params.output_dir_adjcorr)),
            ),
            ProcessorConfig(
                processor=RunTiling(),
                input_func=lambda: toolbox.input_files(self.params.output_dir_glintcorr),
                output_func=lambda: toolbox.output_files(self.params.output_dir_tiling, toolbox.input_files(self.params.output_dir_glintcorr)),
            ),
            ProcessorConfig(
                processor=RunResample(),
                input_func=lambda: toolbox.input_files(self.params.output_dir_glintcorr),
                output_func=lambda: toolbox.output_files(self.params.output_dir_tiling,toolbox.input_files(self.params.output_dir_glintcorr)),
            ),
            ProcessorConfig(
                processor=RunWaterMask(),
                input_func=lambda: toolbox.input_files(self.params.output_dir_tiling),
                output_func=lambda: toolbox.output_files(self.params.output_dir_water_mask,toolbox.input_files(self.params.output_dir_tiling)),
            ),
            ProcessorConfig(
                processor=AQUAVisProductGenerator(),
                input_func=lambda: toolbox.input_files(self.params.output_dir_water_mask),
                output_func=lambda: toolbox.output_files(self.params.output_dir_aquavis,toolbox.input_files(self.params.output_dir_water_mask)),
            )
        ]

    def run(self) -> None:

        processor_ = self._create_processors()

        for config in processor_:

            input_paths = config.input_func()
            output_paths = config.output_func()

            toolbox._process_images(config.processor, input_paths, output_paths)