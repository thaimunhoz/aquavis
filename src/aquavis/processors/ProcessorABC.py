from abc import ABC, abstractmethod
from typing import Dict

class Processors(ABC):
    """
    This class provides a template for processing satellite imagery, including
    atmospheric correction, glint correction, tiling, resampling, water mask generation,
    product generation, and visualization.

    Args:
        params (Dict[str, Any]): Dictionary of parameters for processing.
    """
    @abstractmethod
    def run(self, path_main: str, path_dest: str) -> None:
        """ Run the step algorithm """
        pass