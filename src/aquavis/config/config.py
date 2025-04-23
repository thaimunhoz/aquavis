import os
import json
from typing import Dict, Any

class SatWaterConfig:

    def __init__(self, config_path: str = "config_paths.json"):
        """
        Initialize configuration with paths from JSON file.

        Args:
            config_path: Path to JSON configuration file
        """
        self.config_path = os.path.join(os.path.dirname(__file__), "config_paths.json")
        self.paths = self._load_paths()

    def _load_paths(self) -> Dict[str, Any]:
        """Load paths from JSON configuration file."""
        try:
            with open(self.config_path, 'r') as f:
                return json.load(f)
        except FileNotFoundError:
            raise FileNotFoundError(f"Configuration file not found at {self.config_path}")
        except json.JSONDecodeError:
            raise ValueError(f"Invalid JSON in configuration file {self.config_path}")

    def get_params(self,
                   select_sat: str,
                   tile: str,
                   period_ini: str,
                   period_end: str,
                   output_dir: str,
                   output_type: str,
                   n_cores: int = 12) -> Dict[str, Any]:
        """
        Generate parameters dictionary based on configuration.

        Args:
            All the processing parameters from the original function
            n_cores: Number of cores to use (default: 12)

        Returns:
            Dictionary of parameters for SatWater processing
        """

        # Validate inputs
        if select_sat not in {"landsat", "sentinel"}:
            raise ValueError(f"Invalid satellite selected: {select_sat}. Choose 'landsat' or 'sentinel'.")

        if len(period_ini) != 8 or len(period_end) != 8:
            raise ValueError("Dates must be in the format 'YYYYMMDD'.")

        # Base parameters structure
        params = {
            "aux_info": {
                "period": (period_ini, period_end),
                "n_cores": n_cores,
                "sat_name": select_sat,
                "output_type": output_type,
            },
            "sentinel": {
                "input_dir": self.paths["sentinel_input_dir"],
                "tiles_shp": self.paths["tiles_sentinel"],
            },
            "landsat": {
                "input_dir": self.paths["landsat_input_dir"],
                "generation": "L89",
                "tiles_shp": self.paths["tiles_landsat"],
            },
            "output_dir": output_dir,
        }

        params[select_sat]["tile"] = tile

        return params