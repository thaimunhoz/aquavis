import os
import shutil
from typing import Optional
from src.aquavis.processors.data_class import AquaVisDataLoader
from src.aquavis.processors.atmcor import toolbox_gceratmos as tool
from src.aquavis.processors.atmcor.metadata.metadata_msi import Metadata_MSI_S2
from src.aquavis.processors.atmcor.metadata.metadata_oli import Metadata_OLI_L89
from src.aquavis.processors.atmcor.atm.coefficient_run import RunAtmCoefficients

class AquaVisAtmos:

    """Handles atmospheric correction for satellite imagery."""

    def __init__(self, path_main: str, path_dest: str, satellite: str):
        """Initialize the atmospheric corrector."""
        self.path_main = path_main
        self.path_dest = path_dest
        self.satellite = satellite.upper()

        self.params = AquaVisDataLoader().load_aquavis_data()

    def run(self) -> None:
        """Execute the atmospheric correction process. Select the appropriate method based on the satellite."""
        processor = {
            'MSI_S2': self._process_sentinel,
            'OLI_L8/9': self._process_landsat
        }.get(self.satellite)

        processor()

    def _process_sentinel(self) -> None:
        """Process Sentinel-2 MSI data."""
        print(self.path_main[-65:])

        output_dir = tool.newdirectory(self.path_dest, self.path_main[-65:])
        temp_dir = tool.newdirectory(output_dir, 'tempdir')
        tool._convert_jp2_to_tiff(self.path_main, temp_dir)

        Metadata_MSI_S2(self.satellite, self.path_main).run()

        RunAtmCoefficients(self.path_main, output_dir)._read_atmospheric_coefficients()

        tool._perform_atmospheric_correction(self.path_main, output_dir, temp_dir)
        shutil.rmtree(temp_dir)

    def _process_landsat(self) -> None:
        """Process Landsat 8/9 OLI data."""
        print(self.path_main[-40:])

        output_dir = tool.newdirectory(self.path_dest, self.path_main[-40:])

        Metadata_OLI_L89(self.satellite, self.path_main).run()

        RunAtmCoefficients(self.path_main, output_dir)._read_atmospheric_coefficients()

        tool._perform_atmospheric_correction(self.path_main, output_dir)