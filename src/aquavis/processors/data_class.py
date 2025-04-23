import os
import pickle
import numpy as np
from dataclasses import dataclass
from typing import Dict, Any, List, Tuple, Optional
from src.aquavis.utils import toolbox as toolbox
from src.aquavis.config.config import SatWaterConfig

@dataclass
class AquaVisDataSet:
    select_sat: str
    input_dir: str
    tile: str
    period: Tuple[str, str]
    ncode: str
    output_type: str
    networkdrive_letter: str
    aerosol_type: str
    output_dir_atmcor: str
    output_dir_qaflag: str
    output_dir_adjcorr: str
    output_dir_glintcorr: str
    output_dir_tiling: str
    output_dir_water_mask: str
    output_dir_aquavis: str
    sensor: Optional[str] = None
    geometry: Optional[Dict[int, Dict[str, float]]] = None
    rescale: Optional[Dict[int, Dict[str, float]]] = None
    dict_metadata: Optional[Dict[str, Any]] = None
    type: Optional[str] = None
    bandname: Optional[Dict] = None
    aod: Optional[float] = np.nan
    water_vapour: float = np.nan
    ozone: float = np.nan
    altitude: float = np.nan
    datetime: Optional[Any] = None
    roi: Optional[Any] = None
    values: Optional[Dict[int, Dict[str, float]]] = None
    values_adjcorr: Optional[Dict[int, Dict[str, float]]] = None

class AquaVisDataLoader:


    def load(self, params: Dict = None) -> AquaVisDataSet:
        """Create and return an AquaVisDataSet from the parameters"""
        select_sat = params["aux_info"]["sat_name"]

        output_dirs = {
            'output_dir_atmcor': toolbox.create_dir(os.path.join(params["output_dir"], "atmcor")),
            'output_dir_qaflag': toolbox.create_dir(os.path.join(params["output_dir"], "QA_flags")),
            'output_dir_adjcorr': toolbox.create_dir(os.path.join(params["output_dir"], "adjcorr")),
            'output_dir_glintcorr': toolbox.create_dir(os.path.join(params['output_dir'], 'glintcorr')),
            'output_dir_tiling': toolbox.create_dir(os.path.join(params['output_dir'], 'tiling')),
            'output_dir_water_mask': toolbox.create_dir(os.path.join(params['output_dir'], 'water_mask')),
            'output_dir_aquavis': toolbox.create_dir(os.path.join(params['output_dir'], 'AQUAVis_Product'))
        }

        return AquaVisDataSet(
            select_sat=select_sat,
            output_type=params["aux_info"]["output_type"],
            input_dir=params[select_sat]["input_dir"],
            tile=params[select_sat]["tile"],
            period=params["aux_info"]["period"],
            networkdrive_letter="Z",
            aerosol_type="Maritime",
            ncode='L30' if select_sat == 'landsat' else 'S30',
            **output_dirs
        )

    def save_aquavis_data(self, aquavis_dataset: AquaVisDataSet) -> None:
        with open(os.path.join(os.path.dirname(__file__), "aquavis_dataset.pkl"), "wb") as f:
            pickle.dump(aquavis_dataset, f)

    def load_aquavis_data(self) -> AquaVisDataSet:
        with open(os.path.join(os.path.dirname(__file__), "aquavis_dataset.pkl"), "rb") as f:
            return pickle.load(f)