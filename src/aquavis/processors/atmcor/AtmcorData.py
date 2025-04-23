from dataclasses import dataclass
from typing import Dict, Any, List, Optional
import numpy as np
from datetime import datetime

@dataclass
class MetadataSatellite:
    path_main: str
    path_dest: str
    networkdrive_letter: str
    satellite: str
    aero_type: str
    geometry: Dict[int, Dict[str, float]]
    rescale: Dict[int, Dict[str, float]]
    dict_metadata: Dict[str, Any]
    msi_tile: Optional[str] = None
    type: str = "nan"
    bandname: Optional[Dict] = None
    aod: float = np.nan
    water_vapour: float = np.nan
    ozone: float = np.nan
    altitude: float = np.nan
    datetime: datetime = datetime.now()
    roi: Any = None
    values: Dict[int, Dict[str, float]] = None
    values_adjcorr: Dict[int, Dict[str, float]] = None