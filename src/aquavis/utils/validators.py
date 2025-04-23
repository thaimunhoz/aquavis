from typing import Any, Dict, List, Optional, Union
import re
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

def validate_satellite(satellite: str) -> None:
    """
    Validate satellite name.
    
    Args:
        satellite: Satellite name to validate
        
    Raises:
        ValueError: If satellite name is invalid
    """
    valid_satellites = ['landsat', 'sentinel']
    if satellite not in valid_satellites:
        raise ValueError(f"Invalid satellite: {satellite}. Must be one of {valid_satellites}")

def validate_tile(tile: str, satellite: str) -> None:
    """
    Validate tile name format.
    
    Args:
        tile: Tile name to validate
        satellite: Satellite name
        
    Raises:
        ValueError: If tile name format is invalid
    """
    if satellite == 'landsat':
        # Landsat tile format: path/row (e.g., '123/045')
        if not re.match(r'^\d{3}/\d{3}$', tile):
            raise ValueError(f"Invalid Landsat tile format: {tile}. Expected format: 'path/row'")
    else:  # sentinel
        # Sentinel tile format: MGRS (e.g., '35VLG')
        if not re.match(r'^\d{2}[A-Z]{3}$', tile):
            raise ValueError(f"Invalid Sentinel tile format: {tile}. Expected format: 'MGRS'")

def validate_date(date: str) -> None:
    """
    Validate date format (YYYYMMDD).
    
    Args:
        date: Date string to validate
        
    Raises:
        ValueError: If date format is invalid
    """
    if not re.match(r'^\d{8}$', date):
        raise ValueError(f"Invalid date format: {date}. Expected format: YYYYMMDD")
        
    # Basic date validation
    year = int(date[:4])
    month = int(date[4:6])
    day = int(date[6:8])
    
    if not (1 <= month <= 12):
        raise ValueError(f"Invalid month in date: {date}")
    if not (1 <= day <= 31):  # This is a simple check, could be more precise
        raise ValueError(f"Invalid day in date: {date}")

def validate_config(config: Dict[str, Any], required_keys: List[str]) -> None:
    """
    Validate configuration dictionary.
    
    Args:
        config: Configuration dictionary to validate
        required_keys: List of required configuration keys
        
    Raises:
        ValueError: If required keys are missing
    """
    missing_keys = [key for key in required_keys if key not in config]
    if missing_keys:
        raise ValueError(f"Missing required configuration keys: {missing_keys}")

def validate_path(path: Union[str, Path], must_exist: bool = True) -> Path:
    """
    Validate a file or directory path.
    
    Args:
        path: Path to validate
        must_exist: Whether the path must exist
        
    Returns:
        Path object
        
    Raises:
        FileNotFoundError: If path must exist but doesn't
        ValueError: If path is invalid
    """
    path = Path(path)
    
    if must_exist and not path.exists():
        raise FileNotFoundError(f"Path does not exist: {path}")
        
    return path

def validate_output_type(output_type: str) -> None:
    """
    Validate output type.
    
    Args:
        output_type: Output type to validate
        
    Raises:
        ValueError: If output type is invalid
    """
    valid_types = ['rrs', 'rho']
    if output_type not in valid_types:
        raise ValueError(f"Invalid output type: {output_type}. Must be one of {valid_types}") 