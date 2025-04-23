import os
from pathlib import Path
from typing import List, Optional, Union
import logging

logger = logging.getLogger(__name__)

def ensure_dir(directory: Union[str, Path]) -> None:
    """
    Ensure a directory exists, creating it if necessary.
    
    Args:
        directory: Directory path to ensure exists
    """
    os.makedirs(directory, exist_ok=True)
    logger.debug(f"Ensured directory exists: {directory}")

def list_files(directory: Union[str, Path], pattern: str = "*") -> List[Path]:
    """
    List files in a directory matching a pattern.
    
    Args:
        directory: Directory to search
        pattern: Glob pattern to match files
        
    Returns:
        List of matching file paths
    """
    directory = Path(directory)
    if not directory.exists():
        raise FileNotFoundError(f"Directory not found: {directory}")
        
    files = list(directory.glob(pattern))
    logger.debug(f"Found {len(files)} files matching pattern '{pattern}' in {directory}")
    return files

def validate_file(metadata, scene_name) -> None:

    # check if all bands exists in the path_main folder
    sentinel_bands = ["B02", "B03", "B04", "B8A", "B11"]
    landsat_bands = ["B2", "B3", "B4", "B5", "B6"]

    if metadata.select_sat == "landsat":
        bands = landsat_bands
    else:
        bands = sentinel_bands

    for band in bands:
        if not any(band in file for file in os.listdir(scene_name)):
            raise FileNotFoundError(f"Band {band} not found in the provided path.")

def get_output_path(base_dir: Union[str, Path], subdir: str) -> Path:
    """
    Get a standardized output path for a processing step.
    
    Args:
        base_dir: Base output directory
        subdir: Subdirectory for the processing step
        
    Returns:
        Path to the output directory
    """
    output_dir = Path(base_dir) / subdir
    ensure_dir(output_dir)
    return output_dir

def move_file(src: Union[str, Path], dst: Union[str, Path]) -> None:
    """
    Move a file from source to destination.
    
    Args:
        src: Source file path
        dst: Destination file path
    """
    src = Path(src)
    dst = Path(dst)
    
    if not src.exists():
        raise FileNotFoundError(f"Source file not found: {src}")
        
    ensure_dir(dst.parent)
    src.rename(dst)
    logger.debug(f"Moved file from {src} to {dst}")

def copy_file(src: Union[str, Path], dst: Union[str, Path]) -> None:
    """
    Copy a file from source to destination.
    
    Args:
        src: Source file path
        dst: Destination file path
    """
    src = Path(src)
    dst = Path(dst)
    
    if not src.exists():
        raise FileNotFoundError(f"Source file not found: {src}")
        
    ensure_dir(dst.parent)
    import shutil
    shutil.copy2(src, dst)
    logger.debug(f"Copied file from {src} to {dst}") 