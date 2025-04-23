import os
import glob
from typing import List
from src.aquavis.processors.resample import bandpass
from src.aquavis.processors.resample import resample_toolbox as toolbox

def apply_processing_pipeline(input_path: str, output_base: str, original_band_names: List[str]) -> None:
    """
    Apply resampling and bandpass processing to Sentinel-2 bands.

    Args:
        temp_dir: Temporary directory path
        output_base: Base output directory
        params: Processing parameters dictionary
        original_band_names: Original band names for output naming
    """
    processed_bands = glob.glob(fr'{input_path}\*.tif')

    # Sort bands to maintain original order
    processed_bands = sorted(
        processed_bands,
        key=lambda x: next((i for i, band in enumerate(original_band_names) if band in x), float('inf'))
    )

    # Apply resampling to 30 meters to each band
    for i, band_path in enumerate(processed_bands):

        rescaled_s2 = toolbox.rescale_s2_to_dataarray(band_path, 30)

        band_name = os.path.basename(band_path)
        output_dir = os.path.join(output_base, band_name)
        toolbox.create_dir(output_base)

        bandpass.apply_bandpass(rescaled_s2, original_band_names[i], output_dir)

def gen_resample(sentinel_scene, output_path):
    """Resamples Sentinel-2 bands to match Landsat spatial resolution and applies bandpass correction."""

    SENTINEL_BANDS = ['B02', 'B03', 'B04', 'B8A', 'B11']

    apply_processing_pipeline(
        sentinel_scene,
        output_path,
        SENTINEL_BANDS
    )