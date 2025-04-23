import os
import glob
from typing import List
from src.aquavis.utils import toolbox
from src.aquavis.utils import satwutils
from src.aquavis.utils import io
from src.aquavis.processors.data_class import AquaVisDataLoader
from src.aquavis.processors.brdf import brdf_lee11_QAA_RGB as brdf

def get_landsat_bands(scene_path: str, bands: List[str]) -> List[str]:
    """Get sorted list of Landsat band files for the specified bands."""
    band_files = [
        f for f in glob.glob(os.path.join(scene_path, '*_B*.tif'))
        if any(band in f for band in bands)
    ]
    # Sort bands according to the order in the input list
    return sorted(
        band_files,
        key=lambda x: next((i for i, band in enumerate(bands) if band in x), float('inf'))
    )

def process_band(landsat_band: str, imgtemp_dir: str, sen2_epsg_code: str) -> str:
    """Reproject a single Landsat band to Sentinel projection."""
    temp_img = os.path.join(imgtemp_dir, f"temp_{os.path.basename(landsat_band)}")
    reprojected_data = satwutils.reproject(landsat_band, temp_img, sen2_epsg_code)
    return reprojected_data

def apply_brdf_correction(imgtemp_dir: str, apply_correction: bool) -> List[str]:
    """Apply BRDF correction if enabled and return processed images."""
    if apply_correction:
        brdf.call_brdf_correction(imgtemp_dir, imgtemp_dir, 'landsat')
        return [f for f in glob.glob(fr'{imgtemp_dir}\*.tif') if "brdf_corrected" in f]
    return [f for f in glob.glob(fr'{imgtemp_dir}\*.tif') if "temp" in f]

def gen_tiles(landsat_scene: str, output_path: str) -> None:
    """Generate tiles for a given Landsat scene by reprojecting, resampling, and clipping it to Sentinel tile geometry."""

    print(f"Processing Landsat scene: {landsat_scene}")

    params = AquaVisDataLoader().load_aquavis_data()

    io.validate_file(params, landsat_scene)

    # Define bands to process and get their files
    landsat_bands = ['B2', 'B3', 'B4', 'B5', 'B6']
    landsat_scene_all_bands = get_landsat_bands(landsat_scene, landsat_bands)

    # Create temporary directory for processing
    imgtemp_dir = os.path.join(params.output_dir_tiling, 'temp', f"temp_{os.path.basename(landsat_scene)}")

    sen2_epsg_code = toolbox.tile_epsg(params.tile)
    sen_tile_target_shp = toolbox.get_tile_shp(params.tile, sen2_epsg_code)

    for i, landsat_band in enumerate(landsat_scene_all_bands):

        reprojected_data = process_band(landsat_band, imgtemp_dir, sen2_epsg_code)

        base_dir = os.path.join(output_path, os.path.basename(os.path.dirname(landsat_bands[i])))
        satwutils.create_dir(base_dir)
        out_band = os.path.join(base_dir, os.path.basename(landsat_scene_all_bands[i]))

        satwutils.cut_images_res(reprojected_data, sen_tile_target_shp, out_band, 30)
