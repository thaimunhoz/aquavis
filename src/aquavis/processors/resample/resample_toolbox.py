import os
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
import xarray as xr
import numpy as np
import rioxarray as rio

def rescale_s2_to_dataarray(image_path: str, target_resolution: float) -> xr.DataArray:
    """
    Rescale a Sentinel-2 image to the desired spatial resolution and return it as an xarray.DataArray.

    Args:
        image_path (str): Path to the input Sentinel-2 band image.
        target_resolution (float): Target spatial resolution in meters.

    Returns:
        xarray.DataArray: Rescaled image as an in-memory DataArray.
    """

    da = rio.open_rasterio(image_path)
    resampled_da = da.rio.reproject(
        da.rio.crs,
        resolution=target_resolution,
        resampling=Resampling.bilinear
    )

    # with rasterio.open(image_path) as src:
    #     transform, width, height = calculate_default_transform(
    #         src.crs, src.crs, src.width, src.height, *src.bounds,
    #         dst_width=int(src.width * src.res[0] / target_resolution),
    #         dst_height=int(src.height * src.res[1] / target_resolution)
    #     )
    #
    #     rescaled = np.full((height, width), -9999, dtype=src.dtypes[0])
    #
    #     reproject(
    #         source=rasterio.band(src, 1),
    #         destination=rescaled,
    #         src_transform=src.transform,
    #         src_crs=src.crs,
    #         dst_transform=transform,
    #         dst_crs=src.crs,
    #         nodata=-9999,
    #         resampling=Resampling.bilinear
    #     )
    #
    #     # Wrap into xarray.DataArray
    #     rescaled_da = xr.DataArray(
    #         rescaled,
    #         dims=("y", "x"),
    #         coords={
    #             "y": np.arange(height) * -target_resolution + transform.f,
    #             "x": np.arange(width) * target_resolution + transform.c
    #         },
    #         attrs={
    #             "crs": src.crs.to_string(),
    #             "transform": transform.to_gdal(),
    #             "resolution": (target_resolution, target_resolution),
    #             "nodata": -9999
    #         }
    #     )

    return resampled_da

def create_dir(path):

    if not os.path.exists(path):

        os.makedirs(path)