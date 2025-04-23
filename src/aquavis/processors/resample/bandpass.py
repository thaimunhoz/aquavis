import os.path
import rasterio
import shutil
import xarray as xr

def band2band_sen2land():

    # landsat = (sentinel - intercept) / slope

    dict_bands={}

    dict_bands['blue']={}  # 486
    dict_bands['green']={}  # 560
    dict_bands['red']={}  # 665
    dict_bands['nir']={}  # 865
    dict_bands['swir1'] = {}
    dict_bands['swir2'] = {}

    dict_bands['blue']['multi'] = 1.0161
    dict_bands['blue']['add'] = -0.0005

    dict_bands['green']['multi'] = 0.9959
    dict_bands['green']['add'] = 0.0006

    dict_bands['red']['multi'] = 0.9440
    dict_bands['red']['add'] = -0.0002

    dict_bands['nir']['multi'] = 0.9487
    dict_bands['nir']['add'] = 0.0015

    dict_bands['swir1']['multi'] = 1
    dict_bands['swir1']['add'] = 0

    dict_bands['swir2']['multi'] = 1
    dict_bands['swir2']['add'] = 0

    return dict_bands

def find_name2band(band_wave):

    if band_wave == 'B02':
        band_nm = 'blue'

    elif band_wave == 'B03':
        band_nm = 'green'

    elif band_wave == 'B04':
        band_nm ='red'

    elif band_wave == 'B8A':
        band_nm = 'nir'

    elif band_wave == 'B11':
        band_nm = 'swir1'

    elif band_wave == 'B12':
        band_nm = 'swir2'

    else:

        band_nm = ''

    return band_nm

def apply_bandpass(rescaled_da: xr.DataArray, band_name: str, outp_img: str):
    """
    Apply bandpass correction to a rescaled Sentinel-2 image stored as a DataArray.

    Args:
        rescaled_da (xr.DataArray): Rescaled input image as a DataArray.
        band_name (str): Name of the band (e.g., 'B03') to identify OLI-like band correction.

    Returns:
        xr.DataArray: Corrected image after bandpass adjustment.
    """
    dict_bands = band2band_sen2land()
    band_nm = find_name2band(band_name)

    # Skip bands that don't have OLI equivalents
    if band_nm == '':
        return rescaled_da.copy()

    # Apply the correction
    corrected = (rescaled_da - dict_bands[band_nm]['add']) / dict_bands[band_nm]['multi']

    # Return a new DataArray with updated attributes
    corrected.attrs.update(rescaled_da.attrs)
    corrected.attrs["bandpass_corrected"] = True
    corrected.attrs["band"] = band_name

    corrected.rio.to_raster(outp_img, driver='GTiff', dtype='float32')
