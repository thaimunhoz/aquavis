import os.path
import rasterio
import shutil

def band2band_sen2land():

    # landsat = (sentinel - intercept) / slope

    dict_bands={}

    dict_bands['blue']={}  # 486
    dict_bands['green']={}  # 560
    dict_bands['red']={}  # 665
    dict_bands['nir']={}  # 865

    dict_bands['blue']['multi'] = 1.0161
    dict_bands['blue']['add'] = -0.0005

    dict_bands['green']['multi'] = 0.9959
    dict_bands['green']['add'] = 0.0006

    dict_bands['red']['multi'] = 0.9440
    dict_bands['red']['add'] = -0.0002

    dict_bands['nir']['multi'] = 0.9487
    dict_bands['nir']['add'] = 0.0015

    return dict_bands

def find_wave2band(band_wave):

    if band_wave>=460 and band_wave<=500:
        band_nm = 'blue'
    elif band_wave>=540 and band_wave<=580:
        band_nm = 'green'
    elif band_wave>=640 and band_wave<=680:
        band_nm ='red'
    elif band_wave>=850 and band_wave<=880:
        band_nm = 'nir'  # NIR narrow 8A
    else:
        band_nm = ''
    return band_nm

def find_name2band(band_wave):

    if band_wave == 'B02':
        band_nm = 'blue'

    elif band_wave == 'B03':
        band_nm = 'green'

    elif band_wave == 'B04':
        band_nm ='red'

    elif band_wave == 'B8A':
        band_nm = 'nir'  # NIR narrow 8A

    else:

        band_nm = ''

    return band_nm

def apply_bandpass(input_img, outp_img):

    dict_bands = band2band_sen2land()

    band_name = os.path.basename(input_img).split('_')[2].split('.')[0]
    band_nm = find_name2band(band_name)

    if band_nm == '':  # Sentinel bands without OLI-like bands
        shutil.copy2(input_img, outp_img)

    else:

        src = rasterio.open(input_img)
        raster_arr = src.read()

        raster_arr = (raster_arr[0, :, :] - dict_bands[band_nm]['add'])/dict_bands[band_nm]['multi']

        dst_meta = src.meta.copy()
        dst_meta['compress'] = 'lzw'

        with rasterio.open(outp_img, 'w', **dst_meta) as dst:
            dst.write(raster_arr, 1)