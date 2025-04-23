import os
import glob
import shutil
import rasterio
import numpy as np

'''
BRDF correction for aquatic environments using the approach proposed by Lee et al., 2011.

'''
def solve_2nd_order_poly(A, B, C):

    """ Solve 2nd order polynomial inversion
    where coefficients are xr dataArray
    Take only positive solution, otherwise provide 0.
    """

    # Compute solution according to sign of delta
    delta = B * B - 4 * A * C
    mask = delta > 0

    # By default (and when delta < 0), take value at extremum
    x = - B / (2 * A)

    # When delta > 0, take the largest solution
    x_1 = (-np.where(mask, B, np.nan) - np.sqrt(np.where(mask, delta, np.nan))) / (2 * np.where(mask, A, np.nan))
    x_2 = (-np.where(mask, B, np.nan) + np.sqrt(np.where(mask, delta, np.nan))) / (2 * np.where(mask, A, np.nan))

    x_sol = np.where(x_1 > x_2, x_1, x_2)  # Choose the larger root
    x = np.where(delta > 0, x_sol, x)  # Assign solution where delta > 0

    # Take only positive values
    x = np.where(x > 0, x, 0.0)

    return x

def get_closest_bbw_value(wavelength, bbw_water):
    closest_wavelength = min(bbw_water.keys(), key=lambda x: abs(x - wavelength))
    return bbw_water[closest_wavelength]

# G coefficients for different angles [solar zenith, sensor nadir, sensor azimuth] (ref: Lee et al., 2011)
G_0_w = {
    (0, 0, 0): 0.0604,
    (0, 30, 90): 0.0596,
    (15, 30, 90): 0.0590,
    (30, 30, 90): 0.0584,
    (0, 40, 135): 0.0581,
    (15, 40, 135): 0.0614,
    (30, 40, 135): 0.0624
}

G_1_w = {
    (0, 0, 0): 0.0406,
    (0, 30, 90): 0.0516,
    (15, 30, 90): 0.0562,
    (30, 30, 90): 0.0601,
    (0, 40, 135): 0.0581,
    (15, 40, 135): 0.0524,
    (30, 40, 135): 0.0524
}

G_0_p = {
    (0, 0, 0): 0.0402,
    (0, 30, 90): 0.0408,
    (15, 30, 90): 0.0411,
    (30, 30, 90): 0.0418,
    (0, 40, 135): 0.0414,
    (15, 40, 135): 0.0425,
    (30, 40, 135): 0.0434
}

G_1_p = {
    (0, 0, 0): 0.1310,
    (0, 30, 90): 0.1420,
    (15, 30, 90): 0.1461,
    (30, 30, 90): 0.1492,
    (0, 40, 135): 0.1458,
    (15, 40, 135): 0.1408,
    (30, 40, 135): 0.1406
}

def brdf_correction_lee(Rrs_array, sensor):

    '''
    BRDF correction for aquatic environments using the approach proposed by Lee et al., 2011.

    Args:
        Rrs_array (np.array): Array of remote sensing reflectance values. Shape: (rows, columns, bands)
        sensor (str): Satellite sensor name ("landsat" or "sentinel")
    Returns:
        np.array: Array of corrected remote sensing reflectance values
    '''

    # Ref: Zhang et al., (2009), Mason et al., (2016), Pope and Fry (1997), Pitarch and Vanhellemont (2021)
    if sensor == "landsat":
        aw = {443: 0.00626, 490: 0.01274, 560: 0.06236, 665: 0.371, 865: 2.07}
        bbw_water = {443: 0.0025, 490: 0.001522, 560: 0.000811, 665: 0.000425, 865: 0}

        coeff_X = [-1.15467, -1.17869, -0.24566, -0.06989]
        coeff_Q = [0.022365, 0.548575, 0.167207, 0, 0]
    else:
        aw = {443: 0.00626, 490: 0.01545, 560: 0.0619, 665: 0.429, 865: 2.07}
        bbw_water = {443: 0.0025, 490: 0.001407, 560: 0.000817, 665: 0.000399, 865: 0}

        coeff_X = [-1.09651, -1.33678, -0.35707, -0.08409]
        coeff_Q = [-0.02085, 0.540187, 0.226931, 0.010022, 0]

    sat_bands = [443, 490, 560, 665, 865]

    # Remote Sensing reflectance of subsurface (rrs) (sr^-1)
    rrs_array = np.zeros(Rrs_array.shape)
    for i in np.arange(0, Rrs_array.shape[2]):
        rrs_array[:,:,i] = Rrs_array[:,:,i] / (0.52 + 1.7 * Rrs_array[:,:,i])

    # Absorption at a reference wavelenght (a(λ0)) (m^-1). Here, we considered as reference wavelength 560 nm.
    X = np.log((2 * Rrs_array[:,:,1]) / (Rrs_array[:,:,2] + ((5 * Rrs_array[:,:,3] * Rrs_array[:,:,3]) / Rrs_array[:,:,1])))
    a_lambda_zero = aw[560] + 10 ** (coeff_X[0] + coeff_X[1] * X + coeff_X[2] * (X * X) + coeff_X[3] * (X * X * X))

    # Particle backscattering slope
    BR = (Rrs_array[:,:,2] / Rrs_array[:,:,3])
    q_poly = coeff_Q[0] + coeff_Q[1] * BR + coeff_Q[2] * (BR * BR) + coeff_Q[3] * (BR * BR * BR) + coeff_Q[4] * (BR * BR * BR * BR)
    n_ = 2.0 * (1 - 1.2 * np.exp(-0.9 * q_poly))

    # Backscattering coefficient at reference wavelength
    u_ = np.zeros(rrs_array.shape)
    for i in np.arange(0, Rrs_array.shape[2]):
        u_[:,:,i] = (-0.089 + np.sqrt((0.089**2) + 4 * 0.1245 * rrs_array[:,:,i])) / (2 * 0.1245)

    bbp_lambda_zero = ((u_[:,:,2] * a_lambda_zero) / (1 - u_[:,:,2])) - bbw_water[560]

    # Backscattering coefficient at different wavelenghts (bbp(λ)) (m^-1) - QAA v6
    bbp_ = np.zeros(rrs_array.shape)
    bb = np.zeros(rrs_array.shape) # total backscattering coefficient

    for i in np.arange(0, Rrs_array.shape[2]):
        bbp_[:,:,i] = bbp_lambda_zero * np.power((560/sat_bands[i]), n_)
        bb[:,:,i] = bbp_[:,:,i] + bbw_water[sat_bands[i]]

    # 6. Absorption coefficient at different wavelenght (a(λ)) (m^-1) - using equation 20 from Lee et al., 2011
    absorption = np.zeros(rrs_array.shape)

    for i in np.arange(0, Rrs_array.shape[2]):
        absorption[:, :, i] = (1 - u_[:,:,i]) * (bbw_water[sat_bands[i]] + bbp_[:,:,i]) / u_[:,:,i]

    # 7. Calculation of the normalizer Rrs values, equation 14 from Lee et al., 2011
    Rrs_norm = np.zeros(rrs_array.shape)

    for i in np.arange(0, Rrs_array.shape[2]):

        kk = absorption[:,:,i] + bb[:,:,i]
        Rrs_norm[:,:,i] = (G_0_w[(0, 0, 0)] + G_1_w[(0, 0, 0)] * (bbw_water[sat_bands[i]]/kk)) * (bbw_water[sat_bands[i]]/kk) + (G_0_p[(0, 0, 0)] + G_1_p[(0, 0, 0)] * (bbp_[:,:,i] /kk)) * (bbp_[:,:,i] / kk)

    return Rrs_norm

def call_brdf_correction(img_path, output_path, satellite):

    if satellite == "landsat":

        bands = ["B1", "B2", "B3", "B4", "B5"]
        swir_band = [f for f in glob.glob(fr'{img_path}\*.tif') if "B6" in f][0]
    else:
        bands = ["B01", "B02", "B03", "B04", "B8A"]
        swir_band = [f for f in glob.glob(fr'{img_path}\*.tif') if "B11" in f][0]

    images_list = [f for f in glob.glob(fr'{img_path}\*.tif') if any(band in f for band in bands)]
    images_list = sorted(images_list, key=lambda x: next((i for i, band in enumerate(bands) if band in x), float('inf')))

    img_list = []


    for img in images_list:

        with rasterio.open(img) as src:

            _array = src.read()/np.pi
            _array[_array < 0] = np.nan

        img_list.append(_array)

    Rrs_array = np.stack(img_list, axis=0)[:,0,:,:].transpose(1,2,0)

    corrected_Rrs = brdf_correction_lee(Rrs_array, satellite)

    # Convert Rrs back to reflectance values
    corrected_rho = corrected_Rrs * np.pi

    with rasterio.open(images_list[1]) as ref:
        profile = ref.profile

    profile.update(
        dtype=rasterio.float32,
        count=corrected_rho.shape[2],
        compress="lzw"
    )

    file_names = [f for f in os.listdir(img_path) if any(band in f for band in bands)]

    for i in range(corrected_rho.shape[2]):

        if satellite == "landsat":
            output_filename = os.path.join(output_path, file_names[i].replace("temp", "brdf_corrected"))
        else:
            output_filename = os.path.join(output_path, f"brdf_corrected_{file_names[i]}")

        band_profile = profile.copy()
        band_profile.update(count=1)

        with rasterio.open(output_filename, "w", **band_profile) as dst:

            dst.write(corrected_rho[:, :, i].astype(rasterio.float32), 1)

    # copy swir band as new name
    if satellite == "landsat":
        swir_band_output = os.path.join(output_path, os.path.basename(swir_band).replace("temp", "brdf_corrected"))
    else:
        swir_band_output = os.path.join(output_path, f"brdf_corrected_{os.path.basename(swir_band)}")

    shutil.copy(swir_band, swir_band_output)
