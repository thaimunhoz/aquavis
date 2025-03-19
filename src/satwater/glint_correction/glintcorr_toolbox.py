import numpy as np

def xml_to_dict(element):

    if len(element) == 0:
        return element.text

    result = {}

    for child in element:

        child_dict = xml_to_dict(child)

        if child.tag in result:
            if isinstance(result[child.tag], list):
                result[child.tag].append(child_dict)
            else:
                result[child.tag] = [result[child.tag], child_dict]
        else:
            result[child.tag] = child_dict

    return result


def paramglint(ang: dict, solar_zn: float, view_zn: float, optical_depth_total: float, nw: float) -> dict:

    """
    Returns the Fresnel's reflectance and direct transmittance from atmosphere.
    """

    # Incidence angle:
    raa = abs(ang['solar_az'] - ang['view_az'])

    raa = np.where(raa > 180, 360 - raa, raa)

    cosTheta = np.sqrt((((np.cos(ang['solar_zn'] * (np.pi / 180)) * np.cos(ang['view_zn'] * (np.pi / 180))) + (
                np.sin(ang['solar_zn'] * (np.pi / 180)) * np.sin(ang['view_zn'] * (np.pi / 180)) * np.cos(
            raa * (np.pi / 180)))) + 1) / 2)  # cos2Theta = 2cos^2(Theta) - 1:

    Theta = np.arccos(cosTheta)

    # Transmittance angle:
    Theta_t = np.arcsin((1 / nw) * np.sin(Theta))

    # Fresnel reflectance coefficient:
    rfresnel = 0.5 * (((np.sin(Theta - Theta_t) / np.sin(Theta + Theta_t)) ** 2) + (
                (np.tan(Theta - Theta_t) / np.tan(Theta + Theta_t)) ** 2))

    # Direct Transmittance:
    Tdir = (np.exp(-float(optical_depth_total) / np.cos(float(solar_zn) * (np.pi / 180)))) * (
        np.exp(-float(optical_depth_total) / np.cos(float(view_zn) * (np.pi / 180))))

    return {'rFresnel': rfresnel, 'Tdir': Tdir}