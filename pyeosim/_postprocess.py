"""
Post processing pipeline functions
"""
import xarray


def generate_fp_image(zero_signal, sensor):
    """
    Generate the fixed pattern noise

    Parameters
    ----------
    zero_signal : xarray.DataArray
        array of zeroes of same shape as im
    sensor : pyeosim.sensor.TdiCMOS
        fitted sensor instance

    Returns
    -------
    fpn_image : xarray.DataArray
        xarray of same shape as input
    """
    im = sensor.transform(xarray.concat([zero_signal] * 12, dim='repeat'))
    return im.mean('repeat') - im.mean(['repeat', 'x', 'y'])


def subtract_dark_noise(signal, image_sensor, dark_sensor,
                        fpn_image, fpn_dark):
    """
    Generates a sensor signal and a dark signal and subtracts one from the
    other

    Parameters
    ----------
    signal : xarray.DataArray
        reflectance signal
    image_sensor : pyeosim.sensor.TdiCMOS
        fitted sensor instance for imaging region
    dark_sensor : pyeosim.sensor.TdiCMOS
        fitted sensor instance for non-imaging region
    fpn_image : xarray.DataArray
        fixed pattern noise array for imaging region
    fpn_dark : xarray.DataArray
        fixed pattern noise array for non-imaging region

    Returns
    -------
    noise_corrected_image : xarray.DataArray
        return a noise corrected image
    """
    # generate dark level and image instance
    dark_noise_level = dark_sensor.transform(xarray.ones_like(signal))
    image = image_sensor.transform(signal)
    # calculate mean dark noise level DN
    dark_noise_level = (dark_noise_level - fpn_dark).mean(['x', 'y'])
    return image - fpn_image - dark_noise_level
