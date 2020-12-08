"""
Post processing pipeline functions
"""
import xarray


def generate_ff_image(ref_signal, sensor):
    """
    Generate the fixed pattern noise

    Parameters
    ----------
    ref_signal : xarray.DataArray
        reference signal for flat fielding
    sensor : pyeosim.sensor.TdiCMOS
        fitted sensor instance

    Returns
    -------
    fpn_image : xarray.DataArray
        xarray of same shape as input
    """
    ref_signal = xarray.concat([ref_signal] * 12, dim='repeat')
    im = sensor.transform(ref_signal).astype(float)
    # calculate dark signal
    zero_ref = xarray.zeros_like(ref_signal)
    dark_level = sensor.transform(zero_ref).mean()
    im = im - dark_level
    return (im / im.mean(['x', 'y'])).median('repeat')
    # divide by mean to normalise
    # return (im.mean('repeat') - im.mean(['repeat', 'x', 'y'])).round()


def noise_corrected_signal(signal, image_sensor, dark_sensor, ff_frame):
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
    ff_frame : xarray.DataArray
        fixed pattern noise array for imaging region

    Returns
    -------
    noise_corrected_image : xarray.DataArray
        return a noise corrected image
    """
    # generate dark level and image instance
    dark_frame = dark_sensor.transform(xarray.zeros_like(signal))
    return _flat_field(image_sensor.transform(signal), dark_frame, ff_frame)


def noise_corrected_reflectance(signal, reference_signal, image_sensor,
                                dark_sensor, ff_frame):
    """
    Generates a sensor signal and a dark signal and subtracts one from the
    other

    Parameters
    ----------
    signal : xarray.DataArray
        reflectance signal
    reference_signal : xarray.DataArray
        100% reflectance reference signal
    image_sensor : pyeosim.sensor.TdiCMOS
        fitted sensor instance for imaging region
    dark_sensor : pyeosim.sensor.TdiCMOS
        fitted sensor instance for non-imaging region
    ff_frame : xarray.DataArray
        fixed pattern noise array for imaging region

    Returns
    -------
    noise_corrected_image : xarray.DataArray
        return a noise corrected image
    """
    # generate dark level and image instance for reference signal
    dark_frame = dark_sensor.transform(xarray.zeros_like(reference_signal))

    ref_image = _flat_field(image_sensor.transform(reference_signal),
                            dark_frame, ff_frame)

    # divide signal by reference to give reflectance
    return signal.astype(float)/ref_image.astype(float)


def _flat_field(signal, dark_frame, ff_frame):
    # get the mean dark level for the dark region
    dark_level = dark_frame.mean()
    return (signal - dark_level) / ff_frame
