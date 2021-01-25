"""
Detector simulation stages. These funcs should be chained together by
writing new subclasses of GenericTransformer in 'sensor' submodule.
"""

import numpy
# import xarray
from ._decorators import return_equal_xarray


def add_column_offset(voltage, offset):
    """
    Adds the column offset noise to voltage array

    Parameters
    ----------
    voltage : xarray.DataArray
        voltage array
    offset : xarray.DataArray
        column voltage offset

    Returns
    -------
    voltage : xarray.DataArray
        new voltage
    """
    return voltage + offset


def add_dark_signal(electrons, dark_current, integration_time, dsnu):
    """
    Adds dark current to electron count

    Parameters
    ----------
    electrons : xarray.DataArray
        electron array
    dark_current : float
        dark current in e-/s/pixel
    integration_time : float
        integration period in seconds
    dsnu : xarray.DataArray
        dark current fixed pattern noise distribution
    """
    # Generate dark electron quantity
    n_dsig = dark_current * integration_time
    # make an empty copy of the electrons array to add the shot noise
    n_dsig_shot = electrons.copy()
    # convert mean I_dsig mean electron vals to poisson random process
    n_dsig_shot.values = numpy.random.poisson(n_dsig, size=electrons.shape)
    # add DSNU noise
    n_dark = n_dsig_shot + n_dsig_shot * dsnu
    return (electrons + n_dark).round()


def add_gaussian_noise(electron_count, sigma):
    """
    Adds a zero-mean gaussian to electron count

    Parameters
    ----------
    electron_count : xarray.DataArray
        electron count instance
    sigma : float
        SD of noise (electrons)

    Returns
    -------
    electron_count : xarray.DataArray
        electron count with added dark noise
    """
    gaus_noise = numpy.random.normal(0, sigma,
                                     size=electron_count.shape)
    return (electron_count + gaus_noise).round()


def add_photon_noise(photon_count):
    """
    Converts photon count to a Poisson random process simulating photon shot
    noise.

    Parameters
    ----------
    photon_count : xarray.DataArray
        mean photon counts

    Returns
    -------
    photon_estimate : xarray.DataArray
        photon count instance
    """
    @return_equal_xarray
    def poisson_rv(ar):
        return numpy.random.poisson(ar)
    # convert photon means to a poisson random process
    return photon_count.map_blocks(poisson_rv)


def add_prnu(electrons, prnu):
    """
    Adds Photon Response Non-Uniformity to electrons

    Parameters
    ----------
    electrons : xarray.DataArray
        electrons array
    prnu : xarray.DataArray
        fixed pattern noise array
    """
    # add photo non-uniformity
    return (electrons + electrons * prnu).round()


def CONU(ones, offset_factor):
    """
    Generates column offset non uniformity array

    Parameters
    ----------
    ones : xarray.DataArray
        ones of same shape as image
    offset_factor : float
        column offset factor

    Returns
    -------
    CONU : xarray.DataArray
        Column offset non uniformity array with zero mean
    """
    fpn = ones.isel(band=0) * numpy.random.normal(0, offset_factor,
                                                  size=ones.isel(band=0).shape)
    # drop band and band_name coords if exist
    for var in ['band', 'band_name']:
        try:
            fpn = fpn.drop_vars(var)
        except ValueError:
            pass
    return fpn


def DSNU(ones, dark_current, integration_time, dark_factor):
    """
    Calculates the dark signal non-uniformity array.

    Note that dark_current should be per actual pixel, not per pixel subarray

    This assumes a log-Normal distribution suitable for short integration times
    of less than 100 seconds

    Parameters
    ----------
    ones : xarray.DataArray
        array of ones in same shape as output array
    dark_current : float
        sensor level mean dark current
        in e-/s/pixel
    integration_time : float
        integration time of exposure
    dark_factor : float
        dark current FPN factor range(0...1)

    Returns
    -------
    electron_count : xarray.DataArray
        updated electron count
    """
    sigma = integration_time * dark_current * dark_factor
    # multiply by ones array to convert numpy array to xarray
    fpn = ones * numpy.random.lognormal(0, sigma, size=ones.shape)
    # subtract off the mean to generate a zero-mean distribution
    return fpn - fpn.mean()


def electron_to_voltage(electron_count, v_ref, sense_node_gain,
                        full_well, read_noise):
    """
    Converts electrons (charge) to voltage by simulating a sense node with a
    fixed v_ref

    Parameters
    ----------
    electron_count : xarray.DataArray
        charge or electron count
    v_ref : float
        reference voltage
    sense_node_gain : float
        gain in V/e-
    full_well : int
        max numberof electrons per pixel
    read_noise : int
        read noise in electrons per channel

    Returns
    -------
    voltage : xarray.DataArray
        voltage of sensor
    """
    # truncate at full well
    e = electron_count.where(electron_count < full_well, full_well)
    # add read noise
    e = add_gaussian_noise(e, read_noise).round()
    return v_ref - (e * sense_node_gain)
    # return e * sense_node_gain


def energy_to_quantity(energy):
    """
    Convert energy properties to quanta

    Parameters
    ----------
    ar : xarray.DataArray
        energy in J
    """
    hc = 1.98644586e-25  # J m
    lambda_ = _nm_to_m(energy.wavelength)
    return (lambda_/hc) * energy


def photon_mean(flux, pixel_area, integration_time):
    """
    Convert sensor photon flux density to photon count at sensor

    Parameters
    ----------
    flux : xarray.DataArray
        pixel flux density (photons/second/m2)
    pixel_area : float
        area in microns 2
    integration_time : float
        exposure time in seconds

    Returns
    -------
    photo_count : xarray.DataArray
        per-pixel photon count
    """
    return flux * integration_time * _microns2_to_m2(pixel_area)


def photon_to_electron(photon_count, Q_E):
    """
    Convert photon count to electron count using sensor quantum efficiency.
    Taking average over all wavelengths

    Parameters
    ----------
    photon_count : xarray.DataArray
        photon count instance
    Q_E : float, array-like
        float in range 0-1 specifying photon conversion rate or xarray with
        same wavelength sampling as photon_count
    """
    return (photon_count * Q_E).round()


def PRNU(ones, prnu_factor):
    """
    Calculate the Photo Response Non-Uniformity

    Parameters
    ----------
    ones : xarray.DataArray
        array of ones in same shape as output array
    prnu_factor : float
        Manufacturer quoted PRNU in range 0...1

    Returns
    -------
    prnu : xarray.DataArray
        array of same shape as ones
    """
    return ones * numpy.random.normal(0, prnu_factor, size=ones.shape)


def radiance_to_irradiance_2(radiance, lens_diameter, focal_length):
    """
    Convert radiance (W m-2 sr-1) to irradiance of sensor based on lens.
    Only an approximation for nadir sensor.
    Forshortening is not taken into account.
    Does not take into account cos^4 roll-off

    Parameters
    ----------
    radiance : xarray.DataArray
        at sensor radiance
    lens_diameter : float
        diameter of optic in metres
    focal_length : float
        focal length of optic in metres

    Returns
    -------
    irradiance : xarray.DataArray
        sensor irradiance
    """
    return radiance * (numpy.pi / 4) * ((lens_diameter / focal_length) ** 2)


# def radiance_to_irradiance(radiance, altitude):
#     """
#     Convert radiance (W sr-1). Only an approximation for nadir sensor.
#     Forshortening is not taken into account.
#
#     Parameters
#     ----------
#     radiance : xarray.DataArray
#         at sensor radiance
#     altitude : float
#         sensor altitude (m)
#
#     Returns
#     -------
#     irradiance : xarray.DataArray
#         sensor irradiance
#     """
#     # pixel ground area in m2
#     A_ground = float(radiance.x[1] - radiance.x[0])**2
#     # sensor altitude m
#     # This approach treats the pixel value as a point source with a measured
#     # radiant intensity (integrate out area: I = L A_ground
#     I_ = A_ground * radiance
#     # use inverse square rule to convert to irradiance:
#     # H = (I \cos \theta)/r**2
#     # This version assumes satellite is nadir to pixel (cos 0 = 1)
#     return I_ / (altitude**2)


def voltage_to_DN(voltage, v_ref, bit_depth):
    """
    Model a linear response ADC.

    Parameters
    ----------
    voltage : xarray.DataArray
        voltage
    v_ref : float
        reference voltage of ADC
    bit_depth : int
        ADC depth

    Returns
    -------
    DN : xarray.DataArray
        simulated digital number values
    """
    max_DN = numpy.int(2**bit_depth - 1)
    DN = max_DN * (v_ref - voltage)
    # DN = (voltage * (2**bit_depth))/v_ref
    # DN = (voltage / v_ref) * (2**bit_depth)
    DN = DN.round().where(DN < max_DN, max_DN)
    return DN.round().where(DN > 0, 0)


# Internal functions
def _nm_to_m(x):
    """
    Convert nanometer (nm) to metre (m)
    """
    return x * 1e-9


def _microns2_to_m2(x):
    """
    Convert nanometer
    """
    return x * (1e-6**2)
