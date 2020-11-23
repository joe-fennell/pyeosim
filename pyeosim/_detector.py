"""
Detector simulation functions.
"""

import numpy
from .decorators import return_equal_xarray


def electron_to_DN(electron_count, sensitivity, bit_depth, baseline):
    """
    Model a linear response ADC.

    Parameters
    ----------
    electron_count : xarray.DataArray
        electron count instance
    sensitivity : float, array-like
        sensitivity (ADU/ e-)
    bit_depth : int
        ADC total bits
    baseline : int
        integer constant added to ADC values to prevent negative numbers

    Returns
    -------
    DN : xarray.DataArray
        simulated digital number values
    """
    max_adu = numpy.int(2**bit_depth - 1)
    adu = (electron_count * sensitivity).astype(numpy.int) + int(baseline)
    return adu.where(adu <= max_adu, max_adu)


def add_gaussian_noise(electron_count, sigma):
    """
    Model dark/readout noise as guassian additive noise.

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


def photon_to_electron(photon_count, Q_e):
    """
    Convert photon count to electron count using sensor quantum efficiency

    Parameters
    ----------
    photon_count : xarray.DataArray
        photon count instance
    Q_e : float, array-like
        float in range 0-1 specifying photon conversion rate or xarray with
        same wavelength sampling as photon_count
    """
    return (photon_count * Q_e).round()

@return_equal_xarray
def add_photon_noise(photon_count):
    """
    Sample a poisson distribution using input values as mean

    Parameters
    ----------
    photon_count : xarray.DataArray
        mean photon counts

    Returns
    -------
    photon_estimate : xarray.DataArray
        photon count instance
    """
    return photon_count.map_blocks(numpy.random.poisson)


def photon_mean(irradiance, pixel_area, integration_time):
    """
    Convert sensor irradiance (wavelength in nm, irradiance in W m-2)
    to photon count (poisson lambda coefficient).

    Parameters
    ----------
    irradiance : xarray.DataArray
        pixel irradiance with wavelength dim

    pixel_area : float
        area in microns

    integration_time : float
        exposure time in seconds

    Returns
    -------
    photo_count : xarray.DataArray
        per-pixel photon count
    """
    flux = _irradiance_to_flux(irradiance)
    return flux * integration_time * _microns2_to_m2(pixel_area)


# Internal functions
def _irradiance_to_flux(irradiance):
    """
    Convert xarray spectral object (wavelength in nm, irradiance in W m-2)
    to photon flux

    Parameters
    ----------
    ar : xarray.DataArray
        radiance
    """
    hc = 1.99e-25  # J m
    lambda_ = _nm_to_m(irradiance.wavelength)
    return (lambda_/hc) * irradiance


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
