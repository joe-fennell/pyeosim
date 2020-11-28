"""
Detector simulation stages. These funcs should be chained together by
writing new subclasses of GenericTransformer in 'sensor' submodule.
"""

import numpy
import xarray
from ._decorators import return_equal_xarray


def voltage_to_DN(voltage, v_ref, adc_gain, bit_depth):
    """
    Model a linear response ADC.

    Parameters
    ----------
    voltage : xarray.DataArray
        voltage
    v_ref : float
        reference voltage of ADC
    adc_gain : float
        ADC gain
    bit_depth : int
        ADC depth

    Returns
    -------
    DN : xarray.DataArray
        simulated digital number values
    """
    max_DN = numpy.int(2**bit_depth - 1)
    DN = (adc_gain * (v_ref - voltage)).round()
    DN = DN.where(DN <= max_DN, DN)
    return DN.where(DN > 0, 0)


def electron_to_voltage(electron_count, v_ref, sense_node_gain,
                        full_well):
    """
    Converts electrons (charge) to voltage

    Parameters
    ----------
    electron_count : xarray.DataArray
        charge or electron count
    v_ref : float
        reference voltage
    sense_node_gain : float
        gain in V/e-
    full_well : float
        max numberof electrons per pixel

    Returns
    -------
    voltage : xarray.DataArray
        voltage of sensor
    """
    # truncate at full well
    e = electron_count.round().where(electron_count <= full_well, full_well)
    return v_ref - e * sense_node_gain


def add_poisson_noise(electron_count, Lambda):
    """
    additive poisson noise

    Parameters
    ----------
    electron_count : xarray.DataArray
        electron count instance
    Lambda : float
        Lambda coefficient of poisson distribution

    Returns
    -------
    electron_count : xarray.DataArray
        electron count with added dark noise
    """
    gaus_noise = numpy.random.poisson(Lambda,
                                      size=electron_count.shape)
    return (electron_count + gaus_noise).round()


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
    Convert photon count to electron count using sensor quantum efficiency.
    Taking average over all wavelengths

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


def photon_mean(flux, pixel_area, integration_time):
    """
    Convert sensor irradiance (wavelength in nm, irradiance in W m-2)
    to photon count (poisson lambda coefficient).

    Parameters
    ----------
    flux : xarray.DataArray
        pixel flux (photons/second/m2)
    pixel_area : float
        area in microns
    integration_time : float
        exposure time in seconds

    Returns
    -------
    photo_count : xarray.DataArray
        per-pixel photon count
    """
    return flux * integration_time * _microns2_to_m2(pixel_area)


def irradiance_to_flux(irradiance):
    """
    Convert xarray spectral object (wavelength in nm, irradiance in W m-2)
    to photon flux using hc=1.99e-25

    Parameters
    ----------
    ar : xarray.DataArray
        radiance
    """
    hc = 1.99e-25  # J m
    lambda_ = _nm_to_m(irradiance.wavelength)
    return (lambda_/hc) * irradiance


def radiance_to_irradiance(radiance, altitude):
    """
    Convert radiance (W sr-1). Only an approximation for nadir sensor

    Parameters
    ----------
    radiance : xarray.DataArray
        at sensor radiance
    altitude : float
        sensor altitude (m)

    Returns
    -------
    irradiance : xarray.DataArray
        sensor irradiance
    """
    # ground sample distance in metres
    gsd = float(radiance.x[1] - radiance.x[0])
    # sensor altitude m
    IFOV = numpy.arctan2(gsd, altitude)
    return radiance * IFOV


def band_Qe(SRFs, quantum_efficiency):
    """
    Calculates weighted mean of the quantum efficiency in
    each spectral channel.

    Parameters
    ----------
    SRFs : dict
        Spectral Response Function dictionary
    quantum_efficiency : xarray.DataArray
        Q_e with wavelength coordinate

    Returns
    -------
    Q_e : xarray.DataArray
        Weighted mean quantum efficiency for each band
    """
    bands = numpy.arange(len(SRFs))
    QEs = []
    for srf in SRFs.values():
        weight = srf/srf.integrate('wavelength')
        QEs.append(float((quantum_efficiency * weight).sum()))
    return xarray.DataArray(QEs, [('band', bands)])


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
