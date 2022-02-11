"""Additional post processing tools for image data
"""
import xarray
from scipy.stats import linregress
import numpy as np
import pandas as pd

# applying spatial/spectral downsampling without sensor
def apply_downsampling(arr, sensor, spatial=True, spectral=False,
                       normalise=False):
    """Apply only the spectral and spatial downsampling steps
    of the sensor.

    Args:
        arr: array of radiance or reflectance values
        sensor: a fitted instance of a TdiCmos
        spatial (bool): applies spatial downsampling if True
        spectral (bool): applies spectral downsampling if True
        normalise (bool): Normalises each band by the integral of the spectral
            response function - use if input arr is reflectance

    Returns:
        New array of values
    """
    if spectral:
        arr = sensor.spectral_response.transform(arr, normalise)
    if spatial:
        arr = sensor.steps[1][1](arr,
                                 sensor.psf_fwhm,
                                 sensor.ground_sample_distance)
    return arr

def sensor_correction_experiment(toa_radiance, sensor, mask=None):
    """Performs a calibration experiment based on an input radiance dataset
    and a sensor to give the correction factors needed to convert DN
    to apparent radiance. Note that this assumes a linear relationship
    between DN and radiance.

    Args:
        toa_radiance: At-sensor (top-of-atmosphere) Radiance
        sensor : A valid sensor instance
        mask (bool):

    Returns:
        per-band calibration coefficients (m,c) for linear
        model to convert DN to radiance
    """
    # applying linear regression model
    def find_coeffs(radiance, sensor_output):
        # create vectors
        y = radiance.values.ravel()
        y = y[~np.isnan(y)]
        x = sensor_output.values.ravel()
        x = x[~np.isnan(x)]
        res = linregress(x, y)
        return res.slope, res.intercept

    # converting linregress coefs to xarray
    def corrections_to_xarray(corrs):
        df = pd.DataFrame(corrs).T
        m = xarray.DataArray(df[0], dims='band', name='m')
        c = xarray.DataArray(df[1], dims='band', name='c')
        return xarray.Dataset({'m':m,'c':c})

    # Sensor equivalent radiance
    dn = sensor.fit_transform(toa_radiance)
    L_s = apply_downsampling(toa_radiance,
                            sensor,
                            spatial=True, spectral=False)

    # iterate bands and apply linear regression to each band
    out = {}
    for i, b in enumerate(dn.band_name.values):
        if mask is not None:
            out[i] = find_coeffs(L_s.isel(band=i).where(mask),
                                    dn.isel(band=i).where(mask))
        else:
            out[i] = find_coeffs(L_s.isel(band=i),
                                    dn.isel(band=i))

    return corrections_to_xarray(out)

class LinearRadiometricCorrection(object):
    """Linear radiometric correction using precomputed coefficients
    """
    def __init__(self, coef_array=None, coef_filepath=None):
        """
        Args:
            coef_array (array-like): xarray.Dataset with 'm' and 'c' variables
            coef_filepath (str): path to NetCDF Dataset with 'm' and 'c'
                variables
        """
        if coef_array:
            self.corrections = coef_array
        elif coef_filepath:
            self.corrections = self._load_corrections(coef_filepath)
        else:
            raise ValueError('coef_array or coef_filepath required')

    def _load_corrections(self, path):
        return xarray.load_dataset(path)

    def transform(self, signal):
        """Applies linear correction

        Args:
            signal: Digitial Number array with bands matching those in
                correction coefficient array

        Returns:
            Radiance values
        """
        return (signal * self.corrections['m'] ) + self.corrections['c']
