from .datasets import SRF_SENTINEL_2
from .decorators import spectral_response

import logging
import numpy as np
import pandas
import xarray


# private classes
class _SRF(object):
    """
    A base simulator class - do not use
    """
    __slots__ = 'srfs'

    def __init__(self, reflectance_scale_factor=10000, band_names=None):
        self.srfs = self._load_srfs()
        self.reflectance_scale_factor = reflectance_scale_factor
        self.band_names = band_names

    def fit():
        """
        Method not used but retained for compatibility with SKLearn
        """
        pass

    def transform(self, signal):
        """
        Performs the convolution with the Spectral Response Function only

        Parameters
        ----------
        signal : list
            iterable of DataArrays covering the range of the SRF

        """
        # applies spectral response function to arrays
        def _clip(signal, range):
            out = signal.where((signal.wavelength >= range[0])
                               & (signal.wavelength <= range[1]), drop=True)
            if len(out) < 1:
                raise RuntimeError('Out of range: No values between {}nm and \
                                   {}nm'.format(*range))
            return out

        @spectral_response
        def _get_response(signal, sensor):
            # interpolate sensor to signal wavelengths (linear)
            sensor = sensor.interp(wavelength=signal.wavelength)
            # deal with any float rounding errors
            signal = signal.where(signal >= 0)
            # response
            response = (signal/self.reflectance_scale_factor) * sensor
            # integrate under response spectrum
            # estimate response for a 100% reflectance signal
            norm = sensor.integrate('wavelength')
            out = response.integrate('wavelength')/norm
            out.attrs = signal.attrs
            return out

        def _concatenate_results(results_dict):
            # concatenate bands along a new dimension
            new_ds = {}
            for k, ar in results_dict.items():
                new_ds[k] = ar.expand_dims('band')
            new_ds = xarray.concat(new_ds.values(), dim='band')
            new_ds = new_ds.assign_coords(
                {'band': np.arange(len(results_dict))})
            new_ds = new_ds.assign_coords(
                band_name=("band", list(results_dict.keys())))
            return new_ds

        responses = {}
        if self.band_names is not None:
            band_names = self.band_names
        else:
            band_names = self.srfs.keys()
        # have to iterate rather than map as out of range bands raise error
        for band in band_names:
            # find first dataset with full range coverage
            sensor = self.srfs[band]
            try:
                range = _min_max(sensor)
                responses[band] = _get_response(_clip(signal, range),
                                                _clip(sensor, range))
            except RuntimeError:
                logging.info('signal is not in range for band {}'.format(band))
                pass
        return _concatenate_results(responses)


# private funcs
def _load_S2_spectra(s2_path=SRF_SENTINEL_2, tolerance=0.05):
    # func for splitting col name string
    def _get_names(string):
        split = string.split('_')
        return split[0], split[-1]
    # reads all spectra into a dataset
    out = {'S2A': {}, 'S2B': {}}
    # read single CSV
    srfs = pandas.read_csv(s2_path)
    wlength = srfs.iloc[:, 0]  # first col is wavelength
    for s in srfs.columns[1:]:
        sat, band = _get_names(s)
        # skip out B1, B9, B10
        if band in ['B2', 'B3', 'B4', 'B8', 'B5',
                    'B6', 'B7', 'B8A', 'B11', 'B12']:
            # read as xarray and interpolate
            ar = xarray.DataArray(srfs[s],
                                  coords={'wavelength': wlength},
                                  dims='wavelength',
                                  name=s)
            ar.attrs['Band'] = band
            out[sat][band] = ar
    return out


def _min_max(resp, tolerance=.05):
    # max of spectrum
    max_ = (resp.wavelength[resp > tolerance]).max()
    # min of spectrum
    min_ = (resp.wavelength[resp > tolerance]).min()
    return (float(min_), float(max_))
