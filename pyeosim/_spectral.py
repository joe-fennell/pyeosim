from ._decorators import spectral_response

import logging
import numpy as np
import xarray


# private classes
class _SRF(object):
    """
    A base simulator class - do not use
    """
    __slots__ = 'srfs'

    def __init__(self, band_names=None):
        self.band_names = band_names
        self.band_wavelengths = []
        self.srfs = self._load_srfs()

    def fit(self, signal):
        """
        Method not used but retained for compatibility with SKLearn
        """
        pass

    def transform(self, signal, normalise=False):
        """
        Performs the convolution with the Spectral Response Function only

        Parameters
        ----------
        signal : list
            iterable of DataArrays covering the range of the SRF
        normalise : bool, optional
            if True, response will be divided by integral of the band response

        """
        # applies spectral response function to arrays
        def _clip(signal, range):
            out = signal.sel(wavelength=slice(range[0], range[1]))
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
            response = signal * sensor
            # response = (signal/self.reflectance_scale_factor) * sensor
            # integrate under response spectrum
            # estimate response for a 100% reflectance signal]
            if normalise:
                norm = sensor.integrate('wavelength')
                out = response.integrate('wavelength')/norm

            out = response.integrate('wavelength')
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

    def to_6sv(self):
        """
        Interpolates the bandpass responses to 2.5nm interval and converts to
        micron for 6SV
        """
        def convert_srf(srf):
            def interp(min_wlen, max_wlen):
                # interpolate wavelengths
                return np.arange(min_wlen, max_wlen+ 2.5, 2.5)
            _min, _max = _min_max(srf)
            wlen = interp(_min, _max)
            return _min/1000, _max/1000, list(srf.interp(wavelength=wlen).values)

        out = {}
        for name, srf in self.srfs.items():
            out[name] = convert_srf(srf)
        return out


def bands_from_step_func(step_funcs, min_wavelength=400, max_wavelength=1000):
    """
    Generates a dictionary of spectral response functions using the mean
    value and bandwidth.

    Parameters
    ----------
    step_funcs : dict
        dictionary in format {'band_name': (central_wavelength, band_width,
        transmission)}
    min_wavelength : float
        minimum wavelength to evaluate in nm
    max_wavelength : float
        maximum wavelength to evaluate in nm

    Returns
    -------
    srfs : dict
        spectral response funcs in format {'band_name': xarray.DataArray}
    """
    out = {}
    new_wlen = np.arange(400, 1000)
    for name, params in step_funcs.items():
        try:
            trans = params[2]
        except IndexError:
            trans = 1
        min = params[0] - (params[1] / 2)
        max = params[0] + (params[1] / 2)
        new_r = np.zeros(len(new_wlen))
        # set all values >= min to 1
        new_r[new_wlen >= min] = trans
        # set all values > max back to 0
        new_r[new_wlen > max] = 0
        # make xarray and append to out dict
        out[name] = xarray.DataArray(new_r,
                                     coords=[('wavelength', new_wlen)])
    return out


def band_QE(SRFs, quantum_efficiency):
    """
    Calculates weighted mean of the quantum efficiency in
    each spectral channel.

    Parameters
    ----------
    SRFs : dict
        Spectral Response Function dictionary
    quantum_efficiency : xarray.DataArray, list
        Q_E with wavelength coordinate or a list of Q_Es per band

    Returns
    -------
    Q_E : xarray.DataArray
        Weighted mean quantum efficiency for each band
    """
    # numerical index for each band
    bands = np.arange(len(SRFs))
    QEs = []
    # if QE is supplied as a list per band:
    if type(quantum_efficiency) == list:
        if len(quantum_efficiency) == len(bands):
            QEs = quantum_efficiency
        else:
            raise ValueError('Quantum efficiency list length must equal\
                             number of bands')
    # if QE is given as a f(lambda) spectrum
    else:
        for srf in SRFs.values():
            weight = srf/srf.integrate('wavelength')
            QEs.append(float((quantum_efficiency * weight).sum()))
    return xarray.DataArray(QEs, [('band', bands)])


def _min_max(resp, tolerance=.05):
    # max of spectrum
    max_ = (resp.wavelength[resp > tolerance]).max()
    # min of spectrum
    min_ = (resp.wavelength[resp > tolerance]).min()
    return (float(min_), float(max_))
