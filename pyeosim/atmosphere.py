"""
Module for atmospheric transformers. These can either work from a lookup table
or using the GenericTransformer template
"""

from ._atmosphere import LUT
from .datasets import DATA_PATHS
import numpy as np
import os
import pandas
import  Py6S
import xarray


class Test6S(LUT):
    """
    Test Look-up table for 6SV-style lookup tables
    """

    def __init__(self):
        super().__init__(LUT_path=DATA_PATHS['TEST_LUT'])
        # 6SV output is in W m-2 sr-1 micron-1 so convert to
        # W m-2 sr-1 nanometer-1 by dividing by 1000
        self.LUT = self.LUT / 1000
        self.LUT.attrs = {
            'latitude': 52.04,
            'longitude': 0.76,
            'datetimestring': '2020/06/22T12:00',
            'view_z': 0,
            'view_a': 0,
            'water_vapour': 2.93,
            'ozone': 0.319,
            'AOT550': 0.5,
            'visibility_km': 8.49
        }


class SixSV_atmosphere(object):
    """
    Generates an atmospheric transformer for a specific 6SV atmosphere

    Methods
    -------
    fit(srf)
        Takes a spectral response function and computes the lookup table
    transform(signal)
        Computes the per-band TOA radiance of a signal array
    inverse_transform(signal)
        Computes the per-band BOA reflectance of a signal array
    """
    def __init__(self, SixS, srf=None):
        self.SixS = SixS
        self.srf = srf
        self._coefs = None
        self._fitted = False
        self._SixSV_outputs = {}
        self._6sv_outputs = []

        if srf is not None:
            self.fit(srf)

    def fit(self, srf):
        """
        Parameters
        ----------
        srf : object
            spectral response function object
        """
        def get_correction_coefs(fitted_6s):
            # Calculates the specific coefficients for this SixsV configuration and
            # bandpass combination
            s = fitted_6s
            # direct solar irradiance
            Edir = s.outputs.direct_solar_irradiance
            # diffuse solar irradiance
            Edif = s.outputs.diffuse_solar_irradiance
            # sum for total
            E = Edir + Edif
            # transmissivity
            # absorption transmissivity
            absorb  = s.outputs.trans['global_gas'].upward
            # scattering transmissivity
            scatter = s.outputs.trans['total_scattering'].upward
            # transmissivity (from surface to sensor)
            tau = absorb * scatter
            # Our corrected version:
            # path radiance
            Lp   = s.outputs.atmospheric_intrinsic_radiance #fw
            # correction coefficients for this configuration
            # i.e. surface_reflectance = (L - a) / b,
            #      where, L is at-sensor radiance
            a = Lp
            b = (tau * E) / np.pi
            bw = s.outputs.int_funct_filt
            return a, b, bw

        self.srf = srf
        _srfs = self.srf.to_6sv().values()
        a = []
        b = []
        bw = []
        raw = []
        for srf in _srfs:
            self.SixS.wavelength = Py6S.Wavelength(*srf)
            self.SixS.run()
            _a, _b, _bw = get_correction_coefs(self.SixS)
            a.append(_a)
            b.append(_b)
            bw.append(_bw)

        # convert to xarray
        a = xarray.DataArray(a, coords=[('band', np.arange(len(_srfs)))])
        b = xarray.DataArray(b, coords=[('band', np.arange(len(_srfs)))])
        bw = xarray.DataArray(bw, coords=[('band', np.arange(len(_srfs)))])
        self._coefs = xarray.Dataset({'a':a, 'b':b, 'filter_integral_micron':bw})

    def transform(self, signal):
        """
        Convert BOA reflectance to TOA radiance

        Parameters
        ----------
        signal : xarray.DataArray
            Spectral reflectance dataset
        """
        def reflectance_to_radiance(signal):
            # this returns the spectral radiance in W m-2 sr-1 micron-1
            # correct integral from micron to nm
            S_prime = self._coefs['filter_integral_micron']
            return ((signal * self._coefs['b']) + self._coefs['a']) * S_prime

        # convert spectral reflectance to reflectance across sensor band
        signal_sensor = self.srf.transform(signal, normalise=True)
        return reflectance_to_radiance(signal_sensor)

    def inverse_transform(self, signal):
        """
        Convert TOA radiance to BOA reflectance

        Parameters
        ----------
        signal : xarray.DataArray
            per band radiance dataset

        Outputs
        -------
        reflectance : xarray.DataArray

        """
        def radiance_to_reflectance(signal):
            S_prime = self._coefs['filter_integral_micron']
            return ((signal - self._coefs['a']) / self._coefs['b']) * S_prime

        # convert spectral reflectance to reflectance across sensor band
        signal_sensor = self.srf.transform(signal, normalise=True)
        return radiance_to_reflectance(signal_sensor)


def LUT_from_file(fpath, common_params={}):
    """
    Takes a directory of directories (one per scenario)
    Each file in the directory should be a CSV named with the rho (0...100)
    """
    out = []
    for sim in os.listdir(fpath):
        try:
            path = os.path.join(fpath, sim)
            files = os.listdir(path)
            # get rho from filename
            rhos = [float(x.split('_')[1].split('.')[0])/100 for x in files]
            # open first file
            _first = pandas.read_csv(os.path.join(path, files[0]))
            new = np.empty((len(rhos), len(_first['lambda']), 1))
            # iterate opening files to build array
            for i, f in enumerate(files):
                df = pandas.read_csv(os.path.join(path, f))
                new[i, :, 0] = df['radiance'].values
            # convert to xarray
            # convert to per nm firs
            new = xarray.DataArray(new/1000, coords=[
                ('rho', rhos),
                ('wavelength', 1000 * df['lambda'].values),
                ('scenario', [sim])
            ])
            out.append(new)
        except (IndexError, NotADirectoryError, KeyError):
            pass
    ar = xarray.concat(out, 'scenario').interpolate_na(dim='wavelength').sortby('rho')
    ar.attrs = common_params
    return LUT(ar)
