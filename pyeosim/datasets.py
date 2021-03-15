"""
Datasets and data load funcs used in testing and specific supported sensors
"""
import numpy as np
import os
from os.path import join as pjoin
import xarray
import pandas


def names():
    """
    List all datasets
    """
    return list(DATA_PATHS.keys())


def dload(name):
    """
    Load dataset by name.

    Parameters
    ----------
    name : str
        dataset identifier
    """
    if name in ['SRF_SENTINEL_2']:
        return _load_S2_spectra(DATA_PATHS[name])

    if name in ['SRF_SUPERDOVE']:
        return _load_superdove(DATA_PATHS[name])

    return _load_srf(DATA_PATHS[name])


# private funcs
def _load_S2_spectra(s2_path, tolerance=0.05):
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
                                  name='response')
            ar.attrs['Band'] = band
            out[sat][band] = ar
    return out


def _load_srf(fpath):
    df = pandas.read_csv(fpath, header=None)
    wlen = df[0].astype(float)
    resp = df[1].astype(float)
    new_wlen = np.arange(np.ceil(wlen.min()), np.floor(wlen.max()), 1)
    ar = xarray.DataArray(resp, [('wavelength', wlen)]).interp(
        wavelength=new_wlen)
    ar.name = 'response'
    return ar


def _load_superdove(fpath):
    df = pandas.read_csv(fpath)
    out = {}
    for i in range(0,16,2):
        dat = df.iloc[1:,i:i+2].astype(float)
        name = df.iloc[:,i].name
        dat = dat.drop_duplicates(name)

        new = xarray.DataArray(dat.iloc[:,1].values,
                               coords=[('wavelength', dat.iloc[:,0].values)])
        out[name] = new.dropna(
            'wavelength').interp(wavelength=np.arange(400,1000))
    return out


# paths for all data resources installed alongside
HERE = os.path.abspath(os.path.dirname(__file__))
DATA_PATHS = {
    'TEST_HSI': os.path.abspath(pjoin(HERE, 'data',
                                      'test_hyperspectral_vnir.nc')),
    'TEST_HSI_LARGE': os.path.abspath(pjoin(HERE, 'data',
                                      'test_hyperspectral_vnir_large.nc')),
    'SRF_SENTINEL_2': pjoin(HERE, 'data',
                            'srf_sentinel_2.csv'),
    'SOLAR_SPECTRUM_ASTMG173': pjoin(HERE, 'data',
                                     'solar_spectrum_ASTMG173.csv'),
    'CCD_QE_DD_BACK': pjoin(HERE, 'data',
                            'ccd_qe_dd_back.csv'),
    'CCD_QE_STD_BACK': pjoin(HERE, 'data',
                             'ccd_qe_std_back.csv'),
    'TEST_LUT': pjoin(HERE, 'data', 'test_6s.LUT'),
    'TDI_QE_BACK': pjoin(HERE, 'data', 'teledyne_cmos_qe_back.csv'),
    'SRF_SUPERDOVE': pjoin(HERE, 'data', 'srf_superdove.csv')
    }
