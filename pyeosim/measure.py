"""
Scoring and evaluation metrics
"""
import numpy as np
from ._decorators import spectral_index


def RMSE(x0, x1):
    """
    Calculates the RMSE of 2 xarrays, ignoring NA values

    Parameters
    ----------
    x0 : xarray.DataArray
        first array
    x1 : xarray.DataArray
        second array
    """
    _dif = (x0-x1)**2
    # count all non-NAs
    T = float(_dif.count())
    return np.sqrt(float(_dif.sum())/T)


@spectral_index
def NDVI(da, case):
    """
    Normalized Difference Vegetation Index

    hyperspectral case:
    r is 650-680 mean
    ir is 790-899 mean

    Parameters
    ----------
    da : xarray.DataArray
        reflectance with bands or wavelengths
    case : str
        text flag indicating input type
    """
    if case == 'hyperspectral':
        _r = da.sel(wavelength=slice(650, 680)).mean('wavelength')
        _ir = da.sel(wavelength=slice(790, 899)).mean('wavelength')
    if case == 'multispectral':
        _r = da.sel(band_name='B4')
        _ir = da.sel(band_name='B8')
    return (_ir-_r)/(_ir+_r)
