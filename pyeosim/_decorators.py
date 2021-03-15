"""
Checking decorators for different function types
"""
import functools
import numpy as np
import xarray


def spectral_index(func):
    """
    Adds case handling and limits to spectral index functions
    """
    @functools.wraps(func)
    def wrapper(da, **kwargs):
        if 'band' in da.dims:
            da = da.swap_dims({'band': 'band_name'})
            out = func(da, case='multispectral')
        if 'band_name' in da.dims:
            out = func(da, case='multispectral')
        if 'wavelength' in da.dims:
            out = func(da, case='hyperspectral')
        # mask the results if requested
        if 'min_val' in kwargs:
            out = out.where(out >= kwargs['min_val'])
        if 'max_val' in kwargs:
            out = out.where(out <= kwargs['max_val'])
        return out
    return wrapper


def reflectance_lookup(func):
    """
    Checks first argument has ['wavelength', 'rho'] dims
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # Xarray checks
        signal = args[0]
        if ('wavelength' in signal.dims) & ('rho' in signal.dims):
            return func(*args, **kwargs)
        raise AttributeError('First argument does not have spatial coord')
    return wrapper


def spatial_response(func):
    """
    Checks first argument has ['x','y', ...] dims
    or exactly 2 dims and no named dims
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # Xarray checks
        signal = args[0]
        if ('x' in signal.dims) & ('y' in signal.dims):
            return func(*args, **kwargs)
        raise AttributeError('First argument does not have spatial coord')
    return wrapper


def spectral_response(func):
    """
    Checks first argument has a 'wavelength' dim
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        signal = args[0]
        if 'wavelength' in signal.dims:
            return func(*args, **kwargs)
        raise AttributeError('First argument does not have wavelength coord')
    return wrapper


def sklearn_pixelwise(func, feature_dim=['band']):
    """
    Prepares xarray for sklearn and returns in xarray format
    """
    @functools.wraps(func)
    def wrapper(array, *args, **kwargs):
        # before passing to func make into 2D array
        subset = array.stack(obs=tuple(['y', 'x']),
                             feats=feature_dim).dropna('obs')
        mask = ~subset.isnull().isel(feats=0)
        # make a new copy of array
        new_vals = func(subset, *args, **kwargs)

        if new_vals.ndim > 1:
            new_shape = (subset.shape[0], new_vals.shape[1])
            out = xarray.DataArray(np.empty(new_shape) * np.nan,
                                   coords={'obs': subset.obs,
                                           'output': np.arange(new_shape[1])},
                                   dims=['obs', 'output'])
        else:
            new_shape = subset.shape[0]
            out = xarray.DataArray(np.empty(new_shape) * np.nan,
                                   coords={'obs': subset.obs},
                                   dims='obs')

        out[mask] = new_vals
        return out.unstack()
    return wrapper


def return_equal_xarray(func):
    """
    Transforms filter func for numpy arrays to return dataarrays of same size.
    Note output size must match input size. First arg must be input array
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        out = func(*args, **kwargs)
        # make a new copy of array
        out_ar = args[0].copy()
        out_ar.data = out
        return out_ar
    return wrapper
