"""
Checking decorators for different function types
"""
import functools


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
