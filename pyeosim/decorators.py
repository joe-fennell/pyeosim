"""
Checking decorators for different function types
"""
import functools
import xarray


def spatial_response(func):
    """
    Checks first argument has ['x','y', ...] dims
    or exactly 2 dims and no named dims
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        # Xarray checks
        signal = args[0]
        if ('x' in signal.coords) & ('y' in signal.coords):
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
        if 'wavelength' in signal.coords:
            return func(*args, **kwargs)
        raise AttributeError('First argument does not have wavelength coord')
    return wrapper


def spatial_filter(func):
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
