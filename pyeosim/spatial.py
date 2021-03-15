"""
Spatial response functions
"""
from ._decorators import spatial_response
from scipy import ndimage
import numpy as np
import xarray


@spatial_response
def gaussian_isotropic(signal, psf_fwhm, ground_sample_distance):
    """
    Apply downsampling transform assuming a gaussian PSF and known PSF full
    width at half maximum.

    Parameters
    ----------
    signal : xarray.DataArray
        a 2D xarray raster with valid 'res' in attributes

    psf_fwhm : float
        Full width at half maximum of gaussian kernel. Sentinel 2 10m has a
        value of 22.08m

    ground_sample_distance : float
        ground resolution of instrument

    Returns
    -------
    downsampled_array : xarray.DataArray
        2D xarray raster array at new resolution
    """
    # @return_equal_xarray
    # def apply_gauss_filter(x):
    #     return filters.gaussian(x, sigma=sigma)
    def apply_gaussian(ar):
        # wrapped func
        def gfilter(x):
            return ndimage.gaussian_filter1d(x, sigma)
        # apply in y dim, then in x dim
        ar = xarray.apply_ufunc(gfilter, ar,
                                input_core_dims=[['y']],
                                output_core_dims=[['y']])

        return xarray.apply_ufunc(gfilter, ar,
                                  input_core_dims=[['x']],
                                  output_core_dims=[['x']])

    # interpolate out NAs in each spatial axis
    signal = signal.transpose('y', 'x', ...)
    # get signal resolution
    try:
        res = signal.attrs['res'][0]
    except KeyError:
        try:
            res = signal.attrs['transform'][0]
        except KeyError:
            res = float(np.abs(signal.x[1]-signal.x[0]))
    # calculate sigma from FWHM
    sigma = (psf_fwhm/2.355)/res
    dx = int(ground_sample_distance/res)
    dy = dx
    # signal_new = signal.copy()
    # signal_new = signal.map_blocks(apply_gauss_filter)
    signal_new = apply_gaussian(signal)
    return signal_new.coarsen(x=dx, y=dy, boundary='pad').mean()
