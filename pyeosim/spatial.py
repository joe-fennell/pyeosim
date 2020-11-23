"""
Spatial response functions
"""
from .decorators import spatial_response, return_equal_xarray
from skimage import filters
import numpy as np


class GaussianIsotropic(object):
    """
    Parametric Gaussian PSF
    """

    def __init__(self, psf_fwhm=22.08, ground_sample_distance=10):
        """
        Apply downsampling transform assuming a gaussian PSF and known PSF full
        width at half maximum.

        NOTE: When using in dask lazy mode, it is important that the chunk size
        is large enough in the spatial domain (suggested > 220m ground for a
        PSF FWHM of 22m)

        Parameters
        ----------

        psf_fwhm : float, optional
            Full width at half maximum of gaussian kernel. Sentinel 2 10m has a
            value of 22.08m

        ground_sample_distance : float, optional
            ground resolution of instrument
        """
        self.psf_fwhm = psf_fwhm
        self.ground_sample_distance = ground_sample_distance

    def fit(self, signal):
        """
        Method not used but retained for compatibility with SKLearn
        """
        pass

    def transform(self, signal):
        """
        Apply downsampling transform assuming a gaussian PSF and known PSF full
        width at half maximum.

        Parameters
        ----------
        signal : xarray.DataArray
            a 2D xarray raster with valid 'res' in attributes

        Returns
        -------
        downsampled_array : xarray.DataArray
            2D xarray raster array at new resolution
        """
        return gaussian_isotropic(signal, self.psf_fwhm,
                                  self.ground_sample_distance)


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
    @return_equal_xarray
    def apply_gauss_filter(x):
        return filters.gaussian(x, sigma=sigma, multichannel=True)
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
    signal_new = signal.map_blocks(apply_gauss_filter)
    return signal_new.coarsen(x=dx, y=dy, boundary='pad').mean()
