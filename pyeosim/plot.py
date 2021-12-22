"""
Plotting utilities
"""

import numpy as np


def rgb(open_datafile, return_array=False, **imshow_kwargs):
    """Make an RGB image from a 3 channel array.

    Get values and apply dtype conversion and histogram equalization

    Args:
        open_datafile: raster dataset or DataArray
        return_array (bool): if True, returns a numpy array of image plot
        **imshow_kwargs: passed directly to matplotlib imshow
    """
    ar = open_datafile.transpose('y', 'x', 'band').copy()
    ar.values = _image_histogram_equalization(ar.values)
    if return_array:
        return ar
    ar.plot.imshow(**imshow_kwargs)


def _image_histogram_equalization(image, number_bins=256, cdf=None, bins=None):
    # get image histogram
    if cdf is None:
        image_histogram, bins = np.histogram(image.flatten(), number_bins)
        cdf = image_histogram.cumsum()  # cumulative distribution function
        cdf = 255 * cdf / cdf[-1]  # normalize

    # use linear interpolation of cdf to find new pixel values
    image_equalized = np.interp(image.flatten(), bins[:-1], cdf)
    return image_equalized.reshape(image.shape).astype('uint8')
