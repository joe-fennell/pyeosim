"""
Sentinel 2 simulation tools
"""
from ._spectral import _SRF
from .datasets import dload


class Sentinel2A(_SRF):
    """
    Simulates the Sentinel 2A satellite using SRF values published by ESA (ref:
    COPE-GSEG-EOPG-TN-15-0007) and estimates of the PSF properties
    (https://doi.org/10.3390/rs8060488).
    """

    def _load_srfs(self):
        spx = dload('SRF_SENTINEL_2')
        self.band_wavelengths = {
            'B2': 492.,
            'B3': 560.,
            'B4': 664.,
            'B5': 704.,
            'B6': 740.,
            'B7': 780.,
            'B8': 832.,
            'B8A': 864.
        }
        return spx['S2A']


class Sentinel2B(_SRF):
    """
    Simulates the Sentinel 2B satellite
    """

    def _load_srfs(self):
        spx = dload('SRF_SENTINEL_2')
        self.band_wavelengths = {
            'B2': 492.,
            'B3': 560.,
            'B4': 664.,
            'B5': 704.,
            'B6': 740.,
            'B7': 780.,
            'B8': 832.,
            'B8A': 864.
        }
        return spx['S2B']


class TreeView_1(_SRF):
    """
    Hypothetical satellite with spectral responses equivalent to vis-NIR
    Sentinel 2
    """

    def _load_srfs(self):
        spx = dload('SRF_SENTINEL_2')
        vals = spx['S2A']
        # drop SWIR
        vals.pop('B11')
        vals.pop('B12')
        self.band_wavelengths = {
            'B2': 492.,
            'B3': 560.,
            'B4': 664.,
            'B5': 704.,
            'B6': 740.,
            'B7': 780.,
            'B8': 832.,
            'B8A': 864.
        }
        return vals
