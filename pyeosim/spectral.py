"""
Sentinel 2 simulation tools
"""
from ._spectral import _SRF, bands_from_step_func, band_QE
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
        self.band_names = list(self.band_wavelengths.keys())
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
        self.band_names = list(self.band_wavelengths.keys())
        return spx['S2B']


class Sentinel2VNIR(_SRF):
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
        self.band_names = list(self.band_wavelengths.keys())
        return vals


class SuperDove(_SRF):
    """
    Super Dove estimated spectral response digitised from published plots
    """

    def _load_srfs(self):
        vals = dload('SRF_SUPERDOVE')
        self.band_wavelengths = {
            'CoastalBlue': 443.,
            'Blue': 490.,
            'Green_1': 531.,
            'Green_2': 565.,
            'Yellow': 611.,
            'Red': 666.,
            'RedEdge': 705.,
            'NIR': 865.
        }
        self.band_names = list(self.band_wavelengths.keys())
        return vals

class TreeView_1(_SRF):
    """
    Version 2 based on MRD specifications
    """

    def _load_srfs(self):
        band_defs = {
            'Clouds': (445, 20),
            'Carotenoids': (490, 40),
            'PRI_1': (531, 10),
            'PRI_2': (560, 20),
            'Chlorophyll_1': (620, 20),
            'Chlorophyll_2': (665, 30),
            'RedEdge_1': (700, 15),
            'RedEdge_2': (740, 15),
            'RedEdge_3': (780, 15),
            'NIR': (865, 30)
        }
        self.band_wavelengths = {k: v[0] for (k, v) in band_defs.items()}
        self.band_names = list(self.band_wavelengths.keys())
        return bands_from_step_func(band_defs)

class TreeView_2(_SRF):
    """
    Version 3 based on the a widened bandpass profile in the vis
    """

    def _load_srfs(self):
        band_defs = {
            'Aerosol': (440, 20),
            'Carotenoids': (480, 40),
            'PRI_1': (525, 50),
            'PRI_2': (585, 50),
            'Chlorophyll_1': (630, 40),
            'Chlorophyll_2': (670, 40),
            'RedEdge_1': (700, 15),
            'RedEdge_2': (740, 15),
            'RedEdge_3': (780, 15),
            'NIR': (865, 30)
        }
        self.band_wavelengths = {k: v[0] for (k, v) in band_defs.items()}
        self.band_names = list(self.band_wavelengths.keys())
        return bands_from_step_func(band_defs)
