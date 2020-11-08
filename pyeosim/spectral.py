"""
Sentinel 2 simulation tools
"""
from ._spectral import _SRF, _load_S2_spectra


class Sentinel2A(_SRF):
    """
    Simulates the Sentinel 2A satellite using SRF values published by ESA (ref:
    COPE-GSEG-EOPG-TN-15-0007) and estimates of the PSF properties
    (https://doi.org/10.3390/rs8060488).
    """

    def _load_srfs(self):
        spx = _load_S2_spectra()
        return spx['S2A']


class Sentinel2B(_SRF):
    """
    Simulates the Sentinel 2B satellite
    """

    def _load_srfs(self):
        spx = _load_S2_spectra()
        return spx['S2B']


class TreeView(_SRF):
    """
    Hypothetical satellite with spatial resolution of 5m (4x Sentinel 2) and
    spectral response of Sentinel 2 in VNIR. PSF is based on 4x Sentinel 2
    """

    def _load_srfs(self):
        spx = _load_S2_spectra()
        vals = spx['S2A']
        # drop SWIR
        vals.pop('B11')
        vals.pop('B12')
        return vals
