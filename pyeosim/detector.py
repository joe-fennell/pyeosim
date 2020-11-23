"""
Classes for detector operations
"""

from ._detector import *
from .datasets import dload


class genericCCD(object):
    """
    A generic CCD
    """

    def __init__(self, integration_time, pixel_area=5.3, Q_e='CCD_QE_STD_BACK',
                 readout_noise=2.5, dark_noise=25, sensitivity=.588,
                 bit_depth=16, ADC_baseline=100):
        self.integration_time = integration_time
        self.Q_e = dload(Q_e)
        self.pixel_area = pixel_area
        self.readout_noise = readout_noise
        self.dark_noise = dark_noise
        self.sensitivity = sensitivity
        self.bit_depth = bit_depth
        self.ADC_baseline = ADC_baseline

    def fit(self, signal):
        pass

    def transform(self, signal):
        return signal.pipe(
            photon_mean, self.pixel_area, self.integration_time
        ).pipe(
            add_photon_noise
        ).pipe(
            photon_to_electron, self.Q_e
        ).pipe(
            add_gaussian_noise, self.dark_noise
        ).pipe(
            add_gaussian_noise, self.readout_noise
        ).pipe(
            electron_to_DN, self.sensitivity, self.bit_depth, self.ADC_baseline
        )
