"""
Image post processing pipeline
"""

from ._pipeline import GenericTransformer
from ._postprocess import *
from .sensor import TdiCMOS
from copy import deepcopy


class CmosReflectance(GenericTransformer):
    """
    Wrapper to generate corrected images
    """

    def __init__(self, cmos_sensor, reference_signal):
        """

        Parameters
        ----------
        cmos_sensor : georis.sensor.TdiCMOS
            sensor instance
        reference_signal : xarray.DataArray
            100% reference signal

        """
        super().__init__()
        self.cmos_sensor = cmos_sensor
        self.reference_signal = reference_signal

    def fit(self, signal):
        self.cmos_sensor.fit(signal)
        # make a new TdiCMOS instance to represent the non-image
        # signal - this is already fitted to the original signal
        self.dark_cmos_sensor = deepcopy(self.cmos_sensor)
        # copy column fpn from original to dark region and update
        self.dark_cmos_sensor.column_offset_FPN =\
            self.cmos_sensor.column_offset_FPN
        self.dark_cmos_sensor.set_steps()
        self.set_steps()
        self._fitted = True

    def set_steps(self):
        self.steps = [


        ]
