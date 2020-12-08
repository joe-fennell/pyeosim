"""
Image post processing pipeline
"""

from ._pipeline import GenericTransformer
from ._postprocess import *
from copy import deepcopy


class CmosReflectance(GenericTransformer):
    """
    Wrapper to generate corrected images
    """

    def __init__(self, cmos_sensor, reference_signal, store_steps=False):
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
        self.store_steps = store_steps

    def fit(self, signal):
        self.cmos_sensor.fit(signal)
        # make a new TdiCMOS instance to represent the non-image
        # signal - this is already fitted to the original signal
        self.dark_cmos_sensor = deepcopy(self.cmos_sensor)
        # copy column fpn from original to dark region and update
        self.dark_cmos_sensor.column_offset_FPN =\
            self.cmos_sensor.column_offset_FPN
        self.dark_cmos_sensor.set_steps()
        # generate a flat field image with 50% magnitude of full signal
        self.ff_frame = generate_ff_image(self.reference_signal * .8,
                                          self.cmos_sensor)
        self.set_steps()
        self._fitted = True

    def set_steps(self):
        self.steps = [
            ('Noise-corrected image', noise_corrected_signal,
             {'image_sensor': self.cmos_sensor,
              'dark_sensor': self.dark_cmos_sensor,
              'ff_frame': self.ff_frame}),
            ('Ref_max reflectance image', noise_corrected_reflectance,
             {'reference_signal': self.reference_signal,
              'image_sensor': self.cmos_sensor,
              'dark_sensor': self.dark_cmos_sensor,
              'ff_frame': self.ff_frame})
        ]
