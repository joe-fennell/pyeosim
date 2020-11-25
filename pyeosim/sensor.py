"""
Classes for detector operations
"""

from . import detector as dr
from .datasets import dload
from .spatial import gaussian_isotropic
from .spectral import TreeView_1
from ._pipeline import GenericTransformer


class LinearCCD(GenericTransformer):
    """
    A generic CCD
    """

    def __init__(self, integration_time=.1, pixel_area=5.3,
                 psf_fwhm=11, ground_sample_distance=5, sensor_altitude=5e5,
                 Q_e=.6, dark_noise=25, ccd_vref=1, sense_node_gain=5,
                 full_well=100000, adc_vref=1, adc_gain=1, bit_depth=12,
                 store_steps=False):
        """
        Parameters
        ----------
        integration_time : float
            frame integration time in seconds
        pixel_area : float
            per pixel light-absorbing area in micron2
        psf_fwhm : float
            Point Spread Function Full-Width at Half Maximum
            in ground units (metres)
        ground_sample_distance : float
            Ground sampling distance of sensor
        sensor_altitude : float
            altitude above source in metres
        Q_e : float
            Quantum efficiency of sensor
        dark_noise : float
            dark signal in electrons/second
        ccd_vref : float
            reference voltage of voltage sensor in volts
        sense_node_gain : float
            gain of sense node in microvolts/electron
        full_well : int
            the maximum capacity of a pixel in electrons
        adc_vref : float
            reference voltage of ADC in volts
        adc_gain : float
            ADC gain in volts/DN
        bit_depth : int
            bit depth of the ADC
        store_steps : bool
            if True, all intermediate steps will be stored in the step_outputs
            attribute
        """
        super().__init__()
        self.integration_time = integration_time
        # if string assume a dataset
        if type(Q_e) == str:
            self.Q_e = dload(Q_e)
        else:
            self.Q_e = Q_e
        self.ground_sample_distance = ground_sample_distance
        self.psf_fwhm = psf_fwhm
        self.spectral_response = TreeView_1()
        self.pixel_area = pixel_area
        self.dark_noise = dark_noise
        self.ccd_vref = ccd_vref
        self.sense_node_gain = sense_node_gain * 1e-6  # uV to V
        self.full_well = full_well
        self.adc_vref = adc_vref
        self.adc_gain = adc_gain
        self.bit_depth = bit_depth
        self.sensor_altitude = sensor_altitude
        self.store_steps = store_steps
        self.__set_steps()

    def __set_steps(self):
        self.steps = [
            ('radiance to quanta', dr.photon_mean,
             [self.pixel_area, self.integration_time]),
            ('quanta per original pixel', dr.radiance_to_irradiance,
             [self.sensor_altitude]),
            ('quanta at CCD', self.spectral_response.transform, []),
            ('quanta at resampled pixel', gaussian_isotropic,
             [self.psf_fwhm, self.ground_sample_distance]),
            ('photon noise', dr.add_photon_noise, []),
            ('photon to electron', dr.photon_to_electron, [self.Q_e]),
            ('dark current noise', dr.add_gaussian_noise, [self.dark_noise]),
            ('electron to voltage', dr.electron_to_voltage,
             [self.ccd_vref, self.sense_node_gain, self.full_well]),
            ('voltage to DN', dr.voltage_to_DN, [self.adc_vref,
                                                 self.adc_gain,
                                                 self.bit_depth])
        ]
