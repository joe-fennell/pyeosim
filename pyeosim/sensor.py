"""
Classes for detector operations
"""

from ._sensor import *  # sensor functions stored separately
from .datasets import dload
from .spatial import gaussian_isotropic
from .spectral import TreeView_1
from ._pipeline import GenericTransformer

import numpy
import xarray


class SimpleSensor(GenericTransformer):
    """
    A generic CCD
    """

    def __init__(self, integration_time=.1, pixel_area=5.3,
                 psf_fwhm=11, ground_sample_distance=5, sensor_altitude=5e5,
                 Q_e=.6, dark_noise=25, ccd_vref=1,
                 sense_node_gain=5, full_well=100000, adc_vref=1, adc_gain=1,
                 bit_depth=12, store_steps=False):
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
            dark signal in electrons/second/pixel
        ccd_vref : float
            reference voltage of voltage sensor in volts
        sense_node_gain : float
            gain of sense node in microvolts/electron
        full_well : int
            the maximum capacity of a pixel in electrons
        adc_vref : float
            reference voltage of ADC in volts
        adc_gain : float
            ADC gain in DN/V
        bit_depth : int
            bit depth of the ADC
        store_steps : bool
            if True, all intermediate steps will be stored in the step_outputs
            attribute
        """
        super().__init__()
        self.integration_time = integration_time
        self.ground_sample_distance = ground_sample_distance
        self.psf_fwhm = psf_fwhm
        self.spectral_response = TreeView_1()
        # if string assume a dataset and resample
        if type(Q_e) == str:
            # load dataset from file
            # calculate weighted mean for each band
            self.Q_e = band_Qe(self.spectral_response.srfs, dload(Q_e))
        else:
            self.Q_e = Q_e
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
        self.fit(None)

    def _set_steps(self):
        self.steps = [
            ('irradiance per original pixel', radiance_to_irradiance,
             {'altitude': self.sensor_altitude}),
            ('irradiance to flux', irradiance_to_flux, {}),
            ('flux at CCD', self.spectral_response.transform, {}),
            ('flux at resampled pixel', gaussian_isotropic,
             {'psf_fwhm': self.psf_fwhm,
              'ground_sample_distance': self.ground_sample_distance}),
            ('flux to quanta', photon_mean,
             {'pixel_area': self.pixel_area,
              'integration_time': self.integration_time}),
            ('photon noise', add_photon_noise, {}),
            ('photon to electron', photon_to_electron,
             {'Q_e': self.Q_e}),
            ('dark current noise', add_gaussian_noise,
             {'sigma': self.dark_noise}),
            ('electron to voltage', electron_to_voltage,
             {'v_ref': self.ccd_vref,
              'sense_node_gain': self.sense_node_gain,
              'full_well': self.full_well}),
            ('voltage to DN', voltage_to_DN,
             {'v_ref': self.adc_vref,
              'adc_gain': self.adc_gain,
              'bit_depth': self.bit_depth})
        ]


class TeledyneCMOS(GenericTransformer):
    """
    Teledyne CMOS sensor assuming a monolithic multispectral architecture.

    """

    def __init__(self, integration_time=.1, pixel_area=10**2,
                 psf_fwhm=4, ground_sample_distance=2, sensor_altitude=5e5,
                 Q_e=.6, prnu_factor=.01, dark_current=570,
                 dark_factor=.01, offset_factor=.01, ccd_vref=5,
                 sense_node_gain=5, temperature=293, source_follower_gain=1,
                 full_well=30000, adc_vref=5, adc_gain=1, bit_depth=12,
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
        dark_current : float
            dark current in electrons/second/pixel
        dark_factor : float
            dark current fixed pattern noise factor
        offset_factor : float
            column offset factor
        ccd_vref : float
            reference voltage of voltage sensor in volts
        sense_node_gain : float
            gain of sense node in microvolts/electron
        temperature : float
            temperature of amplifier in kelvin
        source_follower_gain : float
            gain of source follower amplifier in microvolts/electron
        full_well : int
            the maximum capacity of a pixel in electrons
        adc_vref : float
            reference voltage of ADC in volts
        adc_gain : float
            ADC gain in DN/V
        bit_depth : int
            bit depth of the ADC
        store_steps : bool
            if True, all intermediate steps will be stored in the step_outputs
            attribute
        """
        super().__init__()
        self.integration_time = integration_time
        self.ground_sample_distance = ground_sample_distance
        self.psf_fwhm = psf_fwhm
        self.spectral_response = TreeView_1()
        # if string assume a dataset and resample
        if type(Q_e) == str:
            # load dataset from file
            # calculate weighted mean for each band and scale to CMOS-like QE
            self.Q_e = band_Qe(self.spectral_response.srfs, dload(Q_e)) * .8
        else:
            self.Q_e = Q_e
        self.prnu_factor = prnu_factor
        self.pixel_area = pixel_area
        self.dark_current = dark_current
        self.dark_factor = dark_factor
        self.offset_factor = offset_factor
        self.ccd_vref = ccd_vref
        self.sense_node_gain = sense_node_gain * 1e-6  # uV to V
        self.temperature = temperature
        self.source_follower_gain = source_follower_gain
        self.full_well = full_well
        self.adc_vref = adc_vref
        self.adc_gain = adc_gain
        self.bit_depth = bit_depth
        self.sensor_altitude = sensor_altitude
        self.store_steps = store_steps

    def fit(self, signal):
        """
        Sets all steps in the model. Generates fixed pattern noise that is
        preserved until transform is re-run.

        Parameters
        ----------
        signal : xarray.DataArray
            An xarray instance of measurements
        """
        # ones array of same shape as final signal

        def make_ones(self, signal):
            # take first of any additional dims to guarantee 3D arr
            signal = signal.copy()
            for dim in signal.dims:
                if dim not in ['x', 'y', 'wavelength']:
                    signal = signal.isel({dim: 0})
            out = self.spectral_response.transform(signal)
            out = gaussian_isotropic(out, self.psf_fwhm,
                                     self.ground_sample_distance)
            return xarray.ones_like(out)

        ones = make_ones(self, signal)
        # generate a dark current fixed pattern for imaging region
        self.sensor_dc_FPN = DCFPN(ones,
                                   self.dark_current,
                                   self.integration_time,
                                   self.dark_factor).compute()
        # generate a photon response fixed pattern
        self.PRNU = PRNU(ones, self.prnu_factor).compute()
        # generate a column offset fixed pattern
        self.column_offset_FPN = ones.isel(y=0, band=0) * numpy.random.normal(
            0, self.offset_factor**2, size=ones.isel(y=0, band=0).shape)
        self._set_steps()
        self._fitted = True

    def _set_steps(self):
        # called via fit
        self.steps = [
            ('irradiance per original pixel', radiance_to_irradiance,
             {'altitude': self.sensor_altitude}),
            ('irradiance to flux', irradiance_to_flux, {}),
            ('flux at CCD', self.spectral_response.transform, {}),
            ('flux at resampled pixel', gaussian_isotropic,
             {'psf_fwhm': self.psf_fwhm,
              'ground_sample_distance': self.ground_sample_distance}),
            ('flux to quanta', photon_mean,
             {'pixel_area': self.pixel_area,
              'integration_time': self.integration_time}),
            ('photon noise', add_photon_noise, {}),
            ('photon to electron', photon_to_electron,
             {'Q_e': self.Q_e}),
            ('photon FP noise', add_prnu,
             {'prnu': self.PRNU}),
            ('dark current noise', add_dark_noise,
             {'dark_current': self.dark_current,
              'integration_time': self.integration_time,
              'dcfpn': self.sensor_dc_FPN}),
            ('electron to voltage', electron_to_voltage_ktc,
             {'v_ref': self.ccd_vref,
              'sense_node_gain': self.sense_node_gain,
              'full_well': self.full_well,
              'temperature': self.temperature}),
            ('column offset noise', add_column_offset,
             {'offset': self.column_offset_FPN}),
            ('voltage to DN', voltage_to_DN,
             {'v_ref': self.adc_vref,
              'adc_gain': self.adc_gain,
              'bit_depth': self.bit_depth})
        ]
