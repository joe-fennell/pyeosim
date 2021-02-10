"""
Classes for detector operations
"""

from ._sensor import *  # sensor functions stored separately
from .datasets import dload
from .spatial import gaussian_isotropic
from .spectral import *
from ._pipeline import GenericTransformer

import numpy
import xarray


class TeledyneCMOS(GenericTransformer):
    """
    Teledyne CMOS sensor assuming a single multispectral array and perfect
    optics.

    """

    def __init__(self, sensor_altitude=5e5, sensor_ground_speed=7000,
                 ground_sample_distance=2, lens_diameter=.1, psf_fwhm=4,
                 TDI_rows=32, pix_per_row=8000, sensor_width=82.2,
                 pixel_area=10**2, spectral_response=None,
                 quantum_efficiency='TDI_QE_BACK', full_well=3e4,
                 prnu_factor=.01, dark_current=570, dark_factor=.01,
                 offset_factor=.001, ccd_vref=3.1, sense_node_gain=5,
                 read_noise=20, adc_vref=.5, bit_depth=14,
                 store_steps=False):
        """
        Parameters
        ----------
        sensor_altitude : float
            altitude above source in metres
        sensor_ground_speed : float
            relative ground speed of sensor in m/s
        ground_sample_distance : float
            Ground sampling distance of sensor
        lens_diameter : float
            diameter of lens in metres
        psf_fwhm : float
            Point Spread Function Full-Width at Half Maximum
            in ground units (metres)
        TDI_rows : int
            number of rows integrated over per channel
        pix_per_row : int
            number of pixels per row
        sensor_width : float
            width of sensor imaging area in mm
        pixel_area : float
            per pixel light-absorbing area in micron2
        spectral_response : object, None
            object with transform method that simulates the spectral response
            before the sensor. If None, default dataset is used.
        quantum_efficiency : float, str, xarray.DataArray
            Quantum efficiency of sensor as a float, name of internal dataset
            to use or xarray.DataArray with a 'wavelength' dimension covering
            the wavelength range of the spectral_response or 'bands' dimension
            matching the bands in spectral_response
        full_well : int
            the maximum capacity of a pixel in electrons
        prnu_factor : float
            Photo response non uniformity factor
        dark_current : float
            dark current in electrons/second/pixel
        dark_factor : float
            dark current fixed pattern noise factor (DSNU factor)
        offset_factor : float
            column offset factor
        ccd_vref : float
            reference voltage of voltage sensor in volts
        sense_node_gain : float
            gain of sense node in microvolts/electron
        read_noise : int
            Read noise in electrons
        adc_vref : float
            reference voltage of ADC in volts
        bit_depth : int
            bit depth of the ADC
        store_steps : bool
            if True, all intermediate steps will be stored in the step_outputs
            attribute
        """
        super().__init__()
        # Satellite properties
        self.sensor_altitude = sensor_altitude
        self.sensor_ground_speed = sensor_ground_speed
        self.ground_sample_distance = ground_sample_distance
        # Optic properties
        self.lens_diameter = lens_diameter
        self.psf_fwhm = psf_fwhm
        # Sensor geometric properties
        self.TDI_rows = TDI_rows
        self.pix_per_row = pix_per_row
        self.sensor_width = sensor_width
        self.pixel_area = pixel_area
        # Sensor bandpass properties
        if spectral_response is None:
            self.spectral_response = TreeView_1()
        elif type(spectral_response) == str:
            self.spectral_response = eval(spectral_response)
        else:
            self.spectral_response = spectral_response
        # Sensor properties
        self.quantum_efficiency = quantum_efficiency
        self.full_well = full_well
        self.prnu_factor = prnu_factor
        # dark_current is calculated from this later
        self.dark_current = dark_current
        self.dark_factor = dark_factor
        self.offset_factor = offset_factor
        self.ccd_vref = ccd_vref
        self.sense_node_gain = sense_node_gain
        self.read_noise = read_noise
        # ADC properties
        self.adc_vref = adc_vref
        self.bit_depth = bit_depth
        self.PRNU = None  # calculated during fit
        self.DSNU = None  # calculated during fit
        self.column_offset_FPN = None  # calculated during fit
        # Simulation properties
        self.store_steps = store_steps
        # calculate all derived properties
        self.update_derived_params()

    def update_derived_params(self):
        """
        Call after updating any parameter values to recalculate derived params
        """
        # Derived properties
        # integration time is n/line_rate where n is TDI_rows and line_rate
        # is sensor_ground_speed/ground_sample_distance
        self.integration_time = self._integration_time()
        # Retrieve quantum efficiencies
        # if string assume a named dataset so load and resample
        if type(self.quantum_efficiency) == str:
            # load dataset from file
            self.Q_E = band_QE(self.spectral_response.srfs,
                               dload(self.quantum_efficiency))
        elif type(self.quantum_efficiency) == list:
            # assume an ordered list of mean values
            self.Q_E = band_QE(self.spectral_response.srfs,
                               self.quantum_efficiency)
        else:
            # else try to use Q_E in original form
            self.Q_E = self.quantum_efficiency
        # dark current is per pixel so needs to be multiplied by number of
        # pixels in integration
        self._dark_current = self.dark_current * self.TDI_rows
        # calculate derived geometric parameters
        self.swath_width = self.ground_sample_distance * self.pix_per_row
        # calculate angular FoV
        self.afov = 2 * numpy.arctan2((self.swath_width * .5),
                                      self.sensor_altitude)
        # assume square pixels and calculate pix_width in metres
        self.pix_width = numpy.sqrt(self.pixel_area) / 1e6
        # calculate focal length of ideal optic for the FoV and sensor width
        # Note that this applies only to a rectilinear optic (thin optic or
        # mirror) focused at infinity
        self.focal_length = self.sensor_width / numpy.tan(self.afov) / 1e3
        self._sense_node_gain = self.sense_node_gain * 1e-6  # uV to V
        # set steps again to update values
        self._set_steps()

    def fit(self, signal):
        """
        Sets all steps in the model. Generates fixed pattern noise that is
        preserved until transform is re-run.

        Parameters
        ----------
        signal : xarray.DataArray
            An xarray instance of measurements
        """

        # As this is a line scanner, only 1D fixed pattern noise is needed
        # however because each image pixel is integrated over multiple physical
        # pixels, the noise properties destructively add over the rows in each
        # effective channel
        ones = self._make_ones(signal).isel(y=0)
        # generate a dark current fixed pattern for imaging region.
        # Assume the manufacturers quoted values are per physical pixel, so
        # divide by TDI_rows. Note that this is the DSNU, not Dark Signal -
        # this is already accounted for by multiplying the dark signal by the
        # number of TDI rows
        self.DSNU = DSNU(ones,
                         self.dark_current,
                         self.integration_time,
                         self.dark_factor).compute()
        # generate a photon response fixed pattern
        self.PRNU = PRNU(ones, self.prnu_factor).compute()
        # generate a column offset fixed pattern which is constant for all
        # bands as each band is horizontally aligned on the same sensor
        self.column_offset_FPN = CONU(ones, self.offset_factor)
        self._set_steps()
        self._fitted = True

    def _make_ones(self, signal):
        # ones array of same shape as final signal
        # take first of any additional dims to guarantee 3D arr
        signal = signal.copy()
        for dim in signal.dims:
            if dim not in ['x', 'y', 'wavelength']:
                signal = signal.isel({dim: 0})
        out = self.spectral_response.transform(signal)
        out = gaussian_isotropic(out, self.psf_fwhm,
                                 self.ground_sample_distance)
        return xarray.ones_like(out)

    def _integration_time(self):
        # Calculate the necessary line rate
        t = self.sensor_ground_speed / self.ground_sample_distance
        # total integration over all rows
        return self.TDI_rows / t

    def _set_steps(self):
        # called via fit
        self.steps = [
            ('radiant energy to radiant flux', energy_to_quantity, {}),
            ('apply bandpass filters', self.spectral_response.transform, {}),
            ('apply spatial resampling', gaussian_isotropic,
             {'psf_fwhm': self.psf_fwhm,
              'ground_sample_distance': self.ground_sample_distance}),
            ('radiant flux to flux density', radiance_to_irradiance_2,
             {'lens_diameter': self.lens_diameter,
              'focal_length': self.focal_length}),
            ('flux density to flux', photon_mean,
             {'pixel_area': self.pixel_area,
              'integration_time': self.integration_time}),
            ('add photon shot noise', add_photon_noise, {}),
            ('photon to electron', photon_to_electron,
             {'Q_E': self.Q_E}),
            ('add photo response non-uniformity', add_prnu,
             {'prnu': self.PRNU}),
            ('add dark signal', add_dark_signal,
             {'dark_current': self._dark_current,
              'integration_time': self.integration_time,
              'dsnu': self.DSNU}),
            ('electron to voltage', electron_to_voltage,
             {'v_ref': self.ccd_vref,
              'sense_node_gain': self._sense_node_gain,
              'full_well': self.full_well,
              'read_noise': self.read_noise}),
            ('add column offset noise', add_column_offset,
             {'offset': self.column_offset_FPN}),
            ('voltage to DN', voltage_to_DN,
             {'v_ref': self.adc_vref,
              'bit_depth': self.bit_depth})
        ]


class TCMOS_test(TeledyneCMOS):
    """
    Same as TeledyneCMOS but no spatial or spectral resample
    """

    def _set_steps(self):
        super()._set_steps()
        # remove gaussian resample step
        self.steps.pop(1)
        self.steps.pop(2)

    def _make_ones(self, signal):
        # ones array of same shape as final signal
        # take first of any additional dims to guarantee 3D arr
        signal = signal.copy()
        for dim in signal.dims:
            if dim not in ['x', 'y', 'wavelength']:
                signal = signal.isel({dim: 0})
        # out = self.spectral_response.transform(signal)
        return xarray.ones_like(signal)
