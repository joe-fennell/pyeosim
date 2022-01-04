"""Imager and System Simulation
"""

# import xarray
from ._decorators import return_equal_xarray
from .datasets import _dload
from .spatial import gaussian_isotropic
from .spectral import *
from ._pipeline import GenericTransformer

import numpy
import xarray


class TdiCmos(GenericTransformer):
    """Time Delay Integration CMOS generic transformer.

    Note:
        This version does not apply spectral resampling and relies on band
        integration happening in the atmospheric simulation
    """

    def __init__(self, sensor_altitude=5e5, sensor_ground_speed=7000,
                 ground_sample_distance=2, lens_diameter=.1, psf_fwhm=4,
                 TDI_rows=32, pix_per_row=8000, sensor_width=82.2,
                 pixel_area=10**2, spectral_response='TreeView_3()',
                 quantum_efficiency='TDI_QE_BACK', full_well=3e4,
                 prnu_factor=.01, dark_current=570, dark_factor=.01,
                 offset_factor=.001, ccd_vref=3.1, sense_node_gain=5,
                 read_noise=20, adc_vref=.5, bit_depth=14,
                 store_steps=False, apply_spatial_resampling=True):
        """
        Args:
            sensor_altitude (float): altitude above source in metres
            sensor_ground_speed (float): relative ground speed of sensor in m/s
            ground_sample_distance (float): Ground sampling distance of sensor
            lens_diameter (float): diameter of lens in metres
            psf_fwhm (float): Point Spread Function Full-Width at Half Maximum
                in ground units (metres)
            TDI_rows (int): number of rows integrated over per channel
            pix_per_row (int): number of pixels per row
            sensor_width (float): width of sensor imaging area in mm
            pixel_area (float): per pixel light-absorbing area in micron2
            spectral_response: spectral transformer instance
            quantum_efficiency: Quantum efficiency of sensor as a float,
                name of internal dataset to use or xarray.DataArray with a
                'wavelength' dimension covering the wavelength range of the
                spectral_response or 'bands' dimension matching the bands in
                spectral_response
            full_well (int): the maximum capacity of a pixel in electrons
            prnu_factor (float): Photo response non uniformity factor
            dark_current (float): dark current in electrons/second/pixel
            dark_factor (float): dark current fixed pattern noise factor
                (DSNU factor)
            offset_factor (float): column offset factor
            ccd_vref (float): reference voltage of voltage sensor in volts
            sense_node_gain (float): gain of sense node in microvolts/electron
            read_noise (int): Read noise in e-
            adc_vref (float): reference voltage of ADC in volts
            bit_depth (int): bit depth of the ADC
                store_steps (float): if True, all intermediate steps will be
                stored in the step_outputs attribute. Default is False
            apply_spatial_resampling (float): if True, spatial response
                function will be applied. Use False is supplying non-imaging
                data to the simulation
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
        if type(spectral_response) == str:
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
        self.apply_spatial_resampling = apply_spatial_resampling
        # calculate all derived properties
        self.update_derived_params()

    def update_derived_params(self):
        """Re-run all derived parameter calculations.

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
                               _dload(self.quantum_efficiency))
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
        """Precomputes the system constant states.

        This precomputes states for system constants (e.g. the fixed-
        pattern noise).

        Args:
            signal: A Top-Of-Atmosphere radiance dataset integrated over the
                bandpass response

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
        if self.dark_factor == 0:
            self.DSNU = 0
        else:
            self.DSNU = DSNU(ones,
                             self.dark_current,
                             self.integration_time,
                             self.dark_factor).compute()
        # generate a photon response fixed pattern
        if self.prnu_factor == 0:
            self.PRNU = 0
        else:
            self.PRNU = PRNU(ones, self.prnu_factor).compute()
        # generate a column offset fixed pattern which is constant for all
        # bands as each band is horizontally aligned on the same sensor
        if self.offset_factor == 0:
            self.column_offset_FPN = 0
        else:
            self.column_offset_FPN = CONU(ones, self.offset_factor)
        self._set_steps()
        self._fitted = True

    def transform(self, signal):
        """Runs simulation on signal Top-Of-Atmosphere radiance dataset.

        This uses the precomputed states for system constants (e.g. the fixed-
        pattern noise) and generates a temporary random state for other
        parameters (e.g. the dark signal).

        Args:
            signal: A Top-Of-Atmosphere radiance dataset integrated over the
                bandpass response

        Returns:
            sensor digital number output
        """
        return super().transform(signal)

    def fit_transform(self, signal):
        """Runs simulation on signal Top-Of-Atmosphere radiance dataset.

        This both precomputes the system state constants (e.g. the fixed-
        pattern noise) and generates a temporary random state for other
        parameters (e.g. the dark signal).

        Args:
            signal: A Top-Of-Atmosphere radiance dataset integrated over the
                bandpass response

        Returns:
            sensor digital number output
        """
        return super().fit_transform(signal)

    def _make_ones(self, signal):
        # ones array of same shape as final signal
        # take first of any additional dims to guarantee 3D arr
        signal = signal.copy()
        for dim in signal.dims:
            if dim not in ['x', 'y', 'band']:
                signal = signal.isel({dim: 0})
        # Not needed in this version as does takes band, not wavelength
        # out = self.spectral_response.transform(signal)
        # only apply if spatial resampling flag set
        if self.apply_spatial_resampling:
            signal = gaussian_isotropic(signal, self.psf_fwhm,
                                     self.ground_sample_distance)
        return xarray.ones_like(signal)

    def _integration_time(self):
        # Calculate the necessary line rate
        t = self.sensor_ground_speed / self.ground_sample_distance
        # total integration over all rows
        return self.TDI_rows / t

    def _set_steps(self):
        # called via fit
        self.steps = [
            ('radiant energy to radiant flux', energy_to_quantity, {}),
            # The spatial resampling step will be added here
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
        if self.apply_spatial_resampling:
            self.steps.insert(1, ('apply spatial resampling',
                                  gaussian_isotropic,
                                  {'psf_fwhm': self.psf_fwhm,
                                   'ground_sample_distance':
                                       self.ground_sample_distance})
                              )


def add_column_offset(voltage, offset):
    """Adds the column offset noise to voltage array.

    Args:
        voltage: voltage array
        offset: column voltage offset array

    Returns:
        voltage array
    """
    return voltage + offset


def add_dark_signal(electrons, dark_current, integration_time, dsnu):
    """Adds dark current to electron count.

    Args:
        electrons: electron array
        dark_current (float): dark current in e-/s/pixel
        integration_time (float): integration period in seconds
        dsnu : dark current fixed pattern noise array

    Returns:
        electron count array
    """
    # Generate dark electron quantity
    n_dsig = dark_current * integration_time
    # make an empty copy of the electrons array to add the shot noise
    n_dsig_shot = electrons.copy()
    # convert mean I_dsig mean electron vals to poisson random process
    n_dsig_shot.values = numpy.random.poisson(n_dsig, size=electrons.shape)
    # add DSNU noise
    n_dark = n_dsig_shot + n_dsig_shot * dsnu
    return (electrons + n_dark).round()


def add_gaussian_noise(electron_count, sigma):
    """Adds a zero-mean gaussian to electron count

    Args:
        electron_count: electron count array
        sigma: SD of noise in electrons

    Returns:
        electron count array
    """
    gaus_noise = numpy.random.normal(0, sigma,
                                     size=electron_count.shape)
    return (electron_count + gaus_noise).round()


def add_photon_noise(photon_count):
    """Adds photon shot noise.

    Assumes signal is a Poisson random process.

    Args:
        photon_count: mean photon count array

    Returns:
        photon count array
    """
    @return_equal_xarray
    def poisson_rv(ar):
        return numpy.random.poisson(ar)
    # convert photon means to a poisson random process
    return photon_count.map_blocks(poisson_rv)


def add_prnu(electrons, prnu):
    """Adds Photon Response Non-Uniformity to electrons.

    Args:
        electrons : electron count array
        prnu : fixed pattern noise array

    Returns:
        electron count array
    """
    # add photo non-uniformity
    return (electrons + electrons * prnu).round()


def CONU(ones, offset_factor):
    """Generates column offset non uniformity array

    Args:
        ones: array of ones of same shape as image
        offset_factor (float): column offset factor

    Returns:
        Column offset non uniformity array with zero mean
    """
    fpn = ones.isel(band=0) * numpy.random.normal(0, offset_factor,
                                                  size=ones.isel(band=0).shape)
    # drop band and band_name coords if exist
    for var in ['band', 'band_name']:
        try:
            fpn = fpn.drop_vars(var)
        except ValueError:
            pass
    return fpn


def DSNU(ones, dark_current, integration_time, dark_factor):
    """Generates the dark signal non-uniformity array.

    This assumes a log-Normal distribution suitable for short integration times
    of less than 100 seconds

    Note:
        dark_current should be per actual pixel, not per pixel subarray

    Args:
        ones: array of ones in same shape as image
        dark_current (float): sensor level mean dark current in e-/s/pixel
        integration_time (float): integration time of exposure
        dark_factor (float): dark current FPN factor range (0...1)

    Returns:
        electron count array
    """
    sigma = integration_time * dark_current * dark_factor
    # multiply by ones array to convert numpy array to xarray
    fpn = ones * numpy.random.lognormal(0, sigma, size=ones.shape)
    # subtract off the mean to generate a zero-mean distribution
    return fpn - fpn.mean()


def electron_to_voltage(electron_count, v_ref, sense_node_gain,
                        full_well, read_noise):
    """Converts electrons (charge) to voltage.

    Converts electrons (charge) to voltage by simulating a sense node with a
    fixed v_ref

    Args:
        electron_count: electron count array
        v_ref (float): sense node reference voltage
        sense_node_gain (float): gain in V/e-
        full_well (int): max number of electrons per pixel
        read_noise (int): read noise in electrons per channel

    Returns:
        sensor voltage array
    """
    # truncate at full well
    e = electron_count.where(electron_count < full_well, full_well)
    # add read noise
    e = add_gaussian_noise(e, read_noise).round()
    return v_ref * (e * sense_node_gain)


def energy_to_quantity(energy):
    """Convert energy to quanta.

    Args:
        energy: energy in J

    Returns:
        quantity or array of quantities
    """
    hc = 1.98644586e-25  # J m
    lambda_ = _nm_to_m(energy.wavelength)
    return (lambda_/hc) * energy


def photon_mean(flux, pixel_area, integration_time):
    """Convert sensor photon flux density to photon count at sensor.

    Args:
        flux: pixel flux density (photons/second/m2)
        pixel_area (float): area in microns^2
        integration_time (float): integration time in seconds

    Returns:
        photon count array
    """
    return flux * integration_time * _microns2_to_m2(pixel_area)


def photon_to_electron(photon_count, Q_E):
    """Convert photon count to electron count.

    Convert photon count to electron count using sensor quantum efficiency.

    Args:
        photon_count: photon count array
        Q_E: float in range 0-1 specifying photon conversion rate or xarray
            with same wavelength sampling as photon_count
    Returns:
        electron count array
    """
    return (photon_count * Q_E).round()


def PRNU(ones, prnu_factor):
    """Generates the Photo Response Non-Uniformity.

    Args:
        ones: array of ones in same shape as image
        prnu_factor: PRNU in range 0...1

    Returns:
        electron count array
    """
    return ones * numpy.random.normal(0, prnu_factor, size=ones.shape)


def radiance_to_irradiance_2(radiance, lens_diameter, focal_length):
    """Convert radiance (W m-2 sr-1) to irradiance of sensor based on lens.

    Note:
        Only an approximation for nadir sensor.
        Forshortening is not taken into account.
        Does not take into account cos^4 roll-off

    Args:
        radiance: at sensor radiance array
        lens_diameter (float): diameter of optic in metres
        focal_length (float): focal length of optic in metres

    Returns:
        sensor irradiance array
    """
    return radiance * (numpy.pi / 4) * ((lens_diameter / focal_length) ** 2)


def voltage_to_DN(voltage, v_ref, bit_depth):
    """Convert voltage to Digital Number via a linear ADC.

    Model a linear response ADC.

    Args:
        voltage: voltage array
        v_ref (float): reference voltage of ADC
        bit_depth (int): ADC bit depth

    Returns:
        Digital Number array
    """
    max_DN = numpy.int(2**bit_depth - 1)
    # DN = max_DN * (v_ref - voltage)
    DN = max_DN * (voltage / v_ref)
    # DN = (voltage * (2**bit_depth))/v_ref
    # DN = (voltage / v_ref) * (2**bit_depth)
    DN = DN.round().where(DN < max_DN, max_DN)
    return DN.round().where(DN > 0, 0)


# Internal functions
def _nm_to_m(x):
    """
    Convert nanometer (nm) to metre (m)
    """
    return x * 1e-9


def _microns2_to_m2(x):
    """
    Convert nanometer
    """
    return x * (1e-6**2)
