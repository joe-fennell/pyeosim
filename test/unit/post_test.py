from pyeosim.post import apply_downsampling, sensor_correction_experiment, LinearRadiometricCorrection
from pyeosim.datasets import DATA_PATHS
from pyeosim.imager import TdiCmos
from pyeosim.atmosphere import SixSV_atmosphere
from pyeosim.spectral import Sentinel2VNIR

import Py6S

import xarray

reflectance = xarray.load_dataarray(DATA_PATHS['TEST_HSI_LARGE'])
reflectance = reflectance.swap_dims({'band':'wavelength'}) / 10000
reflectance = reflectance.where(reflectance <= 1, other=1).sel(wavelength=slice(400,898))
radiance = SixSV_atmosphere(Py6S.SixS(), Sentinel2VNIR()).transform(reflectance)
sensor = TdiCmos(spectral_response=Sentinel2VNIR())
sensor.fit(radiance)

def test_spatial_downsampling():
    apply_downsampling(reflectance, sensor, spatial=True, spectral=False)

def test_spectral_downsampling():
    apply_downsampling(reflectance, sensor, spatial=False, spectral=True)

def test_normalisation():
    apply_downsampling(reflectance, sensor, spatial=False, spectral=True,
                       normalise=True)

def test_sensor_correction_experiment():
    coefs = sensor_correction_experiment(radiance, sensor)
    corrector = LinearRadiometricCorrection(coef_array=coefs)
    DNs = sensor.transform(radiance)
    corrector.transform(DNs)
