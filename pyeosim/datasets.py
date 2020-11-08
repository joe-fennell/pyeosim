"""
Paths to datasets used in testing and specific supported sensors
"""
import os
from os.path import join as pjoin


# paths for all data resources installed alongside
HERE = os.path.abspath(os.path.dirname(__file__))
TEST_HSI = pjoin(HERE, 'data', 'test_hyperspectral_vnir.nc')
SRF_SENTINEL_2 = pjoin(HERE, 'data', 'srf_sentinel_2.csv')
SOLAR_SPECTRUM_ASTMG173 = pjoin(HERE, 'data', 'solar_spectrum_ASTMG173.csv')
