"""
Module for atmospheric transformers. These can either work from a lookup table
or using the GenericTransformer template
"""

from ._atmosphere import LUT
from .datasets import DATA_PATHS


class Test6S(LUT):
    """
    Test Look-up table for 6SV-style lookup tables
    """

    def __init__(self):
        super().__init__(LUT_path=DATA_PATHS['TEST_LUT'])
        # 6SV output is in W m-2 sr-1 micron-1 so convert to
        # W m-2 sr-1 nanometer-1 by dividing by 1000
        self.LUT = self.LUT / 1000
