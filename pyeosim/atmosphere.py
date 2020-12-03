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
