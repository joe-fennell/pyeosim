"""
User module for accessing atmospheric transform
"""
from .decorators import spectral_response, reflectance_lookup
import xarray


class SixS(object):
    """
    Converts surface reflectance to at-sensor radiance values using the 6S
    lookup tables (generated separately).
    """

    def __init__(self, LUT_object=None, LUT_path=None, parameter_subsets=None,
                 chunks=None):
        """
        Parameters
        ----------
        LUT_path : str
            A path to an xarray lookup table which must have at least
            wavelength and rho
        """

        @reflectance_lookup
        def set_lut(LUT):
            self.LUT = LUT

        if (LUT_object is None) & (LUT_path is None):
            raise ValueError(
                'Either LUT_path or LUT_object needs to be supplied')

        self.parameter_subsets = parameter_subsets
        if LUT_object is not None:
            set_lut(LUT_object)
        else:
            set_lut(xarray.open_dataset(LUT_path, chunks=chunks).sel(
                **parameter_subsets))

    def fit(self, signal):
        pass

    def transform(self, signal):
        return apply_LUT(signal, self.LUT)


@spectral_response
def apply_LUT(signal, LUT):
    if signal.wavelength.min() < LUT.wavelength.min():
        raise RuntimeError('Signal min wavelength less than LUT min.')
    if signal.wavelength.max() > LUT.wavelength.max():
        raise RuntimeError('Signal max wavelength greater than LUT max.')
    return LUT.sel(rho=signal,
                   wavelength=signal.wavelength,
                   method='nearest').drop_vars('rho')
