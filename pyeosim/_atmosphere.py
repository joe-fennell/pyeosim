from ._decorators import spectral_response, reflectance_lookup
import xarray
import json


class LUT(object):
    """Lookup table for atmospheric simulation.

    Converts surface reflectance to at-sensor radiance values using
    lookup tables (generated separately).
    """

    def __init__(self, LUT_object=None, LUT_path=None, parameter_subsets=None,
                 chunks=None):
        """
        Args:
            LUT_object: An xarray lookup table which must
                have at least wavelength and rho.
            LUT_path: A path to an xarray lookup table
                which must have at least wavelength and rho.
            parameter_subsets: A dictionary specifying
                subsets of the xarray object. Refer to xarray docs for more
                information.
            chunks: integer chunking parameter passed to xarray.open_dataset
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
            set_lut(
                xarray.open_dataset(LUT_path, chunks=chunks).pixel_radiance)

        # Apply any subsetting
        if parameter_subsets is not None:
            set_lut(self.LUT.sel(**parameter_subsets))

    def fit(self, signal):
        pass

    def transform(self, signal):
        @spectral_response
        def apply_LUT(signal, LUT):
            if signal.wavelength.min() < LUT.wavelength.min():
                raise RuntimeError('Signal min wavelength less than LUT min.')
            if signal.wavelength.max() > LUT.wavelength.max():
                raise RuntimeError(
                    'Signal max wavelength greater than LUT max.')
            result = LUT.interp(rho=signal,
                                wavelength=signal.wavelength,
                                kwargs={'fill_value':None}).drop_vars('rho')
            return result.where(result>=0, other=0)
        meta = signal.attrs.copy()
        new = apply_LUT(signal, self.LUT)
        new_meta = {
            'input_signal_meta': meta,
            'atmospheric_simulation': json.dumps(self.LUT.attrs,
                                                 default=lambda o: '<Not stored>')
            # 'atmospheric_simulation': str(self.LUT.attrs)
            }
        new.attrs = new_meta
        return new
