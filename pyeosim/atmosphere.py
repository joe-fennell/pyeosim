"""
Module for atmospheric transformers. These can either work from a lookup table
or using the GenericTransformer template
"""

from ._atmosphere import LUT
from .datasets import DATA_PATHS
import numpy as np
import os
import pandas
import xarray


class Test6S(LUT):
    """
    Test Look-up table for 6SV-style lookup tables
    """

    def __init__(self):
        super().__init__(LUT_path=DATA_PATHS['TEST_LUT'])
        # 6SV output is in W m-2 sr-1 micron-1 so convert to
        # W m-2 sr-1 nanometer-1 by dividing by 1000
        self.LUT = self.LUT / 1000
        self.LUT.attrs = {
            'latitude': 52.04,
            'longitude': 0.76,
            'datetimestring': '2020/06/22T12:00',
            'view_z': 0,
            'view_a': 0,
            'water_vapour': 2.93,
            'ozone': 0.319,
            'AOT550': 0.5,
            'visibility_km': 8.49
        }


def LUT_from_file(fpath, common_params={}):
    """
    Takes a directory of directories (one per scenario)
    Each file in the directory should be a CSV named with the rho (0...100)
    """
    out = []
    for sim in os.listdir(fpath):
        try:
            path = os.path.join(fpath, sim)
            files = os.listdir(path)
            # get rho from filename
            rhos = [float(x.split('_')[1].split('.')[0])/100 for x in files]
            # open first file
            _first = pandas.read_csv(os.path.join(path, files[0]))
            new = np.empty((len(rhos), len(_first['lambda']), 1))
            # iterate opening files to build array
            for i, f in enumerate(files):
                df = pandas.read_csv(os.path.join(path, f))
                new[i, :, 0] = df['radiance'].values
            # convert to xarray
            # convert to per nm firs
            new = xarray.DataArray(new/1000, coords=[
                ('rho', rhos),
                ('wavelength', 1000 * df['lambda'].values),
                ('scenario', [sim])
            ])
            out.append(new)
        except IndexError:
            pass
    ar = xarray.concat(out, 'scenario').interpolate_na(dim='wavelength')
    ar.attrs = common_params
    return LUT(ar)
