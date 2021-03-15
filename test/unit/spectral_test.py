from pyeosim import datasets
from pyeosim import spectral

import xarray
import inspect


def test_equivalent():
    test_data = xarray.open_dataset(datasets.TEST_HSI)
    test_data_x = xarray.open_dataset(datasets.TEST_HSI, chunks=2)
    all_classes = _get_public_classes(spectral)
    test_data = test_data.swap_dims({'band': 'wavelength'})
    test_data_x = test_data_x.swap_dims({'band': 'wavelength'})
    assert len(all_classes) > 0
    for spc in all_classes:
        model = spc[1]()
        out = test_data.reflectance.pipe(model.transform)
        outx = test_data_x.reflectance.pipe(model.transform).compute()
        assert out.sum() > 1
        assert outx.sum() > 1
        assert out.sum() == outx.sum()


def _get_public_classes(module):
    clsmembers = inspect.getmembers(module,
                                    inspect.isclass)
    return [x for x in clsmembers if not x[0].startswith('_')]
