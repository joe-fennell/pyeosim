from pyeosim.imager import TdiCmos
from pyeosim.spectral import TreeView_3, Sentinel2VNIR
from pyeosim.atmosphere import SixSV_atmosphere
from pyeosim.datasets import names, DATA_PATHS
from Py6S import SixS
import xarray

def test_standard_simulation():
    # load, set wavelength as dim coords and rescale to [0..1]
    reflectance = xarray.load_dataarray(DATA_PATHS['TEST_HSI_LARGE'])
    reflectance = reflectance.swap_dims({'band':'wavelength'}) / 10000

    # generate a 6SV atmsophere
    tv3 = TreeView_3()
    atm = SixSV_atmosphere(SixS(), tv3)

    # Generate a default TDI CMOS sensor
    TreeView_satellite = TdiCmos(spectral_response=TreeView_3())

    # Perform simulation
    imager_output = TreeView_satellite.fit_transform(atm.transform(reflectance))

    # plot outputs
    imager_output.plot(col='band', col_wrap=3, cmap='gray')

def test_sentinel_simulation():
    # load, set wavelength as dim coords and rescale to [0..1]
    reflectance = xarray.load_dataarray(DATA_PATHS['TEST_HSI_LARGE'])
    reflectance = reflectance.swap_dims({'band':'wavelength'}) / 10000

    # generate a 6SV atmsophere
    vnir = Sentinel2VNIR()
    print(vnir.to_6sv()['B2'])
    atm = SixSV_atmosphere(SixS(), vnir)

    # Generate a default TDI CMOS sensor
    TreeView_satellite = TdiCmos(spectral_response=Sentinel2VNIR())

    # Perform simulation
    imager_output = TreeView_satellite.fit_transform(atm.transform(reflectance))

    # plot outputs
    imager_output.plot(col='band', col_wrap=3, cmap='gray')
