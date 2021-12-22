# Python Earth Observation Simulation
This software was developed for the TreeView project and is intended for the simulation of Earth Observation imagery data.

## Installation
Best option is to use conda:
```
git clone git@github.com:joe-fennell/pyeosim.git
cd pyeosim
conda-build .
conda install --use-local
```

PIP install will work, but you should install Py6S first following the project's instructions
```
git clone git@github.com:joe-fennell/pyeosim.git
pip install .
```

## Getting Started
```python
from pyeosim.imager import TdiCmos
from pyeosim.spectral import TreeView_3
from pyeosim.atmosphere import SixSV_atmosphere
from pyeosim.datasets import names, DATA_PATHS
from Py6S import SixS
import xarray

# load, set wavelength as dim coords and rescale to [0..1]
reflectance = xarray.load_dataarray(DATA_PATHS['TEST_HSI_LARGE'])
reflectance = reflectance.swap_dims({'band':'wavelength'}) / 10000

# generate a 6SV atmsophere
atm = SixSV_atmosphere(SixS(), TreeView_3())

# Generate a default TDI CMOS sensor
TreeView_satellite = TdiCmos()

# Perform simulation
imager_output = TreeView_satellite.fit_transform(atm.transform(reflectance))

# plot outputs
imager_output.plot(col='band', col_wrap=3, cmap='gray')
```
## API Reference
| Submodules |
| --- |
[atmosphere](docs/pyeosim.atmosphere.md)
[datasets](docs/pyeosim.datasets.md)
[imager](docs/pyeosim.imager.md)
[plot](docs/pyeosim.plot.md)
[spatial](docs/pyeosim.spatial.md)
[spectral](docs/pyeosim.spectral.md)
