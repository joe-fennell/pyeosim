<!-- markdownlint-disable -->

# API Overview

## Modules

- [`pyeosim`](./pyeosim.md#module-pyeosim)
- [`pyeosim.atmosphere`](./pyeosim.atmosphere.md#module-pyeosimatmosphere): Atmospheric Simulation
- [`pyeosim.datasets`](./pyeosim.datasets.md#module-pyeosimdatasets): Dataset Management
- [`pyeosim.imager`](./pyeosim.imager.md#module-pyeosimimager): Imager and System Simulation
- [`pyeosim.plot`](./pyeosim.plot.md#module-pyeosimplot): Plotting Utilities
- [`pyeosim.spatial`](./pyeosim.spatial.md#module-pyeosimspatial): Spatial Response Generators
- [`pyeosim.spectral`](./pyeosim.spectral.md#module-pyeosimspectral): Multispectral Spectral Responses and Generators

## Classes

- [`atmosphere.SixSV_atmosphere`](./pyeosim.atmosphere.md#class-sixsv_atmosphere): Generates an atmospheric transformer for a specific 6SV atmosphere
- [`atmosphere.Test6S`](./pyeosim.atmosphere.md#class-test6s): Look-up table test for old-style lookup tables
- [`imager.TdiCmos`](./pyeosim.imager.md#class-tdicmos): Time Delay Integration CMOS generic transformer.
- [`spectral.Sentinel2A`](./pyeosim.spectral.md#class-sentinel2a): Sentinel 2A - all spectral bands.
- [`spectral.Sentinel2B`](./pyeosim.spectral.md#class-sentinel2b): Sentinel 2B - all bands.
- [`spectral.Sentinel2VNIR`](./pyeosim.spectral.md#class-sentinel2vnir): Sentinel 2A - vis-VNIR only.
- [`spectral.SuperDove`](./pyeosim.spectral.md#class-superdove): SuperDove
- [`spectral.TreeView_1`](./pyeosim.spectral.md#class-treeview_1): TreeView Version 1
- [`spectral.TreeView_2`](./pyeosim.spectral.md#class-treeview_2): TreeView Version 2
- [`spectral.TreeView_3`](./pyeosim.spectral.md#class-treeview_3): TreeView Version 3

## Functions

- [`atmosphere.LUT_from_file`](./pyeosim.atmosphere.md#function-lut_from_file): Generates LUT transformer from file.
- [`datasets.names`](./pyeosim.datasets.md#function-names): List all datasets
- [`imager.CONU`](./pyeosim.imager.md#function-conu): Generates column offset non uniformity array
- [`imager.DSNU`](./pyeosim.imager.md#function-dsnu): Generates the dark signal non-uniformity array.
- [`imager.PRNU`](./pyeosim.imager.md#function-prnu): Generates the Photo Response Non-Uniformity.
- [`imager.add_column_offset`](./pyeosim.imager.md#function-add_column_offset): Adds the column offset noise to voltage array.
- [`imager.add_dark_signal`](./pyeosim.imager.md#function-add_dark_signal): Adds dark current to electron count.
- [`imager.add_gaussian_noise`](./pyeosim.imager.md#function-add_gaussian_noise): Adds a zero-mean gaussian to electron count
- [`imager.add_photon_noise`](./pyeosim.imager.md#function-add_photon_noise): Adds photon shot noise.
- [`imager.add_prnu`](./pyeosim.imager.md#function-add_prnu): Adds Photon Response Non-Uniformity to electrons.
- [`imager.electron_to_voltage`](./pyeosim.imager.md#function-electron_to_voltage): Converts electrons (charge) to voltage.
- [`imager.energy_to_quantity`](./pyeosim.imager.md#function-energy_to_quantity): Convert energy to quanta.
- [`imager.photon_mean`](./pyeosim.imager.md#function-photon_mean): Convert sensor photon flux density to photon count at sensor.
- [`imager.photon_to_electron`](./pyeosim.imager.md#function-photon_to_electron): Convert photon count to electron count.
- [`imager.radiance_to_irradiance_2`](./pyeosim.imager.md#function-radiance_to_irradiance_2): Convert radiance (W m-2 sr-1) to irradiance of sensor based on lens.
- [`imager.voltage_to_DN`](./pyeosim.imager.md#function-voltage_to_dn): Convert voltage to Digital Number via a linear ADC.
- [`plot.rgb`](./pyeosim.plot.md#function-rgb): Make an RGB image from a 3 channel array.
- [`spatial.gaussian_isotropic`](./pyeosim.spatial.md#function-gaussian_isotropic): Simulates a gaussian optic and sensor plane.


---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
