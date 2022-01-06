<!-- markdownlint-disable -->

# <kbd>module</kbd> `pyeosim.datasets`
Dataset Management 

Datasets and data load funcs used in testing and specific supported sensors 

The DATA_PATHS dictionary currently contains the following filepaths: 



**Attributes:**
 
 - <b>`TEST_HSI`</b>:  A Vis-VNIR NetCDF4 hyperspectral dataset of 186 spectral bands  and 36 pixels (186x6x6). 
 - <b>`TEST_HSI_LARGE`</b>:  A Vis-VNIR NetCDF4 hyperspectral dataset of 186 spectral  bands and 131,042 pixels (186x362x362). 
 - <b>`SRF_SENTINEL_2`</b>:  per-band Spectral Response Function for the Sentinel 2A and  2B satellites. 
 - <b>`SRF_SUPERDOVE`</b>:  per-band Spectral Response Function for the Planet SuperDove  satellites. 
 - <b>`SOLAR_SPECTRUM_ASTMG173`</b>:  Mean solar surface irradiance spectrum (global).  Wavelength in nm and irradiance in W m-2 nm-1. 
 - <b>`SOLAR_SPECTRUM_ASTME490`</b>:  Mean solar extraterrestrial irradiance spectrum  (global). Wavelength in nm and irradiance in W m-2 nm-1. 
 - <b>`CCD_QE_DD_BACK`</b>:  Published Teledyne sensor Quantum Efficiency for a Dump  Drain back-thinned sensor. 
 - <b>`CCD_QE_STD_BACK`</b>:  Published Teledyne sensor Quantum Efficiency for a  standard back-thinned sensor. 
 - <b>`TDI_QE_BACK`</b>:  Approximate Quantum Efficiency for a TDI back-  thinned sensor. 
 - <b>`TEST_LUT`</b>:  A Top-Of-Atmosphere radiance lookup table for reflectances  at wavelengths between 400 and 900nm. 

**Global Variables**
---------------
- **DATA_PATHS**

---

## <kbd>function</kbd> `names`

```python
names()
```

List all datasets  






---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
