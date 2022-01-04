<!-- markdownlint-disable -->

# <kbd>module</kbd> `pyeosim.imager`
Imager and System Simulation 


---

## <kbd>function</kbd> `add_column_offset`

```python
add_column_offset(voltage, offset)
```

Adds the column offset noise to voltage array. 



**Args:**
 
 - <b>`voltage`</b>:  voltage array 
 - <b>`offset`</b>:  column voltage offset array 



**Returns:**
 voltage array 


---

## <kbd>function</kbd> `add_dark_signal`

```python
add_dark_signal(electrons, dark_current, integration_time, dsnu)
```

Adds dark current to electron count. 



**Args:**
 
 - <b>`electrons`</b>:  electron array 
 - <b>`dark_current`</b> (float):  dark current in e-/s/pixel 
 - <b>`integration_time`</b> (float):  integration period in seconds 
 - <b>`dsnu `</b>:  dark current fixed pattern noise array 



**Returns:**
 electron count array 


---

## <kbd>function</kbd> `add_gaussian_noise`

```python
add_gaussian_noise(electron_count, sigma)
```

Adds a zero-mean gaussian to electron count 



**Args:**
 
 - <b>`electron_count`</b>:  electron count array 
 - <b>`sigma`</b>:  SD of noise in electrons 



**Returns:**
 electron count array 


---

## <kbd>function</kbd> `add_photon_noise`

```python
add_photon_noise(photon_count)
```

Adds photon shot noise. 

Assumes signal is a Poisson random process. 



**Args:**
 
 - <b>`photon_count`</b>:  mean photon count array 



**Returns:**
 photon count array 


---

## <kbd>function</kbd> `add_prnu`

```python
add_prnu(electrons, prnu)
```

Adds Photon Response Non-Uniformity to electrons. 



**Args:**
 
 - <b>`electrons `</b>:  electron count array 
 - <b>`prnu `</b>:  fixed pattern noise array 



**Returns:**
 electron count array 


---

## <kbd>function</kbd> `CONU`

```python
CONU(ones, offset_factor)
```

Generates column offset non uniformity array 



**Args:**
 
 - <b>`ones`</b>:  array of ones of same shape as image 
 - <b>`offset_factor`</b> (float):  column offset factor 



**Returns:**
 Column offset non uniformity array with zero mean 


---

## <kbd>function</kbd> `DSNU`

```python
DSNU(ones, dark_current, integration_time, dark_factor)
```

Generates the dark signal non-uniformity array. 

This assumes a log-Normal distribution suitable for short integration times of less than 100 seconds 



**Note:**

> dark_current should be per actual pixel, not per pixel subarray 
>

**Args:**
 
 - <b>`ones`</b>:  array of ones in same shape as image 
 - <b>`dark_current`</b> (float):  sensor level mean dark current in e-/s/pixel 
 - <b>`integration_time`</b> (float):  integration time of exposure 
 - <b>`dark_factor`</b> (float):  dark current FPN factor range (0...1) 



**Returns:**
 electron count array 


---

## <kbd>function</kbd> `electron_to_voltage`

```python
electron_to_voltage(
    electron_count,
    v_ref,
    sense_node_gain,
    full_well,
    read_noise
)
```

Converts electrons (charge) to voltage. 

Converts electrons (charge) to voltage by simulating a sense node with a fixed v_ref 



**Args:**
 
 - <b>`electron_count`</b>:  electron count array 
 - <b>`v_ref`</b> (float):  sense node reference voltage 
 - <b>`sense_node_gain`</b> (float):  gain in V/e- 
 - <b>`full_well`</b> (int):  max number of electrons per pixel 
 - <b>`read_noise`</b> (int):  read noise in electrons per channel 



**Returns:**
 sensor voltage array 


---

## <kbd>function</kbd> `energy_to_quantity`

```python
energy_to_quantity(energy)
```

Convert energy to quanta. 



**Args:**
 
 - <b>`energy`</b>:  energy in J 



**Returns:**
 quantity or array of quantities 


---

## <kbd>function</kbd> `photon_mean`

```python
photon_mean(flux, pixel_area, integration_time)
```

Convert sensor photon flux density to photon count at sensor. 



**Args:**
 
 - <b>`flux`</b>:  pixel flux density (photons/second/m2) 
 - <b>`pixel_area`</b> (float):  area in microns^2 
 - <b>`integration_time`</b> (float):  integration time in seconds 



**Returns:**
 photon count array 


---

## <kbd>function</kbd> `photon_to_electron`

```python
photon_to_electron(photon_count, Q_E)
```

Convert photon count to electron count. 

Convert photon count to electron count using sensor quantum efficiency. 



**Args:**
 
 - <b>`photon_count`</b>:  photon count array 
 - <b>`Q_E`</b>:  float in range 0-1 specifying photon conversion rate or xarray  with same wavelength sampling as photon_count 

**Returns:**
 electron count array 


---

## <kbd>function</kbd> `PRNU`

```python
PRNU(ones, prnu_factor)
```

Generates the Photo Response Non-Uniformity. 



**Args:**
 
 - <b>`ones`</b>:  array of ones in same shape as image 
 - <b>`prnu_factor`</b>:  PRNU in range 0...1 



**Returns:**
 electron count array 


---

## <kbd>function</kbd> `radiance_to_irradiance_2`

```python
radiance_to_irradiance_2(radiance, lens_diameter, focal_length)
```

Convert radiance (W m-2 sr-1) to irradiance of sensor based on lens. 



**Note:**

> Only an approximation for nadir sensor. Forshortening is not taken into account. Does not take into account cos^4 roll-off 
>

**Args:**
 
 - <b>`radiance`</b>:  at sensor radiance array 
 - <b>`lens_diameter`</b> (float):  diameter of optic in metres 
 - <b>`focal_length`</b> (float):  focal length of optic in metres 



**Returns:**
 sensor irradiance array 


---

## <kbd>function</kbd> `voltage_to_DN`

```python
voltage_to_DN(voltage, v_ref, bit_depth)
```

Convert voltage to Digital Number via a linear ADC. 

Model a linear response ADC. 



**Args:**
 
 - <b>`voltage`</b>:  voltage array 
 - <b>`v_ref`</b> (float):  reference voltage of ADC 
 - <b>`bit_depth`</b> (int):  ADC bit depth 



**Returns:**
 Digital Number array 


---

## <kbd>class</kbd> `TdiCmos`
Time Delay Integration CMOS generic transformer. 



**Note:**

> This version does not apply spectral resampling and relies on band integration happening in the atmospheric simulation 

### <kbd>method</kbd> `__init__`

```python
__init__(
    sensor_altitude=500000.0,
    sensor_ground_speed=7000,
    ground_sample_distance=2,
    lens_diameter=0.1,
    psf_fwhm=4,
    TDI_rows=32,
    pix_per_row=8000,
    sensor_width=82.2,
    pixel_area=100,
    spectral_response='TreeView_3()',
    quantum_efficiency='TDI_QE_BACK',
    full_well=30000.0,
    prnu_factor=0.01,
    dark_current=570,
    dark_factor=0.01,
    offset_factor=0.001,
    ccd_vref=3.1,
    sense_node_gain=5,
    read_noise=20,
    adc_vref=0.5,
    bit_depth=14,
    store_steps=False,
    apply_spatial_resampling=True
)
```



**Args:**
 
 - <b>`sensor_altitude`</b> (float):  altitude above source in metres 
 - <b>`sensor_ground_speed`</b> (float):  relative ground speed of sensor in m/s 
 - <b>`ground_sample_distance`</b> (float):  Ground sampling distance of sensor 
 - <b>`lens_diameter`</b> (float):  diameter of lens in metres 
 - <b>`psf_fwhm`</b> (float):  Point Spread Function Full-Width at Half Maximum  in ground units (metres) 
 - <b>`TDI_rows`</b> (int):  number of rows integrated over per channel 
 - <b>`pix_per_row`</b> (int):  number of pixels per row 
 - <b>`sensor_width`</b> (float):  width of sensor imaging area in mm 
 - <b>`pixel_area`</b> (float):  per pixel light-absorbing area in micron2 
 - <b>`spectral_response`</b>:  spectral transformer instance 
 - <b>`quantum_efficiency`</b>:  Quantum efficiency of sensor as a float,  name of internal dataset to use or xarray.DataArray with a  'wavelength' dimension covering the wavelength range of the  spectral_response or 'bands' dimension matching the bands in  spectral_response 
 - <b>`full_well`</b> (int):  the maximum capacity of a pixel in electrons 
 - <b>`prnu_factor`</b> (float):  Photo response non uniformity factor 
 - <b>`dark_current`</b> (float):  dark current in electrons/second/pixel 
 - <b>`dark_factor`</b> (float):  dark current fixed pattern noise factor  (DSNU factor) 
 - <b>`offset_factor`</b> (float):  column offset factor 
 - <b>`ccd_vref`</b> (float):  reference voltage of voltage sensor in volts 
 - <b>`sense_node_gain`</b> (float):  gain of sense node in microvolts/electron 
 - <b>`read_noise`</b> (int):  Read noise in e- 
 - <b>`adc_vref`</b> (float):  reference voltage of ADC in volts 
 - <b>`bit_depth`</b> (int):  bit depth of the ADC 
 - <b>`store_steps`</b> (float):  if True, all intermediate steps will be stored in the step_outputs attribute. Default is False 
 - <b>`apply_spatial_resampling`</b> (float):  if True, spatial response  function will be applied. Use False is supplying non-imaging  data to the simulation 




---

### <kbd>method</kbd> `fit`

```python
fit(signal)
```

Precomputes the system constant states. 

This precomputes states for system constants (e.g. the fixed- pattern noise). 



**Args:**
 
 - <b>`signal`</b>:  A Top-Of-Atmosphere radiance dataset integrated over the  bandpass response 

---

### <kbd>method</kbd> `fit_transform`

```python
fit_transform(signal)
```

Runs simulation on signal Top-Of-Atmosphere radiance dataset. 

This both precomputes the system state constants (e.g. the fixed- pattern noise) and generates a temporary random state for other parameters (e.g. the dark signal). 



**Args:**
 
 - <b>`signal`</b>:  A Top-Of-Atmosphere radiance dataset integrated over the  bandpass response 



**Returns:**
 sensor digital number output 

---

### <kbd>method</kbd> `transform`

```python
transform(signal)
```

Runs simulation on signal Top-Of-Atmosphere radiance dataset. 

This uses the precomputed states for system constants (e.g. the fixed- pattern noise) and generates a temporary random state for other parameters (e.g. the dark signal). 



**Args:**
 
 - <b>`signal`</b>:  A Top-Of-Atmosphere radiance dataset integrated over the  bandpass response 



**Returns:**
 sensor digital number output 

---

### <kbd>method</kbd> `update_derived_params`

```python
update_derived_params()
```

Re-run all derived parameter calculations. 

Call after updating any parameter values to recalculate derived params 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
