<!-- markdownlint-disable -->

# <kbd>module</kbd> `pyeosim.post`
Additional post processing tools for image data 


---

## <kbd>function</kbd> `apply_downsampling`

```python
apply_downsampling(arr, sensor, spatial=True, spectral=False, normalise=False)
```

Apply only the spectral and spatial downsampling steps of the sensor. 



**Args:**
 
 - <b>`arr`</b>:  array of radiance or reflectance values 
 - <b>`sensor`</b>:  a fitted instance of a TdiCmos 
 - <b>`spatial`</b> (bool):  applies spatial downsampling if True 
 - <b>`spectral`</b> (bool):  applies spectral downsampling if True 
 - <b>`normalise`</b> (bool):  Normalises each band by the integral of the spectral  response function - use if input arr is reflectance 



**Returns:**
 New array of values 


---

## <kbd>function</kbd> `sensor_correction_experiment`

```python
sensor_correction_experiment(toa_radiance, sensor, mask=None)
```

Performs a calibration experiment based on an input radiance dataset and a sensor to give the correction factors needed to convert DN to apparent radiance. Note that this assumes a linear relationship between DN and radiance. 



**Args:**
 
 - <b>`toa_radiance`</b>:  At-sensor (top-of-atmosphere) Radiance 
 - <b>`sensor `</b>:  A valid sensor instance mask (bool): 



**Returns:**
 per-band calibration coefficients (m,c) for linear model to convert DN to radiance 


---

## <kbd>class</kbd> `LinearRadiometricCorrection`
Linear radiometric correction using precomputed coefficients  



### <kbd>method</kbd> `__init__`

```python
__init__(coef_array=None, coef_filepath=None)
```



**Args:**
 
 - <b>`coef_array`</b> (array-like):  xarray.Dataset with 'm' and 'c' variables 
 - <b>`coef_filepath`</b> (str):  path to NetCDF Dataset with 'm' and 'c'  variables 




---

### <kbd>method</kbd> `transform`

```python
transform(signal)
```

Applies linear correction 



**Args:**
 
 - <b>`signal`</b>:  Digitial Number array with bands matching those in  correction coefficient array 



**Returns:**
 Radiance values 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
