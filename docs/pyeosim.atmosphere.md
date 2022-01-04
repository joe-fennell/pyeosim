<!-- markdownlint-disable -->

# <kbd>module</kbd> `pyeosim.atmosphere`
Atmospheric Simulation 

**Global Variables**
---------------
- **DATA_PATHS**

---

## <kbd>function</kbd> `LUT_from_file`

```python
LUT_from_file(fpath, common_params={})
```

Generates LUT transformer from file. 

Takes a directory of directories (one per scenario) 



**Note:**

> Each file in the directory should be a CSV named with the rho (0...100) 
>

**Args:**
 
 - <b>`fpath`</b>:  filepath to a valid CSV response file 
 - <b>`common_params`</b> (dict):  attributes to add to DataArray 



**Returns:**
 LUT transformer instance 


---

## <kbd>class</kbd> `SixSV_atmosphere`
Generates an atmospheric transformer for a specific 6SV atmosphere  



### <kbd>method</kbd> `__init__`

```python
__init__(SixS, srf=None)
```



**Args:**
 
 - <b>`SixS`</b>:  A SixS atmosphere instance 




---

### <kbd>method</kbd> `fit`

```python
fit(srf)
```

Runs 6SV and generates the lookup coefficients 



**Args:**
 
 - <b>`srf`</b>:  spectral response function object 

---

### <kbd>method</kbd> `inverse_transform`

```python
inverse_transform(signal)
```

Convert TOA radiance to BOA reflectance 



**Args:**
 
 - <b>`signal`</b>:  per band radiance dataset 



**Returns:**
 per band BOA reflectance dataset 

---

### <kbd>method</kbd> `transform`

```python
transform(signal)
```

Convert BOA reflectance to per-band TOA radiance 



**Args:**
 
 - <b>`signal`</b>:  Spectral reflectance dataset with 'wavelength' dimension 



**Returns:**
 TOA radiance array with 'bands' spectral dimension 


---

## <kbd>class</kbd> `Test6S`
Look-up table test for old-style lookup tables  



### <kbd>method</kbd> `__init__`

```python
__init__()
```











---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
