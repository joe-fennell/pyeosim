<!-- markdownlint-disable -->

<a href="../pyeosim/atmosphere.py#L0"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

# <kbd>module</kbd> `pyeosim.atmosphere`
Atmospheric Simulation 

This module contains atmospheric simulation transformer classes. 

**Global Variables**
---------------
- **DATA_PATHS**

---

<a href="../pyeosim/atmosphere.py#L144"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../pyeosim/atmosphere.py#L15"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `SixSV_atmosphere`
Generates an atmospheric transformer for a specific 6SV atmosphere  



<a href="../pyeosim/atmosphere.py#L18"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__(SixS, srf=None)
```



**Args:**
 
 - <b>`SixS`</b>:  A SixS atmosphere instance 




---

<a href="../pyeosim/atmosphere.py#L33"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `fit`

```python
fit(srf)
```

Runs 6SV and generates the lookup coefficients 



**Args:**
 
 - <b>`srf`</b>:  spectral response function object 

---

<a href="../pyeosim/atmosphere.py#L106"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../pyeosim/atmosphere.py#L87"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

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

<a href="../pyeosim/atmosphere.py#L122"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

## <kbd>class</kbd> `Test6S`
Look-up table test for old-style lookup tables  



<a href="../pyeosim/atmosphere.py#L126"><img align="right" style="float:right;" src="https://img.shields.io/badge/-source-cccccc?style=flat-square"></a>

### <kbd>method</kbd> `__init__`

```python
__init__()
```











---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
