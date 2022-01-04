<!-- markdownlint-disable -->

# <kbd>module</kbd> `pyeosim.spatial`
Spatial Response Generators 


---

## <kbd>function</kbd> `gaussian_isotropic`

```python
gaussian_isotropic(signal, psf_fwhm, ground_sample_distance)
```

Simulates a gaussian optic and sensor plane. 

Apply downsampling transform assuming a gaussian PSF and known PSF full width at half maximum. 



**Args:**
 
 - <b>`signal`</b>:  a 2D xarray raster with valid 'res' in attributes 
 - <b>`psf_fwhm`</b> (float):  Full width at half maximum of gaussian kernel. 
 - <b>`ground_sample_distance`</b> (float):  ground resolution of instrument 



**Returns:**
 2D xarray raster array at new resolution 




---

_This file was automatically generated via [lazydocs](https://github.com/ml-tooling/lazydocs)._
