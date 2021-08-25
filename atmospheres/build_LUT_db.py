import concurrent.futures
import json
import logging
import numpy as np
import os
import pandas
from Py6S import *


# reflectance values
rho_range = np.linspace(0.01, 1, 100)
# rho_range = np.linspace(0.01, 1, 3)
# wavelengths [0.4...0.9] in microns
wlens = np.linspace(.4, .9, 501)

def run_lambertian(rho, lat, lon, dt_string, view_z, view_a, save_path, output_json=False):
    # Create a SixS object called s (used as the standard name by convention)
    s = SixS()
    s.ground_reflectance = GroundReflectance.HomogeneousLambertian(rho)
    geom = Geometry.User()
    geom.from_time_and_location(lat=lat,
                                lon=lon,
                                datetimestring=dt_string,
                                view_z=0,
                                view_a=0)
    s.geometry = geom
    s.altitudes.set_sensor_satellite_level()
    # Run the 6S simulation defined by this SixS object across the
    # whole VNIR range
    _lambda, _res = SixSHelpers.Wavelengths.run_wavelengths(s, wlens,
                                                            output_name='pixel_radiance')
    fname = os.path.join(save_path, '_{:.0f}.csv'.format(rho*100))
    # result to pandas and save
    pandas.DataFrame({'lambda': _lambda,
                      'radiance': _res}).set_index('lambda').to_csv(fname)
    if output_json:
        # run just once and save output to file
        s.run()
        with open('{}.json'.format(save_path), 'w') as f:
            f.write(json.dumps(s.outputs.values))
    return fname

def main():
    # master save path
    main_path = './atmospheres_2/'
    # add new inputs here but be careful to change folder name or will
    # overwrite
    inputs =[
        ('2020/22/06T10:30', 57.47, 4.22, 'inv_20200622_1030_nadir'),
        ('2020/22/06T12:00', 57.47, 4.22, 'inv_20200622_1200_nadir'),
        ('2020/01/03T10:30', 57.47, 4.22, 'inv_20200301_1030_nadir'),
        ('2020/01/10T10:30', 57.47, 4.22, 'inv_20201001_1030_nadir'),
        ('2020/01/08T10:30', 57.47, 4.22, 'inv_20200801_1030_nadir'),
        ('2020/01/05T10:30', 57.47, 4.22, 'inv_20200501_1030_nadir'),
        ('2020/22/06T10:30', 52.04, 0.76, 'mk_20200622_1030_nadir'),
        ('2020/22/06T12:00', 52.04, 0.76, 'mk_20200622_1200_nadir'),
        ('2020/01/03T10:30', 52.04, 0.76, 'mk_20200301_1030_nadir'),
        ('2020/01/10T10:30', 52.04, 0.76, 'mk_20201001_1030_nadir'),
        ('2020/01/08T10:30', 52.04, 0.76, 'mk_20200801_1030_nadir'),
        ('2020/01/05T10:30', 52.04, 0.76, 'mk_20200501_1030_nadir'),
    ]
    # iterate
    for params in inputs:
        dt_string = params[0]
        lat = params[1]
        lon = params[2]
        save_path = os.path.join(main_path, params[3])
        os.makedirs(save_path, exist_ok=True)
        for _rho in rho_range:
            if _rho == rho_range[-1]:
                save = True
            else:
                save = False
            _fname = run_lambertian(_rho, lat, lon, dt_string, 0, 0, save_path, save)
            logging.info('{} saved'.format(_fname))

if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    main()