shear-stacking-tests
====================

This repository contains python scripts to calculate stacked shear profiles and tests based upon them, e.g. consistency for different slices of lensed background galaxies. An example application is given in section 4.4 of [this paper](http://arxiv.org/abs/1405.4285).

Configuration file format
-------------------------

Virtually all aspects of the script can be controlled from a config file, in json format. One example is given below:

```json
{
        "coords": "angular",
        "maxrange": 1.1,
        "lens_catalog": "/catalogs/redmapper/sva1_gold_1.0.2_redmapper_v6.2.12_lgt5_desformat_catalog.fit",
        "lens_cuts": [],
        "lens_z_key": "Z_LAMBDA",
        "shape_dir": "/catalogs/im3shapev7_ngmix009/",
        "shape_z_key": "ZP",
        "shape_ra_key": "ALPHAWIN_J2000_R",
        "shape_dec_key": "DELTAWIN_J2000_R",
        "shape_e1_key": "im3shape_r_e1",
        "shape_e2_key": "im3shape_r_e2",
        "shape_weight_key": "SNR_weight",
        "shape_cuts": [
                "MODEST_CLASS == 1",
                "im3shape_r_exists == 1",
                "im3shape_r_error_flag == 0"
        ],
        "splittings": {
                "FLAGS_I": [0, 1, 2, 4],
	        "ZP": [0.7, 0.9, 1.1, 1.5]
        }
}
```

Coordinates can be either `angular` or `physical` and refer to the units of the radial profiles. `maxrange` is the maximum extent of the profile, in deg or Mpc/h (depending on `coords`).

`lens_cuts` and `shape_cuts` affect which objects are loaded from either lensing or shape catalogs; filtering can be used on all columns available in either fits files.

`splitting` denote the kind of slices of the shape catalogs. The keys of this dictionary specify either a column in the shape catalog or a function in the global scope of the scripts. For the latter, it is most convenient to add this function to `common.py` (The same can be done for `shape_weight_key`). The values denote the limits of the slices for each key, with the upper limit being excluded, e.g. `"FLAGS_I": [0, 1, 2, 4]` creates three slices:
```
FLAGS_I in [0,1), [1,2), [2,4) 
```

There is no formal limit on how many different categories/keys can be done, but the number of slices should not exceed 5 (otherwise the plots get rather busy). 

Dependencies
------------

* [GalSim](https://github.com/GalSim-developers/GalSim)
* [Erin Sheldon's esutil](https://code.google.com/p/esutil/)
* [Erin Sheldon's fitsio](https://github.com/esheldon/fitsio)
* matplotlib
