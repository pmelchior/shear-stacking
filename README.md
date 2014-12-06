shear-stacking-tests
====================

This repository contains python scripts to calculate stacked shear profiles and tests based upon them, e.g. consistency for different slices of lensed background galaxies. The basic concept is that the lensing signal in terms of surface mass density (instead of shear) should be entirely determined by the properties of the lens sample and have no dependence on source galaxy properties. An example application is given in section 4.4 of [this paper](http://arxiv.org/abs/1405.4285).

Usage
-----

```
python stack_slices.py test/config.json
python plot_slices.py test/config.json
```

The first command opens the config file, creates FITS tables for each of the shape catalogs listed in the config file, and stores them in the directory `test/`. The second command, uses the output of the first and generates stacked shear profiles (in either angular or physical units) for each of the source slices and stores them in the directory `test/`. After that, all `*.fits` files in `test/` can be deleted, or the script can be rerun with modified binning/colors...

```
python stack_slices.py test/config.json shape_catalog.fits
```
Overrides the default shape catalog selection in the config file and only works on `shape_catalog.fits`.

Configuration file format
-------------------------

Virtually all aspects of the script can be controlled from a config file, in json format. This way, the scripts do not have to be altered to adjust to the pecularities of the lens or shape catalogs.

One example is given below:

```json
{
        "coords": "angular",
        "maxrange": 1.1,
        "lens_catalog": "/catalogs/redmapper/redmapper_catalog.fit",
        "lens_cuts": [],
        "lens_z_key": "Z_LAMBDA",
        "shape_dir": "/catalogs/im3shapev7_ngmix009/",
        "shape_z_key": "ZP",
        "shape_ra_key": "ALPHAWIN_J2000_R",
        "shape_dec_key": "DELTAWIN_J2000_R",
        "shape_e1_key": "im3shape_r_e1",
        "shape_e2_key": "im3shape_r_e2",
        "shape_weight_key": "weight_im3shape",
        "shape_sensitivity_key": "nbc_m"
        "shape_cuts": [
                "MODEST_CLASS == 1",
                "im3shape_r_exists == 1",
                "im3shape_r_error_flag == 0"
        ],
        "splittings": {
                "FLAGS_I": [0, 1, 2, 4],
                "ZP": [0.7, 0.9, 1.1, 1.5],
                "B_D": [0.0, 0.3, 0.7, 1.0]
        },
        "functions": {
                "weight_im3shape": "0.2**2/(0.2**2 + (0.1*20/s['im3shape_r_snr'])**2)",
                "nbc_m": "1 + s['im3shape_r_nbc_m']",
                "B_D": "s['im3shape_r_bulge_flux'] / s['im3shape_r_disc_flux']"
        }
}
```

Coordinates can be either `angular` or `physical` and refer to the units of the radial profiles. `maxrange` is the maximum extent of the profile, in deg or Mpc/h (depending on `coords`).

`lens_cuts` and `shape_cuts` affect which objects are loaded from either lensing or shape catalogs; filtering can be used on all columns available in either fits files.

`splitting` denote the kind of slices of the shape catalogs. The keys of this dictionary specify either a column in the shape catalog or an entry in the `functions` list. For the latter, any function based on the shape catalog (denoted as `s`) can be implemented. The values denote the limits of the slices for each key, with the upper limit being excluded, e.g. `"FLAGS_I": [0, 1, 2, 4]` creates three slices:
```
FLAGS_I in [0,1), [1,2), [2,4) 
```

There is no formal limit on how many different categories/keys can be done, but the number of slices should not exceed 5 (otherwise the plots get rather busy).

User-defined `functions` can be used for any `_key` entry in the config file, both for lens and for shape catalogs.

Dependencies
------------

* [GalSim](https://github.com/GalSim-developers/GalSim)
* [Erin Sheldon's esutil](https://code.google.com/p/esutil/)
* [Erin Sheldon's fitsio](https://github.com/esheldon/fitsio)
* matplotlib
