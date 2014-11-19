#!/bin/env python

import os, json
from sys import argv
from glob import glob
from time import sleep

if len(argv) < 2:
    print "usage: " + argv[0] + " <config file> [outdir]"
    exit(0)

configfile = argv[1]
if os.path.exists(configfile) is False:
    print "configfile " + configfile + " does not exist!"
    exit(0)

print "opening configfile " + configfile
fp = open(configfile)
config = json.load(fp)
fp.close()

shapedir = config['shape_dir']
lensfile = config['lens_catalog']

if len(argv) > 2:
    outdir = argv[2]
else:
    outdir = os.path.dirname(configfile)
if outdir[-1] != '/':
    outdir += '/'

shapefiles = glob(shapedir + "/" + "*.fits.gz")
thisdir = os.path.dirname(os.path.realpath(__file__))
for shapefile in shapefiles:
    basename = os.path.basename(shapefile)
    basename = basename.split(".")[0]
    stackfile = outdir + basename + '_DeltaSigma.fits'
    lockfile = outdir + basename + ".lock"
    if os.path.exists(stackfile) is False and os.path.exists(lockfile) is False:
        print basename
        os.system('hostname > ' + lockfile)
        os.popen('python ' + thisdir +'/stack_slices.py ' + lensfile + " " + shapefile + " " + outdir + " 1>&2 >> " + lockfile + " && rm " + lockfile)

