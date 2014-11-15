#!/bin/env python

import numpy as np
import esutil as eu
import os
from sys import argv
from common import *
import math 
import fitsio
import json

def getValues(data, key):
    # what values are used for the slices
    if key in globals():
        return eval(key)(data)
    else:
        return data[key]

def getSliceMask(values, lower, upper, return_num=False):
    if return_num is False:
        return (values >= lower) & (values < upper)
    else:
        return sum((values >= lower) & (values < upper))

# get separation in deg for distance L in Mpc/h at redshift z
# uses c/H0 = 3000 Mpc/h
def Dist2Ang(L, z):
    global cosmo
    return L / cosmo.Da(z) / 3000. * 180./math.pi

def Ang2Dist(theta, z):
    global cosmo
    return theta * cosmo.Da(z) * 3000. / 180. * math.pi

if __name__ == '__main__':
    if len(argv) < 4:
        print "usage: " + argv[0] + " <lens catalog> <shape catalog> <output dir>"
        exit(1)

    lensfile = argv[1]
    shapefile = argv[2]
    outdir = argv[3]
    os.system('mkdir -p ' + outdir)
    if outdir[-1] != '/':
        outdir += '/'
    
    basename = os.path.basename(shapefile)
    basename = basename.split(".")[0]
    matchfile = outdir + basename + '_matches.bin'
    stackfile = outdir + basename + '_DeltaSigma.fits'
    configfile = outdir + 'config.json'

    if os.path.exists(configfile) is False:
        print "configfile " + configfile + " does not exist!"
        exit(0)
    else:
        print "opening configfile " + configfile
        fp = open(configfile)
        config = json.load(fp)
        fp.close()

    if config['coords'] == "physical":
        maxrange = 5. # Mpc/h
    else:
        maxrange = 1.1  # deg

    if os.path.exists(stackfile) is False:
        # open lens catalog, apply selection if desired
        hdu = fitsio.FITS(lensfile)
        if len(config['lens_cuts']) == 0:
            lenses = hdu[1][:]
        else:
            cuts = " && ".join(config['lens_cuts'])
            mask = hdu[1].where(cuts)
            lenses = hdu[1][mask]
            del mask
        print "lens sample: %d" % lenses.size
        
        if config['coords'] == "physical":
            maxrange = Dist2Ang(maxrange, lenses[lens_z_key])

        # open shapes, apply post-run selections
        shdu = fitsio.FITS(shapefile)
        if len(config['shape_cuts']) == 0:
            shapes = hdu[1][:]
        else:
            cuts = " && ".join(config['shape_cuts'])
            mask = shdu[1].where(cuts)
            shapes = shdu[1][mask]
            del mask
        print "shape sample: %d" % shapes.size

        # find all galaxies in shape catalog within maxrange arcmin 
        # of each lens center
        print "matching lens and source catalog..."
        if os.path.exists(matchfile) is False:
            # CAVEAT: make sure to have enough space where you put the match file
            # it has 24 byte per match, which quickly becomes Gb's of data 
            h = eu.htm.HTM(8)
            h.match(lenses['RA'], lenses['DEC'], shapes[config['shape_ra_key']], shapes[config['shape_dec_key']], maxrange, maxmatch=-1, file=matchfile)
            del h
        else:
            print "  re-using existing matchfile", matchfile

        htmf = HTMFile(matchfile)
        Ngal = htmf.n_matches
        print "  found ", Ngal, "matches"

        # iterate over all lenses, write DeltaSigma, r, weight into file
        print "stacking lenses..."
        fits = fitsio.FITS(stackfile, 'rw')
        specz_calib = getSpecZCalibration()
        counter = 0
        done = 0
        for m1, m2, d12 in htmf.matches():
            lens = lenses[m1]
            shapes_lens = shapes[m2]
            n_gal = shapes_lens.size

            # compute effective Sigma_crit
            z_phot, cz = getSigmaCritCorrection(specz_calib, lens[config['lens_z_key']])
            sigma_crit = getSigmaCritEffective(z_phot, cz, shapes_lens[config['shape_z_key']])
            # determine extent in DeltaSigma array
            data = np.empty(shapes_lens.size, dtype=[('radius', 'f8'), ('DeltaSigma', 'f8'), ('DeltaSigma_x', 'f8'), ('weight', 'f8'), ('slices', 'S%d' % len(config['splittings']))])
            data['DeltaSigma'], data['DeltaSigma_x'] = sigma_crit * tangentialShear(shapes_lens[config['shape_ra_key']], shapes_lens[config['shape_dec_key']], shapes_lens[config['shape_e1_key']], -shapes_lens[config['shape_e2_key']], lens['RA'], lens['DEC'], computeB=True)
            if config['coords'] == "physical":
                data['radius'] = Ang2Dist(np.array(d12), lens[lens_z_key])
            else:
                data['radius'] = d12
            data['weight'] = getValues(shapes_lens, config['shape_weight_key'])/sigma_crit**2

            # get indices for all sources in each slice
            for key, limit in config['splittings'].iteritems():
                slice_str = np.repeat('-', n_gal)
                values = getValues(shapes_lens, key)
                for s in xrange(len(limit)-1):
                    mask = getSliceMask(values, limit[s], limit[s+1])
                    slice_str[mask] = '%d' % s
                    del mask
                if key == config['splittings'].keys()[0]:
                    data['slices'] = slice_str
                else:
                    data['slices'] = np.core.defchararray.add(data['slices'], slice_str)
                del values, slice_str
            
            # write out results
            if counter == 0:
                fits.write(data)
            else:
                fits[-1].append(data)
            counter += 1
            done += n_gal

        fits.close()
        print "done. Created " + stackfile

    else:
        print "stackfile " + stackfile + " already exists."
        print "Delete it or use different label to rerun this script."
        exit(0)


