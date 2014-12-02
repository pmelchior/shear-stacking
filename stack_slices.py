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
    stackfile = outdir + basename + '_DeltaSigma.fits'
    tmpstackfile = '/tmp/' + basename + '_DeltaSigma.fits'
    matchfile = '/tmp/' + basename + '_matches.bin'
    configfile = outdir + 'config.json'

    if os.path.exists(configfile) is False:
        print "configfile " + configfile + " does not exist!"
        exit(0)

    print "opening configfile " + configfile
    fp = open(configfile)
    config = json.load(fp)
    fp.close()

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
        
        maxrange = config['maxrange']
        if config['coords'] == "physical":
            maxrange = Dist2Ang(maxrange, lenses[config['lens_z_key']])

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
        h = eu.htm.HTM(8)
        h.match(lenses['RA'], lenses['DEC'], shapes[config['shape_ra_key']], shapes[config['shape_dec_key']], maxrange, maxmatch=-1, file=matchfile)
        del h
        htmf = HTMFile(matchfile)
        Nmatch = htmf.n_matches
        print "  found ", Nmatch, "matches"

        # iterate over all lenses, write DeltaSigma, r, weight into file
        if Nmatch:
            print "stacking lenses..."
            fits = fitsio.FITS(tmpstackfile, 'rw')
            data = np.empty(Nmatch, dtype=[('radius_angular', 'f4'), ('radius_physical', 'f4'), ('DeltaSigma', 'f8'), ('DeltaSigma_x', 'f8'), ('weight', 'f8'), ('slices', '%di1' % len(config['splittings']))])
            specz_calib = getSpecZCalibration()
            done = 0
            for m1, m2, d12 in htmf.matches():
                lens = lenses[m1]
                shapes_lens = shapes[m2]
                n_gal = shapes_lens.size

                # compute effective Sigma_crit
                z_phot, cz = getSigmaCritCorrection(specz_calib, lens[config['lens_z_key']])
                sigma_crit = getSigmaCritEffective(z_phot, cz, shapes_lens[config['shape_z_key']])
                # determine extent in DeltaSigma array
                
                data['DeltaSigma'][done:done+n_gal], data['DeltaSigma_x'][done:done+n_gal] = sigma_crit * tangentialShear(shapes_lens[config['shape_ra_key']], shapes_lens[config['shape_dec_key']], shapes_lens[config['shape_e1_key']], -shapes_lens[config['shape_e2_key']], lens['RA'], lens['DEC'], computeB=True)
                data['radius_angular'][done:done+n_gal] = d12
                data['radius_physical'][done:done+n_gal] = Ang2Dist(np.array(d12), lens[config['lens_z_key']])
                data['weight'][done:done+n_gal] = getValues(shapes_lens, config['shape_weight_key'])/sigma_crit**2

                # get indices for all sources in each slice
                i = 0
                for key, limit in config['splittings'].iteritems():
                    data['slices'][done:done+n_gal][:,i] = -1 # null value
                    values = getValues(shapes_lens, key)
                    for s in xrange(len(limit)-1):
                        mask = getSliceMask(values, limit[s], limit[s+1])
                        data['slices'][done:done+n_gal][:,i][mask] = s
                        del mask
                    i += 1
                    del values
                done += n_gal

            fits.write(data)
            fits.close()
            os.system('mv ' + tmpstackfile + ' ' + stackfile)
            print "done. Created " + stackfile
        os.system('rm ' + matchfile)
        hdu.close()
        shdu.close()
    else:
        print "stackfile " + stackfile + " already exists."
