#!/bin/env python

import numpy as np
import esutil as eu
import os, errno
from sys import argv
from shear_stacking import *
import math 
import fitsio
import json
from glob import glob

def getValues(s, key, functions):
    # what values are used for the slices: functions are direct columns?
    if key in functions.keys():
        return eval(functions[key])
    else:
        return s[key]

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
    # parse inputs
    try:
        configfile = argv[1]
        startindex = int(argv[2])
    except IndexError:
        print "usage: " + argv[0] + " <config file> <start index>"
        raise SystemExit
    try:
        fp = open(configfile)
        print "opening configfile " + configfile
        config = json.load(fp)
        fp.close()
    except IOError:
        print "configfile " + configfile + " does not exist!"
        raise SystemExit

    lensfile = config['lens_catalog']
    outdir = os.path.dirname(configfile)
    if outdir[-1] != '/':
        outdir += '/'
    try:
        os.makedirs(outdir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(outdir):
            pass
        else: raise

    if config['coords'] not in ['angular', 'physical']:
        print "config: specify either 'angular' or 'physical' coordinates"
        raise SystemExit
    
    # open shapes, apply post-run selections
    print "opening " + shapefile
    shapefile = config['shape_file']
    chunk_size = config['shape_chunk_size']
    shdu = fitsio.FITS(shapefile)
    if len(config['shape_cuts']) == 0:
        shapes = shdu[1][startindex : startindex+chunk_size]
    else:
        # that's not really elegant since the .where runs on entire table
        cuts = " && ".join(config['shape_cuts'])
        mask = shdu[1].where(cuts)[startindex : startindex+chunk_size]
        if mask.size:
            shapes = shdu[1][mask]
        else:
            shapes = np.array([])
        del mask
    print "shape sample: %d" % shapes.size

    
    if shapes.size:
        basename = os.path.basename(shapefile)
        basename = basename.split(".")[0]
        stackfile = outdir + basename + '_DeltaSigma_%d.fits' % startindex
        matchfile = '/tmp/' + basename + '_matches_%d.bin' % startindex

        # open lens catalog, apply selection if desired
        hdu = fitsio.FITS(lensfile)
        if len(config['lens_cuts']) == 0:
            lenses = hdu[1][:]
        else:
            cuts = " && ".join(config['lens_cuts'])
            mask = hdu[1].where(cuts)
            if mask.size:
                lenses = hdu[1][mask]
            else:
                lenses = np.array([])
            del mask
        print "lens sample: %d" % lenses.size

        # find all galaxies in shape catalog within maxrange arcmin 
        # of each lens center
        if lenses.size and shapes.size:
            print "matching lens and source catalog..."
            maxrange = float(config['maxrange'])
            if config['coords'] == "physical":
                maxrange = Dist2Ang(maxrange, lenses[config['lens_z_key']])

            h = eu.htm.HTM(8)
            matchfile = matchfile.encode('ascii') # htm.match expects ascii filenames
            h.match(lenses['RA'], lenses['DEC'], shapes[config['shape_ra_key']], shapes[config['shape_dec_key']], maxrange, maxmatch=-1, file=matchfile)
            htmf = HTMFile(matchfile)
            Nmatch = htmf.n_matches
            print "  found ", Nmatch, "matches"
        else:
            Nmatch = 0

        if Nmatch:
            print "stacking lenses..."
            fits = fitsio.FITS(tmpstackfile, 'rw')
            data = np.empty(Nmatch, dtype=[('radius_angular', 'f4'), ('radius_physical', 'f4'), ('DeltaSigma', 'f8'), ('DeltaSigma_x', 'f8'), ('sensitivity', 'f8'), ('weight', 'f8'), ('slices', '%di1' % len(config['splittings']))])
            specz_calib = getSpecZCalibration()
            done = 0

            # iterate over all lenses, write DeltaSigma, r, weight into file
            for m1, m2, d12 in htmf.matches():
                lens = lenses[m1]
                shapes_lens = shapes[m2]
                n_gal = shapes_lens.size

                # compute effective Sigma_crit
                z_phot, cz = getSigmaCritCorrection(specz_calib, lens[config['lens_z_key']])
                sigma_crit = getSigmaCritEffective(z_phot, cz, shapes_lens[config['shape_z_key']])
                # compute tangential and cross shear
                gt, gx = tangentialShear(shapes_lens[config['shape_ra_key']], shapes_lens[config['shape_dec_key']], getValues(shapes_lens, config['shape_e1_key'], config['functions']), getValues(shapes_lens, config['shape_e2_key'], config['functions']), lens['RA'], lens['DEC'], computeB=True)

                data['DeltaSigma'][done:done+n_gal] = sigma_crit * gt
                data['DeltaSigma_x'][done:done+n_gal] = sigma_crit * gx
                data['radius_angular'][done:done+n_gal] = d12
                data['radius_physical'][done:done+n_gal] = Ang2Dist(np.array(d12), lens[config['lens_z_key']])

                # compute sensitivity and weights
                data['sensitivity'][done:done+n_gal] = getValues(shapes_lens, config['shape_sensitivity_key'], config['functions'])
                data['weight'][done:done+n_gal] = getValues(shapes_lens, config['shape_weight_key'], config['functions'])/sigma_crit**2

                # get indices for all sources in each slice
                i = 0
                for key, limit in config['splittings'].iteritems():
                    data['slices'][done:done+n_gal][:,i] = -1 # null value
                    values = getValues(shapes_lens, key, config['functions'])
                    for s in xrange(len(limit)-1):
                        mask = getSliceMask(values, limit[s], limit[s+1])
                        data['slices'][done:done+n_gal][:,i][mask] = s
                        del mask
                    i += 1
                    del values
                done += n_gal
                del lens, shapes_lens, z_phot, cz, sigma_crit, gt, gx
            os.system('rm ' + matchfile)

            fits.write(data)
            fits.close()
            print "done. Created " + stackfile

        hdu.close()
    shdu.close()

