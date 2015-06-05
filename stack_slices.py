#!/bin/env python

import os, errno, json
import numpy as np
import esutil as eu
from sys import argv
from glob import glob
from shear_stacking import *

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

def getQuadrantMask(lens, shapes, config):
    quad_flags = lens['quad_flags']
    # no quadrant pairs OK
    if quad_flags <= 1:
        return np.array([], dtype='bool')
    # all OK
    if quad_flags == (2**0 + 2**1 + 2**2 + 2**3 + 2**4):
        return np.ones(shapes.size, dtype='bool')

    # not all quadrant pairs are OK
    # FIXME: needs to be defined from angles on the curved sky
    # NOTE: very restrictive due to strict ordering of quadrants
    # e.g. if all sources are in the top half, but quad_flags & 2 > 0,
    # then no source will be selected even if quad_flags & 8 > 0
    if quad_flags & 2 > 0:
        return shapes[config['shape_dec_key']] > lens[config['lens_dec_key']]
    if quad_flags & 4 > 0:
        return shapes[config['shape_ra_key']] < lens[config['lens_ra_key']]
    if quad_flags & 8 > 0:
        return shapes[config['shape_dec_key']] < lens[config['lens_dec_key']]
    if quad_flags & 16 > 0:
        return shapes[config['shape_ra_key']] > lens[config['lens_ra_key']]
        
if __name__ == '__main__':
    # parse inputs
    try:
        configfile = argv[1]
        chunk_index = int(argv[2])
    except IndexError:
        print "usage: " + argv[0] + " <config file> <chunk index>"
        raise SystemExit
    try:
        fp = open(configfile)
        print "opening configfile " + configfile
        config = json.load(fp)
        fp.close()
    except IOError:
        print "configfile " + configfile + " does not exist!"
        raise SystemExit

    if config['coords'] not in ['angular', 'physical']:
        print "config: specify either 'angular' or 'physical' coordinates"
        raise SystemExit

    # open shape catalog
    outdir = os.path.dirname(configfile) + "/"
    shapefile = config['shape_file']
    shapes = getShapeCatalog(config, verbose=True, chunk_index=chunk_index)

    if shapes.size:
        basename = os.path.basename(shapefile)
        basename = ".".join(basename.split(".")[:-1])
        stackfile = outdir + basename + '_DeltaSigma_%d.fits' % chunk_index
        matchfile = '/tmp/' + basename + '_matches_%d.bin' % chunk_index

        # open lens catalog
        lenses = getLensCatalog(config, verbose=True)
        # do we have the column for the quadrant check?
        do_quadrant_check = 'quad_flags' in lenses.dtype.names

        # find all galaxies in shape catalog within maxrange arcmin 
        # of each lens center
        if lenses.size:
            print "matching lens and source catalog..."
            maxrange = float(config['maxrange'])
            if config['coords'] == "physical":
                maxrange = Dist2Ang(maxrange, lenses[config['lens_z_key']])

            h = eu.htm.HTM(8)
            matchfile = matchfile.encode('ascii') # htm.match expects ascii filenames
            h.match(lenses[config['lens_ra_key']], lenses[config['lens_dec_key']], shapes[config['shape_ra_key']], shapes[config['shape_dec_key']], maxrange, maxmatch=-1, file=matchfile)
            htmf = HTMFile(matchfile)
            Nmatch = htmf.n_matches
            print "  found ", Nmatch, "matches"
        else:
            Nmatch = 0

        if Nmatch:
            print "stacking lenses..."
            fits = fitsio.FITS(stackfile, 'rw', clobber=True)
            data = np.empty(Nmatch, dtype=[('lens_index', 'i8'), ('radius_angular', 'f4'), ('radius_physical', 'f4'), ('DeltaSigma', 'f8'), ('DeltaSigma_x', 'f8'), ('sensitivity', 'f8'), ('weight', 'f8'), ('slices', '%di1' % len(config['splittings']))])
            
            # get effective lensing weights (w * Sigma_crit ^-1 or -2)
            wz1 = getWZ(power=1)
            wz2 = getWZ(power=2)

            # iterate over all lenses, write DeltaSigma, r, weight into file
            done = 0
            for m1, m2, d12 in htmf.matches():
                lens = lenses[m1]
                shapes_lens = shapes[m2]
                # check which sources around a lens we can use
                if do_quadrant_check:
                    mask = getQuadrantMask(lens, shapes_lens, config)
                    shapes_lens = shapes_lens[mask]
                    d12 = np.array(d12)[mask]
                n_gal = shapes_lens.size

                if n_gal:
                    # compute tangential and cross shear
                    gt, gx = tangentialShear(shapes_lens[config['shape_ra_key']], shapes_lens[config['shape_dec_key']], getValues(shapes_lens, config['shape_e1_key'], config['functions']), getValues(shapes_lens, config['shape_e2_key'], config['functions']), lens[config['lens_ra_key']], lens[config['lens_dec_key']], computeB=True)

                    data['lens_index'][done:done+n_gal] = m1
                    data['DeltaSigma'][done:done+n_gal] = gt
                    data['DeltaSigma_x'][done:done+n_gal] = gx
                    data['radius_angular'][done:done+n_gal] = d12
                    data['radius_physical'][done:done+n_gal] = Ang2Dist(np.array(d12), lens[config['lens_z_key']])

                    # compute sensitivity and weights: with the photo-z bins, we use
                    # DeltaSigma = wz1 < gt> / (wz2 <s>),
                    # where wz1 and wz2 are the effective weights at given lens z
                    zs_bin = getValues(shapes_lens, config['shape_z_key'], config['functions'])
                    sensitivity = getValues(shapes_lens, config['shape_sensitivity_key'], config['functions'])
                    for b in np.unique(zs_bin):
                        wz1_ = extrap(lens[config['lens_z_key']], wz1['z'], wz1['bin%d' % b])
                        wz2_ = extrap(lens[config['lens_z_key']], wz2['z'], wz2['bin%d' % b])
                        mask = zs_bin == b
                        data['weight'][done:done+n_gal][mask] = wz1_
                        # in the plot_slices script, the sensitivity will be
                        # multiplied with the weight, so we need to divide wz1 out
                        data['sensitivity'][done:done+n_gal][mask] = sensitivity[mask] * wz2_ / wz1_
                        del mask

                    # get indices for all sources in each slice
                    i = 0
                    for key, limit in config['splittings'].iteritems():
                        data['slices'][done:done+n_gal][:,i] = -1 # null value
                        if config['split_type'] == 'shape':
                            values = getValues(shapes_lens, key, config['functions'])
                            for s in xrange(len(limit)-1):
                                mask = getSliceMask(values, limit[s], limit[s+1])
                                data['slices'][done:done+n_gal][:,i][mask] = s
                                del mask
                            del values
                        elif config['split_type'] == 'lens':
                            value = getValues(lens, key, config['functions'])
                            for s in xrange(len(limit)-1):
                                # each lens can only be in one slice per key
                                if getSliceMask(value, limit[s], limit[s+1]):
                                    data['slices'][done:done+n_gal][:,i] = s
                                    break
                        i += 1

                    done += n_gal
                    del shapes_lens, gt, gx, zs_bin, sensitivity

            # finish up
            os.system('rm ' + matchfile)
            fits.write(data[:done])
            fits.close()
            print "done. Created " + stackfile
