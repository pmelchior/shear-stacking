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
    
    # open shapes file(s)
    shapefile = config['shape_file']
    chunk_size = config['shape_chunk_size']
    shdu = fitsio.FITS(shapefile)
    extra = None
    print "opening shapefile %s (%d entries)" % (shapefile, shdu[1].get_nrows())

    if len(config['shape_cuts']) == 0:
        shapes = shdu[1][chunk_index*chunk_size : (chunk_index+1)*chunk_size]
        try:
            ehdu = fitsio.FITS(config['shape_file_extra'])
            print "opening extra shapefile " + config['shape_file_extra']
            extra = ehdu[1][chunkindex*chunk_size : (chunkindex+1)*chunk_size]
            ehdu.close()
        except KeyError:
            pass
    else:
    # apply shape cuts: either on the file itself of on the extra file
    # since we're working with FITS type selections, we can't apply it
    # directly to the shapes array, but need to go back to the catalogs.
    # that's not really elegant since the .where runs on entire table
        cuts = " && ".join(config['shape_cuts'])
        try:
            ehdu = fitsio.FITS(config['shape_file_extra'])
            print "opening extra shapefile " + config['shape_file_extra']
            mask = ehdu[1].where(cuts)
            print "selecting %d shapes" % mask.size
            mask = mask[chunk_index*chunk_size : (chunk_index+1)*chunk_size]
            shapes = shdu[1][mask]
            extra = ehdu[1][mask]
            ehdu.close()
        except KeyError:
            mask = shdu[1].where(cuts)
            print "selecting %d shapes" % mask.size
            shapes = shdu[1][mask[chunk_index*chunk_size : (chunk_index+1)*chunk_size]]
        del mask
    print "shape sample (this chunk): %d" % shapes.size
    
    # if there's an extra file: join data with shapes
    if extra is not None:
        from numpy.lib import recfunctions
        columns = shapes.dtype.names
        ecolumns = []
        for col in extra.dtype.names:
            if col not in columns:
                ecolumns.append(col)
        shapes = recfunctions.rec_append_fields(shapes, ecolumns, [extra[c] for c in ecolumns])

    if shapes.size:
        basename = os.path.basename(shapefile)
        basename = basename.split(".")[0]
        stackfile = outdir + basename + '_DeltaSigma_%d.fits' % chunk_index
        matchfile = '/tmp/' + basename + '_matches_%d.bin' % chunk_index

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
            fits = fitsio.FITS(stackfile, 'rw', clobber=True)
            data = np.empty(Nmatch, dtype=[('radius_angular', 'f4'), ('radius_physical', 'f4'), ('DeltaSigma', 'f8'), ('DeltaSigma_x', 'f8'), ('sensitivity', 'f8'), ('weight', 'f8'), ('slices', '%di1' % len(config['splittings']))])
            
            # get effective lensing weights (w * Sigma_crit ^-1 or -2)
            wz1 = getWZ(power=1)
            wz2 = getWZ(power=2)

            # iterate over all lenses, write DeltaSigma, r, weight into file
            done = 0
            for m1, m2, d12 in htmf.matches():
                lens = lenses[m1]
                shapes_lens = shapes[m2]
                n_gal = shapes_lens.size

                # compute tangential and cross shear
                gt, gx = tangentialShear(shapes_lens[config['shape_ra_key']], shapes_lens[config['shape_dec_key']], getValues(shapes_lens, config['shape_e1_key'], config['functions']), getValues(shapes_lens, config['shape_e2_key'], config['functions']), lens['RA'], lens['DEC'], computeB=True)

                data['DeltaSigma'][done:done+n_gal] = gt
                data['DeltaSigma_x'][done:done+n_gal] = gx
                data['radius_angular'][done:done+n_gal] = d12
                data['radius_physical'][done:done+n_gal] = Ang2Dist(np.array(d12), lens[config['lens_z_key']])

                # compute sensitivity and weights: with the photo-z bins, we use
                # DeltaSigma = wz1 < gt> / (wz2 <s>),
                # where wz1 and wz2 are the effective weights at given lens z
                zs_bin = getValues(shapes_lens, config['shape_z_key'], config['functions'])
                sensitivity = getValues(shapes_lens, config['shape_sensitivity_key'], config['functions'])
                for b in xrange(3): # 3 bins
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
                del lens, shapes_lens, gt, gx, zs_bin, sensitivity
            os.system('rm ' + matchfile)

            fits.write(data)
            fits.close()
            print "done. Created " + stackfile

        hdu.close()
    shdu.close()

