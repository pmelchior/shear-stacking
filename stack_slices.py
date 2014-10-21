#!/bin/env python

import numpy as np
import esutil as eu
import pyfits
from sys import argv
from os.path import exists
from os import system
from common import *
from glob import glob
from math import pi
import fitsio

def getValues(data, key):
    # what values are used for the slices
    if callable(key):
        return key(data)
    else:
        return data[key]

def getSliceMask(values, lower, upper, return_num=False):
    if return_num is False:
        return (values >= lower) & (values < upper)
    else:
        return sum((values >= lower) & (values < upper))

def saveStack(outdir, DeltaSigma, DeltaSigma_cross, radius, weight, keys, splittings, slices, last_element):
    np.save(outdir + "/DeltaSigma.npy", DeltaSigma[:last_element['all']])
    np.save(outdir + "/DeltaSigma_cross.npy", DeltaSigma_cross[:last_element['all']])
    np.save(outdir + "/weight.npy", weight[:last_element['all']])
    if coords == "physical":
        np.save(outdir + "/radius_physical.npy", radius[:last_element['all']])
    else:
        np.save(outdir + "/radius_angular.npy", radius[:last_element['all']])
    keynames = []
    for k in keys:
        if callable(k):
            key_name = k.__name__
        else:
            key_name = k
        keynames.append(key_name)
        for vv in slices[key_name].keys():
            filename = key_name + "_%d" % vv + ".npy"
            np.save(outdir + "/" + filename, slices[key_name][vv][:last_element[key_name][vv]])
    np.save(outdir + "/splittings.npy", splittings)
    np.savetxt(outdir + "/keynames.txt", keynames, fmt='%s')

# get separation in deg for distance L in Mpc/h at redshift z
# uses c/H0 = 3000 Mpc/h
def Dist2Ang(L, z):
    global cosmo
    return L / cosmo.Da(z) / 3000. * 180./math.pi

def Ang2Dist(theta, z):
    global cosmo
    return theta * cosmo.Da(z) * 3000. / 180. * math.pi

def B_D(data):
    return data['im3shape_' + band.lower() + '_bulge_flux'] / data['im3shape_' + band.lower() + '_disc_flux']

if __name__ == '__main__':
    if len(argv) < 5:
        print "usage: " + argv[0] + " <lens catalog> <shape catalog> <band> <output label> [tmpdir]"
        exit(1)

    lensfile = argv[1]
    shapefile = argv[2]
    band = argv[3]
    label = argv[4]
    if len(argv) > 5:
        tmpdir = argv[5]
    else:
        tmpdir = "/tmp/"
    
    coords = "angular"
    lens_z_key = 'Z_LAMBDA'
    shape_z_key = 'ZP'
    shape_ra_key = 'ALPHAWIN_J2000_' + band.upper()
    shape_dec_key = 'DELTAWIN_J2000_' + band.upper()
    
    keys = [shape_z_key, 'im3shape_' + band.lower() + '_snr', 'im3shape_' + band.lower() + '_radius', 'im3shape_' + band.lower() + '_stamp_size', 'im3shape_' + band.lower() + '_bulge_flux', 'im3shape_' + band.lower() + '_disc_flux', B_D, 'im3shape_' + band.lower() + '_mask_fraction']
    splittings = [[0.7, 0.9, 1.1, 1.5], [20,40,60,1000], [0.263,0.789,26.3], [32,48,64,128], 3, 3, 3, [0., 0.2, 0.4, 1]]
    
    #keys = ['im3shape_' + band.lower() + '_info_flag']
    #splittings = [[0,1,8,128,1024, 2**21]]

    outdir = tmpdir + "/" + label
    system('mkdir -p ' + outdir)

    if coords == "physical":
        maxrange = 5. # Mpc/h
    else:
        maxrange = 1.1  # deg

    matchfile = outdir + '/matches_' + band.lower() + '.bin'
    stackfile = outdir + "/DeltaSigma.fits"

    if exists(stackfile) is False:
        # open lens catalog, apply selection if desired
        hdu = pyfits.open(lensfile)
        lenses = hdu[1].data
        #good_cl = (lenses[lens_z_key] < 0.6)
        #lenses = lenses[good_cl]
        
        if coords == "physical":
            maxrange = Dist2Ang(maxrange, lenses[lens_z_key])
        print "lens sample: %d" % lenses.size

        # open shapes, apply post-run selections
        shdu = pyfits.open(shapefile)
        good_sh = ModestSG(shdu[1].data) & (shdu[1].data['im3shape_' + band.lower() + '_exists'] == 1) & (shdu[1].data['im3shape_' + band.lower() + '_error_flag'] == 0) & (shdu[1].data['FLAGS_' + band.upper()] == 0)
        shapes = shdu[1].data[good_sh]
        print "shape sample: %d" % shapes.size

        # find all galaxies in shape catalog within maxrange arcmin 
        # of each lens center
        print "matching lens and source catalog..."
        if exists(matchfile) is False:
            # CAVEAT: make sure to have enough space where you put the match file
            # it has 24 byte per match, which quickly becomes Gb's of data 
            h = eu.htm.HTM(8)
            h.match(lenses['RA'], lenses['DEC'], shapes[shape_ra_key], shapes[shape_dec_key], maxrange, maxmatch=-1, file=matchfile)
            del h
        else:
            print "  re-using existing matchfile", matchfile

        htmf = HTMFile(matchfile)
        Ngal = htmf.n_matches
        print "  found ", Ngal, "matches"

        # define the slice splittings
        # can be any key in the shape catalog or a function thereof
        print "determining slice ranges..."
        keynames = []
        for i in range(len(keys)):
            if callable(keys[i]):
                key_name = keys[i].__name__
            else:
                key_name = keys[i]
            keynames.append(key_name)
            if hasattr(splittings[i], '__iter__') is False:
                # remove 2.5% on either side to reduce impact of outliers
                delta = 95./splittings[i]
                ranges = [2.5 + k*delta for k in range(splittings[i]+1)]
                # FIXME: there is no selection on the shapes here, e.g. no BG cut
                values = getValues(shapes, keys[i])
                splittings[i] = percentile(values, ranges)
                del ranges, values
            print "  " + key_name + ":", splittings[i]
        np.save(outdir + "/splittings.npy", splittings)
        np.savetxt(outdir + "/keynames.txt", keynames, fmt='%s')

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
            z_phot, cz = getSigmaCritCorrection(specz_calib, lens[lens_z_key])
            sigma_crit = getSigmaCritEffective(z_phot, cz, shapes_lens[shape_z_key])
            # determine extent in DeltaSigma array
            data = np.empty(shapes_lens.size, dtype=[('radius', 'f8'), ('DeltaSigma', 'f8'), ('DeltaSigma_x', 'f8'), ('weight', 'f8'), ('slices', 'S%d' % len(keys))])
            data['DeltaSigma'], data['DeltaSigma_x'] = sigma_crit * tangentialShear(shapes_lens[shape_ra_key], shapes_lens[shape_dec_key], shapes_lens['im3shape_' + band.lower() + '_e1'], -shapes_lens['im3shape_' + band.lower() + '_e2'], lens['RA'], lens['DEC'], computeB=True)
            if coords == "physical":
                data['radius'] = Ang2Dist(np.array(d12), lens[lens_z_key])
            else:
                data['radius'] = d12
            data['weight'] = 0.2/(0.2**2 + (0.1*20/shapes_lens['im3shape_' + band.lower() + '_snr'])**2)**0.5/sigma_crit**2

            # get indices for all sources in each slice
            for i in xrange(len(keys)):
                slice_str = np.repeat('-', n_gal)
                values = getValues(shapes_lens, keys[i])
                for s in xrange(len(splittings[i])-1):
                    mask = getSliceMask(values, splittings[i][s], splittings[i][s+1])
                    slice_str[mask] = '%d' % s
                    del mask
                if i == 0:
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

            # output status update and do backup
            if counter % 100 == 0:
                print '  lens %d, matched %.2f%%' % (counter, done*100./htmf.n_matches)
            del data, shapes_lens, z_phot, cz, sigma_crit

        fits.close()

    else:
        print "stackfile " + stackfile + " already exists."
        print "Delete it or use different label to rerun this script."
        exit(0)


