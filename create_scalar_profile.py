#!/bin/env python

import fitsio, numpy as np, json
import esutil as eu
from esutil import numpy_util
from shear_stacking import *
from sys import argv
from multiprocessing import Pool, current_process, cpu_count
from glob import glob
import pylab as plt


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

def createProfile(config):
    n_jack = config['n_jack']
    if config['coords'] == "physical":
        bins =  np.exp(0.3883*np.arange(-10, 10))
    else:
        bins = np.arange(1,11,1)
    
    # create profile for all data and for each slice defined
    pnames = ['all']
    for key, limit in config['splittings'].iteritems():
        for s in xrange(len(limit)-1):
            pnames.append("%s_%d" % (key, s))

    profile = {}
    # each slice get a binned profile for all data and each jackknife region
    for pname in pnames:
        profile[pname] = [BinnedScalarProfile(bins) for i in xrange(n_jack + 1)]
    return profile

# BinnedScalarProfile has a += operator, so we just have to call this for
# each slices/jackknifed profile
def appendToProfile(profile, profile_):
    for pname in profile.keys():
        for i in xrange(len(profile[pname])):
            profile[pname][i] += profile_[pname][i]

def insertIntoProfile(profile, pname, radius, q, w, region=-1, mask=None):
    if mask is None:
        for i in xrange(len(profile[pname])):
            if i != region:
                profile[pname][i].insert(radius, q, w)
    else:
        for i in xrange(len(profile[pname])):
            if i != region:
                profile[pname][i].insert(radius[mask], q[mask], w[mask])

def stackShapes(shapes, lenses, config, regions):
    chunk_size = config['shape_chunk_size']
    shapefile = config['shape_file']
    thread_id = current_process()._identity
    basename = os.path.basename(shapefile)
    basename = ".".join(basename.split(".")[:-1])
    matchfile = '/tmp/' + basename + '_matches_%d.bin' % thread_id

    # do we have the column for the quadrant check?
    do_quadrant_check = 'quad_flags' in lenses.dtype.names

    # find all galaxies in shape catalog within maxrange arcmin 
    # of each lens center
    maxrange = float(config['maxrange'])
    if config['coords'] == "physical":
        maxrange = Dist2Ang(maxrange, lenses[config['lens_z_key']])

    h = eu.htm.HTM(8)
    matchfile = matchfile.encode('ascii') # htm.match expects ascii filenames
    h.match(lenses[config['lens_ra_key']], lenses[config['lens_dec_key']], shapes[config['shape_ra_key']], shapes[config['shape_dec_key']], maxrange, maxmatch=-1, file=matchfile)
    htmf = HTMFile(matchfile)
    Nmatch = htmf.n_matches

    # profile container
    profile = createProfile(config)

    if Nmatch:
        # iterate over all lenses, write scalar value, r, weight into file
        for m1, m2, d12 in htmf.matches():
            lens = lenses[m1]
            region = regions[m1]
            shapes_lens = shapes[m2]
            
            # check which sources around a lens we can use
            if do_quadrant_check:
                mask = getQuadrantMask(lens, shapes_lens, config)
                shapes_lens = shapes_lens[mask]
                d12 = np.array(d12)[mask]
                del mask
            n_gal = shapes_lens.size

            if n_gal:
                q = getValues(shapes_lens, config['shape_scalar_key'], config['functions'])
                w = getValues(shapes_lens, config['shape_weight_key'], config['functions'])
                radius = np.array(d12)
                if config['coords'] == "physical":
                    radius = Ang2Dist(d12, lens[config['lens_z_key']])

                # save unsliced profile first
                insertIntoProfile(profile, 'all', radius, q, w, region=region)
                    
                # find out in which slice each pair falls
                for key, limit in config['splittings'].iteritems():
                    if config['split_type'] == 'shape':
                        values = getValues(shapes_lens, key, config['functions'])
                        for s in xrange(len(limit)-1):
                            pname = "%s_%d" % (key, s)
                            mask = getSliceMask(values, limit[s], limit[s+1])
                            insertIntoProfile(profile, pname, radius, q, w, region=region, mask=mask)
                            del mask
                        del values
                        
                    elif config['split_type'] == 'lens':
                        value = getValues(lens, key, config['functions'])
                        for s in xrange(len(limit)-1):
                            pname = "%s_%d" % (key, s)
                            # each lens can only be in one slice per key
                            if getSliceMask(value, limit[s], limit[s+1]):
                                mask = None
                                insertIntoProfile(profile, pname, radius, q, w, region=region, mask=mask)
                                break
                            
                del shapes_lens, q, w, radius

    # finish up
    os.system('rm ' + matchfile)
    return profile

def getJackknifeRegions(config, lenses, outdir):
    # if jacknife errors are desired: create jackknife regions from
    # the lens file by k-means clustering and assign each lens
    # to the nearest k-means center
    # If reuse_jack is specified: reload previously generated centers
    # to use fixed regions
    if config['n_jack']:
        import kmeans_radec
        jack_file = outdir + "km_centers.npy"
        radec = np.dstack((lenses[config['lens_ra_key']], lenses[config['lens_dec_key']]))[0]
        if not os.path.exists(jack_file):
            print "defining %d jackknife regions" % n_jack
            maxiter = 100
            tol = 1.0e-5
            km = kmeans_radec.kmeans_sample(radec, n_jack, maxiter=maxiter, tol=tol)
            if not km.converged:
                raise RuntimeError("k means did not converge")

            # save result for later
            try:
                os.makedirs(outdir_jack)
            except OSError:
                pass
            np.save(jack_file, km.centers)

        else:
            print "reusing jackknife regions from " + jack_file
            centers_ = np.load(jack_file)
            km = kmeans_radec.KMeans(centers_)
        
        # define regions: ids of lenses assigned to each k-means cluster
        regions  = km.find_nearest(radec)
    else:
        # do not use regions: -1 will never be selected for jackknifes
        regions = -1 * np.ones(len(lenses), dtype='int8')
    return regions

def computeMeanStdForProfile(profile):
    n_jack = len(profile)-1
    # use build-in method to calculate in-bin means and dispersions
    if n_jack == 0:
        mean_r, n, mean_q, std_q, sum_w = profile.getProfile()
        mask = (n>0)
        return {"mean_r": mean_r[mask], "n": n[mask], "mean_q": mean_q[mask], "std_q": std_q[mask], "sum_w": sum_w[mask], "n_jack": n_jack }
    
    else: # jackknife
        q = []
        missing = []
        for i in xrange(n_jack):
            r_, n_, q_, std_q, sum_w = profile[i].getProfile()
            missing.append(n_ == 0)
            q.append(q_)
        missing = np.array(missing)
        q = np.ma.masked_array(q, mask=missing)
        mean_q = q.mean(axis=0)

        # result for normal/non-jackknife profile
        mean_r, n, mean0, std_q, sum_w = profile[-1].getProfile()
        mask = (n>0)

        # variance and bias-corrected mean needs number of actual jackknifes:
        # to be corrected for available data in each radial bin
        n_avail = n_jack - missing.sum(axis=0)
        mean_q = n_avail*mean0 - (n_avail - 1)*mean_q
        std_q = ((n_avail - 1.)/n_avail * ((q - mean_q)**2).sum(axis=0))**0.5

    return {"mean_r": mean_r[mask], "n": n[mask], "mean_q": mean_q.data[mask], "std_q": std_q.data[mask], "sum_w": sum_w[mask], "n_jack": n_jack }
                
def collapseJackknifes(profile):
    for pname in profile.keys():
        profile[pname] = computeMeanStdForProfile(profile[pname])


if __name__ == '__main__':
    # parse inputs
    try:
        configfile = argv[1]
    except IndexError:
        print "usage: " + argv[0] + " <scalar config file>"
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


    outdir = os.path.dirname(configfile) + "/"
    scalarfiles = outdir + 'scalar_' + config['shape_scalar_key'] + '*.npz'

    # only do something if there are no profiles present
    if len(glob(scalarfiles)) == 0:
    
        # open shape catalog
        shapefile = config['shape_file']
        shapes_all = getShapeCatalog(config, verbose=True, chunk_index=None)

        # open lens catalog
        lenses = getLensCatalog(config, verbose=True)

        if shapes_all.size and lenses.size:

            # container to hold profiles
            profile = createProfile(config)

            # get the jackknife regions (if specified in config)
            regions = getJackknifeRegions(config, lenses, outdir)
            
            # cut into manageable junks and distribute over cpus
            print "running lens-source stacking ..."
            n_processes = cpu_count()
            pool = Pool(processes=n_processes)
            chunk_size = config['shape_chunk_size']
            splits = shapes_all.size/chunk_size
            results = [pool.apply_async(stackShapes, (shapes, lenses, config, regions)) for shapes in np.array_split(shapes_all, splits)]
            i = 0
            for r in results:
                profile_ = r.get()
                appendToProfile(profile, profile_)
                print "  job %d/%d done" % (i, splits)
                i+=1

            # collapse jackknifes into means and stds
            print "aggregating results..."
            collapseJackknifes(profile)
            
            # save profiles
            for pname in profile.keys():
                filename = outdir + 'scalar_' + config['shape_scalar_key'] + '_' + pname + '.npz'
                print "writing " + filename
                np.savez(filename, **(profile[pname]))

    else:
        print "Scalar profiles " + scalarfiles + " already exist."
        
        

                
                






