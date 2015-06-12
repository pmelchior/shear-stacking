#!/bin/env python

import os, errno, json, fitsio, copy
import numpy as np
from sys import argv
from shear_stacking import *
from glob import glob
from multiprocessing import Pool, cpu_count

def computeMeanStdForProfile(profile, name, s=None, n_jack=0):
    # use build-in method to calculate in-bin means and dispersions
    if n_jack == 0:
        if s is None:
            mean_r, n, mean_q, std_q, sum_w = profile[name].getProfile()
        else:
            mean_r, n, mean_q, std_q, sum_w = profile[name][s].getProfile()
        mask = (n>0)
        return mean_r[mask], n[mask], mean_q[mask], std_q[mask], sum_w[mask]
    else: # jackknife
        q = []
        missing = []
        for i in xrange(n_jack):
            if s is None:
                r_, n_, q_, std_q, sum_w = profile[i][name].getProfile()
            else:
                r_, n_, q_, std_q, sum_w = profile[i][name][s].getProfile()
            missing.append(n_ == 0)
            q.append(q_)
        missing = np.array(missing)
        q = np.ma.masked_array(q, mask=missing)
        mean_q = q.mean(axis=0)

        # result for normal/non-jackknife profile
        if s is None:
            mean_r, n, mean0, std_q, sum_w = profile[-1][name].getProfile()
        else:
            mean_r, n, mean0, std_q, sum_w = profile[-1][name][s].getProfile()
        mask = (n>0)

        # variance and bias-corrected mean needs number of actual jackknifes:
        # to be corrected for available data in each radial bin
        n_avail = n_jack - missing.sum(axis=0)
        mean_q = n_avail*mean0 - (n_avail - 1)*mean_q
        std_q = ((n_avail - 1.)/n_avail * ((q - mean_q)**2).sum(axis=0))**0.5
        return mean_r[mask], n[mask], mean_q.data[mask], std_q.data[mask], sum_w[mask]

def insertIntoProfile(data, profile, config):
    # total profiles (E and B mode)
    profile['E'].insert(data['radius_' + config['coords']], data['DeltaSigma'], data['weight'], data['sensitivity'])
    profile['B'].insert(data['radius_' + config['coords']], data['DeltaSigma_x'], data['weight'], data['sensitivity'])
    # get index list of matching objects from each splitting
    i = 0
    for key, limit in config['splittings'].iteritems():
        for s in xrange(len(limit)-1):
            mask = (data['slices'][:,i] == s)
            profile[key][s].insert(data['radius_' + config['coords']][mask], data['DeltaSigma'][mask], data['weight'][mask], data['sensitivity'][mask])
            del mask
        i += 1

def jackknife_masks(ids, region, n_jack):
    # jackknife region for each cluster
    rids = np.zeros(ids.size, dtype='uint8')
    start = 0
    stop = 0
    current_id = ids[0]
    for i in xrange(1, ids.size):
        stop = i
        if ids[i] != current_id:
            rids[start:stop] = region[current_id]
            current_id = ids[i]
            start = i
    # last segment cannot complete
    rids[start:] = region[current_id]

    # allow everything but region i
    masks = []
    for i in xrange(n_jack):
        masks.append(rids != i)
    return masks

def readIntoProfile(stackfile, profile, config, n_jack=0, region=None):
    print "opening file " + stackfile
    try:
        fits = fitsio.FITS(stackfile)
        data = fits[1].read()

        if config['coords'] == "angular":
            data['radius_angular'] *= 60

        if n_jack < 2:
            insertIntoProfile(data, profile, config)
        else:
            # select only data from lenses in given region
            masks = jackknife_masks(data['lens_index'], region, n_jack=n_jack)
            for i in xrange(n_jack):
                insertIntoProfile(data[masks[i]], profile[i], config)
            # for bias correction: add normal profile, no jackknifing
            insertIntoProfile(data, profile[-1], config)
            del masks

        # clean up
        fits.close()
        del data
    except:
        print "  trouble opening! skipping..."
    return profile

def createProfiles(bins, config, n_jack=0):
    # create empty profiles for each of the jackknife regions
    # add one extra profile for the normal sample (all pairs without jackknife)
    if n_jack > 1:
        profile = []
        for i in xrange(n_jack+1):
            profile.append(createProfiles(bins, config))
        return profile

    # create profile for E/B mode of all data
    # and for each slice defined
    profile = {'E': BinnedScalarProfile(bins),
               'B': BinnedScalarProfile(bins)
           }
    for key, limit in config['splittings'].iteritems():
        profile[key] = []
        for s in xrange(len(limit)-1):
            profile[key].append(BinnedScalarProfile(bins))
    return profile

# simple helper function: for each key in profile structure, append profile2
def appendProfile(profile, profile2, config, n_jack=0):
    # for jackknifes: split into regions
    if n_jack > 1:
        for i in xrange(n_jack+1):
            appendProfile(profile[i], profile2[i], config)
    else:
        profile['E'] += profile2['E']
        profile['B'] += profile2['B']
        for key, limit  in config['splittings'].iteritems():
            for s in xrange(len(limit)-1):
                profile[key][s] += profile2[key][s] 

if __name__ == '__main__':
    # parse inputs
    try:
        configfile = argv[1]
        coords = argv[2]
    except IndexError:
        print "usage: " + argv[0] + " <config file> <angular/physical> [jackknife regions]"
        raise SystemExit
    try:
        fp = open(configfile)
        print "opening configfile " + configfile
        config = json.load(fp)
        fp.close()
    except IOError:
        print "configfile " + configfile + " does not exist!"
        raise SystemExit
    
    if coords not in ['angular', 'physical']:
        print "specify either angular or physical coordinates"
        raise SystemExit

    if len(argv) > 3:
        n_jack = int(argv[3])
    else:
        n_jack = 0

    indir = os.path.dirname(configfile) + "/"
    outdir = indir
    stackfiles = glob(indir + '/*_DeltaSigma*.fits')
    if len(stackfiles) == 0:
        print "run stack_slices.py before!"
        raise SystemExit

    if coords == "physical":
        bins =  np.exp(0.3883*np.arange(-10, 12))
    else:
        bins = np.arange(1,11,1)

    # if jacknife errors are desired: create jackknife regions from
    # the lens file by k-means clustering and assign each lens
    # to the nearest k-means center
    region = None
    if n_jack:
        print "defining %d jackknife regions" % n_jack
        import kmeans_radec
        lenses = getLensCatalog(config, verbose=True)
        radec = np.dstack((lenses[config['lens_ra_key']], lenses[config['lens_dec_key']]))[0]
        maxiter=100
        tol=1.0e-5
        km=kmeans_radec.kmeans_sample(radec, n_jack, maxiter=maxiter, tol=tol)
        if not km.converged:
            raise RuntimeError("k means did not converge")
        
        # define regions: ids of lenses assigned to each k-means cluster
        region  = km.find_nearest(radec)

    # set up containers
    profile = createProfiles(bins, config, n_jack=n_jack)
    # need a deep copy of the empty structure for each of the threads
    initprofile = copy.deepcopy(profile)

    # iterate thru all DeltaSigma files: distribute them over given radial bins
    n_processes = min(6, cpu_count())
    if n_jack:
        n_processes = cpu_count()
    pool = Pool(processes=n_processes)
    results = [pool.apply_async(readIntoProfile, (stackfile, initprofile, config, n_jack, region)) for stackfile in stackfiles]
    for r in results:
        profile_ = r.get()
        appendProfile(profile, profile_, config, n_jack=n_jack)

    # save profiles to npz
    filename = outdir + "shear_profile_%s_EB.npz" % coords
    mean_r, n, mean_e, std_e, sum_w = computeMeanStdForProfile(profile, 'E', n_jack=n_jack)
    mean_r, n, mean_b, std_b, sum_w = computeMeanStdForProfile(profile, 'B', n_jack=n_jack)
    kwargs = {"mean_r": mean_r, "n": n, "mean_e": mean_e, "std_e": std_e, "mean_b": mean_b, "std_b": std_b, "sum_w": sum_w, "n_jack": n_jack }
    np.savez(filename, **kwargs)
    print "writing " + filename
    
    for key, limit  in config['splittings'].iteritems():
        for s in xrange(len(limit)-1):
            mean_r, n, mean_e, std_e, sum_w = computeMeanStdForProfile(profile, key, s=s, n_jack=n_jack)
            filename = outdir + "shear_profile_%s_%s_%d.npz" % (coords, key, s)
            kwargs = {"mean_r": mean_r, "n": n, "mean_e": mean_e, "std_e": std_e, "sum_w": sum_w, "n_jack": n_jack }
            np.savez(filename, **kwargs)
            print "writing " + filename
