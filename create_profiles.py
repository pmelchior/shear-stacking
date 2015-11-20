#!/bin/env python

import fitsio, numpy as np, json
import esutil as eu
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
    l = config['minrange']
    u = config['maxrange']
    n_bins = config['n_bins']
    bin_type = config['bin_type']

    if bin_type == "linear":
        bins = np.linspace(l, u, n_bins+1)
    elif bin_type == "log":
        dlogx = (np.log(u) - np.log(l))/n_bins
        bins = np.exp(np.log(l) + dlogx * np.arange(n_bins+1))
    else:
        raise NotImplementedError("bin_type %s not in ['linear', 'log']" % bin_type)
    
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

def insertIntoProfile(profile, pname, R, Q, W, S, region=-1, mask=None):
    if mask is None:
        for i in xrange(len(profile[pname])):
            if i != region:
                profile[pname][i].insert(R, Q, W, S=S)
    else:
        for i in xrange(len(profile[pname])):
            if i != region:
                if S is not None:
                    profile[pname][i].insert(R[mask], Q[mask], W[mask], S=S[mask])
                else:
                    profile[pname][i].insert(R[mask], Q[mask], W[mask], S=S)

def getShearValues(shapes_lens, lens, config):
    global wz1, wz2
    
    # compute tangential shear
    gt = tangentialShear(shapes_lens[config['shape_ra_key']], shapes_lens[config['shape_dec_key']], getValues(shapes_lens, config['shape_e1_key'], config['functions']), getValues(shapes_lens, config['shape_e2_key'], config['functions']), lens[config['lens_ra_key']], lens[config['lens_dec_key']], computeB=False)

    # compute DeltaSigma from source redshift
    # we use DeltaSigma = wz2 * wz1**-1 < gt> / (wz2 <s>),
    # where wz1 and wz2 are the effective inverse Sigma_crit
    # weights at given lens z
    """W = getValues(shapes_lens, config['shape_weight_key'], config['functions'])
    z_l = getValues(lens, config['lens_z_key'], config['functions'])
    z_s = getValues(shapes_lens, config['shape_z_key'], config['functions'])
    Sigma_crit = getSigmaCrit(z_l, z_s)
    mask = z_s > z_l
    gt[mask] *= Sigma_crit[mask]**-1
    W[mask] *= Sigma_crit[mask]**-2
    W[mask == False] = 0
    """

    # compute sensitivity and weights: with the photo-z bins,
    # the precomputed wz1 already contains the measurement weights,
    # we just need to apply the effective Sigma_crit to gt and
    # replace W with wz1**2 (dropping the measurement weight which would otherwise
    # be counted twice).
    # See Sheldon et al., 2004, AJ, 127, 2544 (eq. 19)
    W = np.zeros(gt.size) # better safe than sorry
    zs_bin = getValues(shapes_lens, config['shape_z_key'], config['functions'])
    for b in np.unique(zs_bin):
        wz1_ = extrap(lens[config['lens_z_key']], wz1['z'], wz1['bin%d' % b])
        mask = zs_bin == b
        gt[mask] /= wz1_
        W[mask] = wz1_**2
    
    S = getValues(shapes_lens, config['shape_sensitivity_key'], config['functions'])
    return gt, W, S

def stackShapes(shapes, lenses, profile_type, config, regions):
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

                # define the profile quantities: radius, q, weight, sensitivity
                if config['coords'] == "physical":
                    R = Ang2Dist(d12, lens[config['lens_z_key']])
                else:
                    R = np.array(d12)
                
                if profile_type == "scalar":
                    Q = getValues(shapes_lens, config['shape_scalar_key'], config['functions'])
                    W = getValues(shapes_lens, config['shape_weight_key'], config['functions'])
                    S = None
                if profile_type == "shear":
                    Q, W, S = getShearValues(shapes_lens, lens, config)

                # save unsliced profile first
                insertIntoProfile(profile, 'all', R, Q, W, S, region=region)
                    
                # find out in which slice each pair falls
                for key, limit in config['splittings'].iteritems():
                    if config['split_type'] == 'shape':
                        values = getValues(shapes_lens, key, config['functions'])
                        for s in xrange(len(limit)-1):
                            pname = "%s_%d" % (key, s)
                            mask = getSliceMask(values, limit[s], limit[s+1])
                            insertIntoProfile(profile, pname, R, Q, W, S, region=region, mask=mask)
                            del mask
                        del values
                        
                    elif config['split_type'] == 'lens':
                        value = getValues(lens, key, config['functions'])
                        for s in xrange(len(limit)-1):
                            pname = "%s_%d" % (key, s)
                            # each lens can only be in one slice per key
                            if getSliceMask(value, limit[s], limit[s+1]):
                                mask = None
                                insertIntoProfile(profile, pname, R, Q, W, S, region=region, mask=mask)
                                break
                            
                del shapes_lens, R, Q, W, S

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
        n_jack = config['n_jack']
        import kmeans_radec
        jack_file = outdir + "n_jack/km_centers.npy"
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
                os.makedirs(outdir + "n_jack")
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
        return {"mean_r": mean_r[mask], "n": n[mask], "mean_q": mean_q[mask], "std_q": std_q[mask], "sum_w": sum_w[mask]}
    
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

    return {"mean_r": mean_r[mask], "n": n[mask], "mean_q": mean_q.data[mask], "std_q": std_q.data[mask], "sum_w": sum_w[mask]}
                
def collapseJackknifes(profile):
    for pname in profile.keys():
        profile[pname] = computeMeanStdForProfile(profile[pname])


if __name__ == '__main__':
    # parse inputs
    try:
        configfile = argv[1]
        profile_type = argv[2]
    except IndexError:
        print "usage: " + argv[0] + " <config file> <shear/scalar>"
        raise SystemExit
    try:
        fp = open(configfile)
        print "opening configfile " + configfile
        config = json.load(fp)
        fp.close()
    except IOError:
        print "configfile " + configfile + " does not exist!"
        raise SystemExit

    if profile_type not in ['shear', 'scalar']:
        print "specify profile_type from ['shear', 'scalar']"
        raise SystemExit    
    
    if config['coords'] not in ['angular', 'physical']:
        print "config: specify either 'angular' or 'physical' coordinates"
        raise SystemExit


    outdir = os.path.dirname(configfile) + "/"
    if profile_type == "shear":
        name = "shear_"
    if profile_type == "scalar":
        name = "scalar_" + config['shape_scalar_key'] + "_"  
    profile_files = outdir + name + '*.npz'

    # only do something if there are no profiles present
    if len(glob(profile_files)) == 0:
    
        # open shape catalog
        shapes_all = getShapeCatalog(config, verbose=True)
        if shapes_all.size == 0:
            print "Shape catalog empty"
            raise SystemExit
        
        # open lens catalog
        lenses = getLensCatalog(config, verbose=True)
        if lenses.size == 0:
            print "Lens catalog empty"
            raise SystemExit

        # container to hold profiles
        profile = createProfile(config)

        # get the jackknife regions (if specified in config)
        regions = getJackknifeRegions(config, lenses, outdir)

        # load lensing weights (w * Sigma_crit ^-1 or -2) for shear profiles
        if profile_type == "shear":
            wz1 = getWZ(power=1)

        # cut into manageable junks and distribute over cpus
        print "running lens-source stacking ..."
        n_processes = cpu_count()
        pool = Pool(processes=n_processes)
        chunk_size = config['shape_chunk_size']
        splits = len(shapes_all)/chunk_size
        if len(shapes_all) % chunk_size != 0:
            splits += 1
        results = [pool.apply_async(stackShapes, (shapes_all[i*chunk_size: (i+1)*chunk_size], lenses, profile_type, config, regions)) for i in range(splits)]
        j = 0
        for r in results:
            profile_ = r.get()
            appendToProfile(profile, profile_)
            r, n, mean_q, std_q, sum_w = profile['all'][-1].getProfile()
            print "  job %d/%d (n_pairs = %.3fe9) done" % (j, splits, n.sum() / 1e9)
            j+=1

        # save jackknife region results
        if config['n_jack']:
            print "saving jackknife profiles..."
            for pname in profile.keys():
                for i in xrange(len(profile[pname])):
                    filename = outdir + 'n_jack/' + name + pname + '_%d.npz' % i
                    profile[pname][i].save(filename)

        # collapse jackknifes into means and stds
        print "aggregating results..."
        collapseJackknifes(profile)

        # save profiles
        for pname in profile.keys():
            filename = outdir + name + pname + '.npz'
            print "writing " + filename
            np.savez(filename, **(profile[pname]))

        # print all profile to stdout
        p = profile['all']
        print "\nALL profile:"
        print "{0:>8s} | {1:>12s} | {2:>12s} | {3:>12s} +- {4:>12s}".format("RADIUS", "NUMBER", "SUM(W)/AREA", "MEAN", "STD")
        print "-" * 70
        for i in xrange(len(p['n'])):
            print "{0:8.2f} | {1:12g} | {2:12g} | {3:12g} +- {4:12g}".format(p['mean_r'][i], p['n'][i], p['sum_w'][i], p['mean_q'][i], p['std_q'][i])

    else:
        print "Profiles " + profile_files + " already exist."
        
        

                
                






