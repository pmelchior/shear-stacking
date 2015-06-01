#!/bin/env python

import os, errno, json, fitsio, copy
import matplotlib
matplotlib.use('agg')
import pylab as plt
import numpy as np
import esutil as eu
from sys import argv
from shear_stacking import *
from glob import glob
from multiprocessing import Pool

# colors based on blue/white/red divergent colormap
# from Kevin Moreland:
# http://www.sandia.gov/~kmorel/documents/ColorMaps/
# To emphasize the mid-range, I used a darker midpoint of 0.33 instead of 0.88
# split is the length of the splitting list
def getColors(split):
    colors = [(0.23137254901960785, 0.29803921568627451, 0.75294117647058822, 1.0), (0.70588235294117652, 0.015686274509803921, 0.14901960784313725, 1.0)]
    if split < 3:
        raise AssertionError("Splitting must at least have two separate bins") 
    if split == 4:
        colors.insert(1, (0.7803921568627451, 0.7803921568627451, 0.7803921568627451, 1.0))
    if split == 5:
        colors.insert(1, (0.71372549019607845, 0.70196078431372544, 0.90588235294117647, 1.0))
        colors.insert(2, (0.92941176470588238, 0.65490196078431373, 0.63137254901960782, 1.0))
    if split == 6:
        colors.insert(1, (0.62745098039215685, 0.61568627450980395, 0.91764705882352937, 1.0))
        colors.insert(2, (0.7803921568627451, 0.7803921568627451, 0.7803921568627451, 1.0))
        colors.insert(3, (0.93333333333333335, 0.53725490196078429, 0.50980392156862742, 1.0))
    if split > 6:
        raise NotImplementedError("Splittings > 5 are not implemented")
    return colors

def makeEBProfile(ax, key_name, profile_E, profile_B, coords, lw=1):
    mean_r, n, mean_q, std_q = profile_B.getProfile()
    ax.errorbar(mean_r, mean_q, yerr=std_q, c='r', marker='.', label='B-mode', lw=lw)
    mean_r, n, mean_q, std_q = profile_E.getProfile()
    ax.errorbar(mean_r, mean_q, yerr=std_q, c='k', marker='.', label='E-mode', lw=lw)

    xlim = ax.get_xlim()
    if coords == "angular":
        ax.plot(xlim, [0,0], 'k:')
        ax.set_xlim(xlim)
    ax.legend(loc='upper right', numpoints=1, frameon=False, fontsize='small')
    ax.set_ylabel(r'$\Delta\Sigma\ [10^{14}\ \mathrm{M}_\odot \mathrm{Mpc}^{-2}]$')
    if coords == "physical":
        ax.set_xlabel('Radius [Mpc/$h$]')
        ax.set_xscale('symlog', linthreshx=1e-2)
        ax.xaxis.set_minor_locator(matplotlib.ticker.LogLocator(subs=np.arange(2, 10)))
        ax.set_yscale('symlog', linthreshy=1e-3)
        ax.yaxis.set_minor_locator(matplotlib.ticker.LogLocator(subs=np.arange(2, 10)))
    else:
        ax.set_xlabel('Radius [arcmin]')
    return mean_r, n, mean_q, std_q
    

def makeSlicedProfile(ax, key_name, profiles, limits, coords, xlim, ylim, lw=1):
    if coords == "angular":
        ax.plot(xlim, [0,0], 'k:')
    if plt.matplotlib.rcParams['text.usetex']:
        title = r'\texttt{' + key_name.replace("_", "\_") + '}'
    else:
        title = key_name

    # make the profile for each split
    colors = getColors(len(limits))
    for s in range(len(limits)-1):
        label = '$\in ['
        if isinstance(limits[s], (int, long)):
            label += '%d, ' % limits[s]
        else:
            label += '%.2f, ' % limits[s]
        if isinstance(limits[s+1], (int, long)):
            label += '%d)$' % limits[s+1]
        else:
            label += '%.2f)$' % limits[s+1]
        mean_r_, n_, mean_q_, std_q_ = profiles[s].getProfile()
        ax.errorbar(mean_r_, mean_q_, yerr=std_q_, c=colors[s], marker='.', label=label, lw=lw)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    legend = ax.legend(loc='upper right', numpoints=1, title=title, frameon=False, fontsize='small')
    plt.setp(legend.get_title(),fontsize='small')
    ax.set_ylabel(r'$\Delta\Sigma\ [10^{14}\ \mathrm{M}_\odot \mathrm{Mpc}^{-2}]$')
    if coords == "physical":
        ax.set_xlabel('Radius [Mpc/$h$]')
        ax.set_xscale('symlog', linthreshx=1e-2)
        ax.xaxis.set_minor_locator(matplotlib.ticker.LogLocator(subs=np.arange(2, 10)))
        ax.set_yscale('symlog', linthreshy=1e-3)
        ax.yaxis.set_minor_locator(matplotlib.ticker.LogLocator(subs=np.arange(2, 10)))
    else:
        ax.set_xlabel('Radius [arcmin]')

def readNbin(stackfile, profile, coords):
    print "opening file " + stackfile
    fits = fitsio.FITS(stackfile)
    data = fits[1].read()

    if coords == "angular":
        data['radius_angular'] *= 60

    # total profiles (E and B mode)
    profile['all_E'].insert(data['radius_' + config['coords']], data['DeltaSigma'], data['weight'], data['sensitivity'])
    profile['all_B'].insert(data['radius_' + config['coords']], data['DeltaSigma_x'], data['weight'], data['sensitivity'])
    # get index list of matching objects from each splitting
    i = 0
    for key, limit in config['splittings'].iteritems():
        for s in xrange(len(limit)-1):
            mask = (data['slices'][:,i] == s)
            profile[key][s].insert(data['radius_' + config['coords']][mask], data['DeltaSigma'][mask], data['weight'][mask], data['sensitivity'][mask])
            del mask
        i += 1

    # clean up
    fits.close()
    del data
    return profile

if __name__ == '__main__':
    # parse inputs
    try:
        configfile = argv[1]
        coords = argv[2]
    except IndexError:
        print "usage: " + argv[0] + " <config file> <angular/physical> [outdir]"
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

    indir = os.path.dirname(configfile)
    if len(argv) > 3:
        outdir = argv[3]
    else:
        outdir = indir
    if outdir[-1] != '/':
        outdir += '/'
    try:
        os.makedirs(outdir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(outdir):
            pass
        else: raise

    stackfiles = glob(indir + '/*_DeltaSigma*.fits')
    if len(stackfiles) == 0:
        print "run stack_slices.py before!"
        raise SystemExit

    if coords == "physical":
        bins =  np.exp(0.3883*np.arange(-12, 12))
    else:
        bins = np.arange(1,11,1) #config['maxrange']*60, 2)

    # set up containers
    profile = {'all_E': BinnedScalarProfile(bins),
               'all_B': BinnedScalarProfile(bins)
               }
    for key, limit in config['splittings'].iteritems():
        profile[key] = []
        for s in xrange(len(limit)-1):
            profile[key].append(BinnedScalarProfile(bins))
    initprofile = copy.deepcopy(profile)
    keys = config['splittings'].keys()

    # iterate thru all DeltaSigma files
    n_processes = min(len(stackfiles), 8)
    pool = Pool(processes=n_processes)
    results = [pool.apply_async(readNbin, (stackfile, initprofile, coords)) for stackfile in stackfiles]
    for r in results:
        thisprofile = r.get()
        profile['all_E'] += thisprofile['all_E']
        profile['all_B'] += thisprofile['all_B']
        for key, limit  in config['splittings'].iteritems():
            for s in xrange(len(limit)-1):
                profile[key][s] += thisprofile[key][s] 
        
    # plot generation: E/B profile
    plotfile = outdir + 'shear_stack_EB_' + coords + '.png'
    print plotfile
    setTeXPlot(sampling=2)
    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(111)
    mean_r, n, mean_q, std_q = makeEBProfile(ax, 'all', profile['all_E'], profile['all_B'], coords)
    pivot = (mean_q + std_q/2)[n > 0].max()
    if coords == "physical":
        xlim = (1e-2, mean_r[n > 0].max()*2)
        ylim = (1e-3, 2*pivot)
    else:
        xlim = ax.get_xlim()
        ylim = (-0.15*pivot, 1.25*pivot)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    fig.subplots_adjust(wspace=0, hspace=0, left=0.16, bottom=0.13, right=0.98, top=0.97)
    fig.savefig(plotfile)

    # sliced profile plots
    for key, limit in config['splittings'].iteritems():
        print "  " + key
        plotfile = outdir + 'shear_stack_' + key + '_' + coords + '.png'
        fig = plt.figure(figsize=(5, 4))
        ax = fig.add_subplot(111)
        makeSlicedProfile(ax, key, profile[key], config['splittings'][key], coords, xlim, ylim)
        fig.subplots_adjust(wspace=0, hspace=0, left=0.16, bottom=0.13, right=0.98, top=0.97)
        fig.savefig(plotfile)

