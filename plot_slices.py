#!/bin/env python

import os, json, fitsio
import matplotlib
matplotlib.use('agg')
import pylab as plt
import numpy as np
import esutil as eu
from sys import argv
from common import *
from glob import glob
from matplotlib.ticker import NullFormatter

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
    if split > 5:
        raise NotImplementedError("Splittings > 5 are not implemented")
    return colors

def makeEBProfile(ax, bins, key_name, profile_E, profile_B, coords, lw=1):
    xlim = [bins[0], bins[-1] + 2]
    ax.plot(xlim, [0,0], 'k:')
    if plt.matplotlib.rcParams['text.usetex']:
        label = r'\texttt{' + key_name.replace("_", "\_") + '}'
    else:
        label = key_name
    mean_r, n, mean_q, std_q = profile_B.getProfile()
    ax.errorbar(mean_r, mean_q, yerr=std_q, c='r', marker='.', label='E-mode', lw=lw)
    mean_r, n, mean_q, std_q = profile_E.getProfile()
    ax.errorbar(mean_r, mean_q, yerr=std_q, c='k', marker='.', label='B-mode', lw=lw)  
    ax.legend(loc='upper right', numpoints=1, frameon=False, fontsize='x-small')
    ax.set_ylabel(r'$\Delta\Sigma\ [10^{14}\ \mathrm{M}_\odot \mathrm{Mpc}^{-2}]$')
    if coords == "physical":
        ax.set_xlabel('Radius [Mpc/$h$]')
    else:
        ax.set_xlabel('Radius [arcmin]')
    return mean_r, n, mean_q, std_q
    

def makeSlicedProfile(ax, bins, key_name, profiles, limits, coords, ylim, lw=1):
    xlim = [bins[0], bins[-1] + 2]
    ax.plot(xlim, [0,0], 'k:')
    if plt.matplotlib.rcParams['text.usetex']:
        label = r'\texttt{' + key_name.replace("_", "\_") + '}'
    else:
        label = key_name

    # make the profile for each split
    colors = getColors(len(limits))
    for s in range(len(limits)-1):
        string = '$\in ['
        if isinstance(limits[s], (int, long)):
            string += '%d, %d)$'
        else:
            string += '%.2f, %.2f)$'
        label = string % (limits[s], limits[s+1])
        mean_r_, n_, mean_q_, std_q_ = profiles[s].getProfile()
        ax.errorbar(mean_r_, mean_q_, yerr=std_q_, c=colors[s], marker='.', label=label, lw=lw)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.legend(loc='upper right', numpoints=1, frameon=False, fontsize='x-small')
    ax.set_ylabel(r'$\Delta\Sigma\ [10^{14}\ \mathrm{M}_\odot \mathrm{Mpc}^{-2}]$')
    if coords == "physical":
        ax.set_xlabel('Radius [Mpc/$h$]')
    else:
        ax.set_xlabel('Radius [arcmin]')

if __name__ == '__main__':
    if len(argv) < 2:
        print "usage: " + argv[0] +  "<config file> [outdir]"
        exit(0)

    configfile = argv[1]
    if os.path.exists(configfile) is False:
        print "configfile " + configfile + " does not exist!"
        exit(0)

    print "opening configfile " + configfile
    fp = open(configfile)
    config = json.load(fp)

    if len(argv) > 2:
        outdir = argv[2]
    else:
        outdir = os.path.dirname(configfile)
    if outdir[-1] != '/':
        outdir += '/'

    stackfiles = glob(outdir + '*_DeltaSigma.fits')[::2]
    if len(stackfiles) == 0:
        print "run stack_slices.py before!"
        exit(0)

    if config['coords'] == "physical":
        bins = np.arange(0, config['maxrange'], 0.5)
    else:
        bins = np.arange(1, config['maxrange']*60, 5)

    # set up containers
    profile = {'all_E': BinnedScalarProfile(bins),
               'all_B': BinnedScalarProfile(bins)
               }
    for key, limit in config['splittings'].iteritems():
        profile[key] = []
        for s in xrange(len(limit)-1):
            profile[key].append(BinnedScalarProfile(bins))
    keys = config['splittings'].keys()

    # iterate thru all DeltaSigma files
    for stackfile in stackfiles:
        print "opening file " + stackfile
        fits = fitsio.FITS(stackfile)
        data = fits[1].read()

        if config['coords'] == "angular":
            data['radius'] *= 60
        
        # total profiles (E and B mode)
        profile['all_E'].insert(data['radius'], data['DeltaSigma'], data['weight'])
        profile['all_B'].insert(data['radius'], data['DeltaSigma_x'], data['weight'])
        # get index list of matching objects from each splitting
        i = 0
        for key, limit in config['splittings'].iteritems():
            for s in xrange(len(limit)-1):
                mask = (data['slices'][:,i] == s)
                profile[key][s].insert(data['radius'][mask], data['DeltaSigma'][mask], data['weight'][mask])
                del mask
            i += 1

        # clean up
        fits.close()
        del data
        
    # Plot generation
    plotfile = outdir + 'shear_stack_EB.pdf'
    setTeXPlot(sampling=2)
    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(111)
    mean_r, n, mean_q, std_q = makeEBProfile(ax, bins, 'all', profile['all_E'], profile['all_B'], config['coords'])
    pivot = (mean_q + std_q/2).max()
    ylim = (-0.15*pivot, 1.25*pivot)
    ax.set_ylim(ylim)
    fig.subplots_adjust(wspace=0, hspace=0, left=0.13, bottom=0.13, right=0.97, top=0.97)
    fig.savefig(plotfile)
    for key, limit in config['splittings'].iteritems():
        print "  " + key
        plotfile = outdir + 'shear_stack_' + key + '.pdf'
        fig = plt.figure(figsize=(5, 4))
        ax = fig.add_subplot(111)
        makeSlicedProfile(ax, bins, key, profile[key], config['splittings'][key], config['coords'], ylim)
        fig.subplots_adjust(wspace=0, hspace=0, left=0.13, bottom=0.13, right=0.97, top=0.97)
        fig.savefig(plotfile)

