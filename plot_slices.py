#!/bin/env python

import matplotlib
matplotlib.use('agg')
import pylab as plt
import numpy as np
import esutil as eu
import pyfits
from sys import argv
from os.path import exists
from common import *
from glob import glob
from matplotlib.ticker import NullFormatter

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

def makeSlicedProfilePlot(ax, bins, radius, DeltaSigma, weight, slices, splittings, key_name, lw=1, mean_profile=None):
    xlim = [bins[0], bins[-1] + 2]
    ax.plot(xlim, [0,0], 'k:')
    if plt.matplotlib.rcParams['text.usetex']:
        label = r'\texttt{' + key_name.replace("_", "\_") + '}'
    else:
        label = key_name
    if mean_profile is None:
        mean_r, n, mean_q, std_q = scalarProfile(bins, radius, DeltaSigma, weight)
    else:
        mean_r, n, mean_q, std_q = mean_profile
    ax.errorbar(mean_r, mean_q, yerr=std_q, c='k', marker='.', label=label, lw=lw)

    # make the profile for each split
    colors = getColors(len(splittings))
    for s in range(len(splittings)-1):
        string = '$\in ['
        if isinstance(splittings[s], (int, long)):
            string += '%d, %d)$'
        else:
            string += '%.2f, %.2f)$'
        label = string % (splittings[s], splittings[s+1])
        if slices[s].size:
            mean_r_, n_, mean_q_, std_q_ = scalarProfile(bins, radius[slices[s]], DeltaSigma[slices[s]], weight[slices[s]])
            ax.errorbar(mean_r_, mean_q_, yerr=std_q_, c=colors[s], marker='.', label=label, lw=lw)
    pivot = (mean_q + std_q/2).max()
    ax.set_ylim(-0.15*pivot, 1.25*pivot)
    ax.set_xlim(xlim)
    ax.legend(loc='upper right', numpoints=1, frameon=False, fontsize='small')
    return mean_r, n, mean_q, std_q

if __name__ == '__main__':
    if len(argv) < 3:
        print "usage: " + argv[0] + " <band> <output label>"
        exit(1)

    band = argv[1]
    label = argv[2]
    matchfile = '/tmp/matches_' + band.lower() + '.bin'
    stackfile = '/tmp/shear_stack_' + band.lower() + '_' + label + '.npz'
    plotfile = 'shear_stack_' + band.lower() + '_' + label + '.pdf'

    if exists(stackfile) is False:
        print "run stack_slices.py before!"
        exit(0)
    else:
        # load from previously saved file
        print "loading from file", stackfile
        X = np.load(stackfile)
        DeltaSigma = X['DeltaSigma']
        weight = X['weight']
        radius = X['radius']
        maxrange = np.ceil(radius.max())
        splittings = X['splittings']
        slices = {}
        for i in range(len(keys)):
            if callable(keys[i]):
                key_name = keys[i].__name__
            else:
                key_name = keys[i]
            slices[key_name] = {}
            if iterable(splittings[i]):
                for s in range(len(splittings[i])-1):
                    slices[key_name][s] = X[key_name + "_%d" % s]
            else:
                for s in range(splittings[i]):
                    slices[key_name][s] = X[key_name + "_%d" % s]
        
    # Plot generation
    setTeXPlot(sampling=2)
    print "generating plots..."
    maxcol = 4 # how many columns per row
    rows = (len(keys)-1)/maxcol + 1
    fig = plt.figure(figsize=(12, 4*rows))
    bins = np.arange(1, maxrange, 2)
    mean_profile = None
    for i in range(len(keys)):
        if callable(keys[i]):
            key_name = keys[i].__name__
        else:
            key_name = keys[i]
        ax = fig.add_subplot(rows, min(maxcol, len(keys)), i+1)
        if mean_profile is None:
            mean_profile = makeSlicedProfilePlot(ax, bins, radius, DeltaSigma, weight, slices[key_name], splittings[i], key_name)
        else:
            makeSlicedProfilePlot(ax, bins, radius, DeltaSigma, weight, slices[key_name], splittings[i], key_name, mean_profile=mean_profile)
        if i/maxcol == rows-1:
            ax.set_xlabel('Radius [arcmin]')
        if i%maxcol == 0:
            ax.set_ylabel(r'$\Delta\Sigma\ [10^{14} \mathrm{M}_\odot \mathrm{Mpc}^{-2}]$')
        else:
            ax.yaxis.set_major_formatter(NullFormatter())

    # size of the bottom margin
    if rows == 1:
        bottom = 0.13
    else:
        bottom = 0.07
    plt.subplots_adjust(wspace=0, hspace=0, left=0.07, bottom=bottom, right=0.99, top=0.97)
    plt.savefig(plotfile)
    print "done. open", plotfile

