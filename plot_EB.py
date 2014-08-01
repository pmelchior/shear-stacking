#!/bin/env python

import matplotlib
matplotlib.use('agg')
import pylab as plt
import numpy as np
import pyfits
from sys import argv
from os.path import exists
from common import *
from matplotlib.ticker import NullFormatter


def makeEBProfilePlot(ax, bins, radius, DeltaSigma, DeltaSigma_cross, weight, lw=1):
    xlim = [bins[0], bins[-1] + 2]
    ax.plot(xlim, [0,0], 'k:')
    mean_r, n, mean_q, std_q = scalarProfile(bins, radius, DeltaSigma, weight)
    ax.errorbar(mean_r, mean_q, yerr=std_q, c='k', marker='.', label='E-mode', lw=lw)
    pivot = (mean_q + std_q/2).max()
    mean_r, n, mean_q, std_q = scalarProfile(bins, radius, DeltaSigma_cross, weight)
    ax.errorbar(mean_r, mean_q, yerr=std_q, c='r', marker='.', label='B-mode', lw=lw)

    ax.set_ylim(-0.1*pivot, 1.1*pivot)
    ax.set_xlim(xlim)
    ax.legend(loc='upper right', numpoints=1, frameon=False, fontsize='small')

if __name__ == '__main__':
    if len(argv) < 3:
        print "usage: " + argv[0] + " <band> <output label>"
        exit(1)

    band = argv[1]
    label = argv[2]
    matchfile = '/tmp/matches_' + band.lower() + '.bin'
    stackfile = '/tmp/shear_stack_' + band.lower() + '_' + label + '.npz'
    plotfile = 'shear_stack_EB_' + band.lower() + '_' + label + '.pdf'

    if exists(stackfile) is False:
        print "run stack_slices.py before!"
        exit(0)
    else:
        # load from previously saved file
        print "loading from file", stackfile
        X = np.load(stackfile)
        DeltaSigma = X['DeltaSigma']
        DeltaSigma_cross = X['DeltaSigma_cross']
        weight = X['weight']
        radius = X['radius']
        maxrange = np.ceil(radius.max())
        
    # Plot generation
    from TeXConfig import *
    setTeXPlot(sampling=2)
    print "generating plot..."
    fig = plt.figure(figsize=(5, 4))
    bins = np.arange(1, maxrange, 2)
    ax = fig.add_subplot(111)
    makeEBProfilePlot(ax, bins, radius, DeltaSigma, DeltaSigma_cross, weight)
    ax.set_xlabel('Radius [arcmin]')
    ax.set_ylabel(r'$\Delta\Sigma\ [10^{14}\ \mathrm{M}_\odot \mathrm{Mpc}^{-2}]$')
    plt.subplots_adjust(wspace=0, hspace=0, left=0.16, bottom=0.13, right=0.99, top=0.97)
    plt.savefig(plotfile)
    print "done. open", plotfile

