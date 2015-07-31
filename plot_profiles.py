#!/bin/env python

import os, errno, json, fitsio, copy
import matplotlib
matplotlib.use('agg')
import pylab as plt
import numpy as np
from sys import argv
from shear_stacking import *
from glob import glob
from multiprocessing import Pool, cpu_count

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

def makeEBProfile(ax, profile, plot_type, lw=1):
    colors = getColors(3)
    if plot_type == "shear":
        ax.errorbar(profile['mean_r'], profile['mean_e'], yerr=profile['std_e'], c='k', marker='.', label='E-mode', lw=lw, zorder=10)
        ax.errorbar(profile['mean_r'], profile['mean_b'], yerr=profile['std_b'], c=colors[1], marker='.', label='B-mode', lw=lw, zorder=9)
    else:
        ax.errorbar(profile['mean_r'], profile['sum_w'], yerr=None, c='k', marker='.', label='all', lw=lw, zorder=10)
    
    xlim = ax.get_xlim()
    if coords == "angular":
        ax.plot(xlim, [0,0], 'k:')
        ax.set_xlim(xlim)
    n_pairs = profile['n'].sum()
    print n_pairs
    title = r'$n_\mathrm{pair} = %.2f\cdot 10^9$' % (n_pairs/1e9)
    legend = ax.legend(loc='upper right', numpoints=1, title=title, frameon=False, fontsize='small')
    
def makeSlicedProfile(ax, key_name, profile, plot_type, limits, xlim, ylim, lw=1):
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
        if plot_type == "shear":
            ax.errorbar(profile[s]['mean_r'], profile[s]['mean_e'], yerr=profile[s]['std_e'], c=colors[s], marker='.', label=label, lw=lw)
        else:
            ax.errorbar(profile[s]['mean_r'], profile[s]['sum_w'], yerr=None, c=colors[s], marker='.', label=label, lw=lw)
    ax.set_xlim(xlim)
    #ax.set_ylim(ylim)
    legend = ax.legend(loc='upper right', numpoints=1, title=title, frameon=False, fontsize='small')
    plt.setp(legend.get_title(),fontsize='small')

def makeAxisLabels(ax, coords, plot_type):
    if plot_type == "shear":
        ax.set_ylabel(r'$\Delta\Sigma\ [10^{14}\ \mathrm{M}_\odot \mathrm{Mpc}^{-2}]$')
    else:
        # FIXME: make sure bins are given in Mpc
        # even when physical: they are currently in Mpc/h
        # otherwise they are in arcmin
        ax.set_ylabel(r'$\sum_\mathrm{pairs}{\langle\Sigma_\mathrm{crit}^{-2}\rangle}_w$')
    if coords == "physical":
        ax.set_xlabel('Radius [Mpc/$h$]')
        ax.set_xscale('symlog', linthreshx=1e-2)
        ax.xaxis.set_minor_locator(matplotlib.ticker.LogLocator(subs=np.arange(2, 10)))
        if plot_type == "shear":
            ax.set_yscale('symlog', linthreshy=1e3)
        else:
            ax.set_yscale('log')
        ax.yaxis.set_minor_locator(matplotlib.ticker.LogLocator(subs=np.arange(2, 10)))
    else:
        ax.set_xlabel('Radius [arcmin]')


if __name__ == '__main__':
    # parse inputs
    try:
        configfile = argv[1]
        coords = argv[2]
    except IndexError:
        print "usage: " + argv[0] + " <config file> <angular/physical> [shear/boost]"
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

    try:
        plot_type = argv[3]
    except IndexError:
        plot_type = "shear"

    if plot_type not in ['shear', 'boost']:
        print "specify either shear or boost as plot_type"
        raise SystemExit
    
    indir = os.path.dirname(configfile) + "/"
    outdir = indir

    # load profiles
    profiles = {}
    name = "shear_profile_"
    if plot_type == "boost":
        name = "boost_factor_"
    
    # E/B
    profiles['EB'] = np.load(indir + 'shear_profile_%s_EB.npz' % coords)
    # splits
    for key, limit  in config['splittings'].iteritems():
        profiles[key] = []
        for s in xrange(len(limit)-1):
            filename = indir + 'shear_profile_%s_%s_%d.npz' % (coords, key, s)
            profiles[key].append(np.load(filename))

    # plot generation: E/B profile
    setTeXPlot(sampling=2)
    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(111)
    makeEBProfile(ax, profiles['EB'], plot_type)
    present = profiles['EB']['n'] > 0
    pivot = (profiles['EB']['mean_e'] + profiles['EB']['std_e'])[present].max()
    pivot_ = (profiles['EB']['mean_e'])[present].min()
    if coords == "physical":
        xlim = (1e-2, profiles['EB']['mean_r'][present].max()*1.5)
        ylim = (1e3, 1e7)#max(0, pivot_/1.25), pivot*1.5)
    else:
        xlim = ax.get_xlim()
        ylim = (-0.15*pivot, 1.25*pivot)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    makeAxisLabels(ax, coords, plot_type)
    fig.subplots_adjust(wspace=0, hspace=0, left=0.17, bottom=0.13, right=0.98, top=0.98)
    plotfile = outdir + name + '%s_EB.png' % coords
    fig.savefig(plotfile)

    # sliced profile plots
    for key in config['splittings'].keys():
        print "  " + key
        fig = plt.figure(figsize=(5, 4))
        ax = fig.add_subplot(111)
        makeSlicedProfile(ax, key, profiles[key], plot_type, config['splittings'][key], xlim, ylim)
        makeAxisLabels(ax, coords, plot_type)
        fig.subplots_adjust(wspace=0, hspace=0, left=0.16, bottom=0.13, right=0.98, top=0.97)
        plotfile = outdir + name + "%s_%s.png" % (coords, key)
        fig.savefig(plotfile)

