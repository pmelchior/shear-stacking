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
    else:
        ax.set_xlabel('Radius [arcmin]')

def makeProfileRatio(profile1, profile2):
    p = {}
    for k, v in profile1.iteritems():
        if k in ['mean_e', 'sum_w']:
            try:
                p[k] = v/extrap(profile1['mean_r'], profile2['mean_r'], profile2[k])
            except KeyError:
                pass
        else:
            p[k] = v
    return p
    
if __name__ == '__main__':
    # parse inputs
    try:
        configfile1 = argv[1]
        configfile2 = argv[2]
        coords = argv[3]
    except IndexError:
        print "usage: " + argv[0] + " <config file 1> <config file 2> <angular/physical> [shear/boost]"
        raise SystemExit
    try:
        fp = open(configfile1)
        print "opening configfile " + configfile1
        config1 = json.load(fp)
        fp.close()
    except IOError:
        print "configfile " + configfile1 + " does not exist!"
        raise SystemExit
    try:
        fp = open(configfile2)
        print "opening configfile " + configfile2
        config2 = json.load(fp)
        fp.close()
    except IOError:
        print "configfile " + configfile2 + " does not exist!"
        raise SystemExit
    
    if coords not in ['angular', 'physical']:
        print "specify either angular or physical coordinates"
        raise SystemExit

    try:
        plot_type = argv[4]
    except IndexError:
        plot_type = "shear"

    if plot_type not in ['shear', 'boost']:
        print "specify either shear or boost as plot_type"
        raise SystemExit
    
    indir1 = os.path.dirname(configfile1) + "/"
    indir2 = os.path.dirname(configfile2) + "/"
    outdir = indir1

    # load profiles
    profilesr = {}
    profiles1 = {}
    profiles2 = {}
    name = "shear_profile_ratio_"
    if plot_type == "boost":
        name = "boost_factor_ratio_"
    
    # E/B
    profiles1['EB'] = np.load(indir1 + 'shear_profile_%s_EB.npz' % coords)
    profiles2['EB'] = np.load(indir2 + 'shear_profile_%s_EB.npz' % coords)
    profilesr['EB'] = makeProfileRatio(profiles1['EB'], profiles2['EB'])
        
    # splits
    for key, limit  in config1['splittings'].iteritems():
        profiles1[key] = []
        profiles2[key] = []
        profilesr[key] = []
        for s in xrange(len(limit)-1):
            filename = indir1 + 'shear_profile_%s_%s_%d.npz' % (coords, key, s)
            profiles1[key].append(np.load(filename))
            filename = indir2 + 'shear_profile_%s_%s_%d.npz' % (coords, key, s)
            profiles2[key].append(np.load(filename))
            profilesr[key].append(makeProfileRatio(profiles1[key][-1], profiles2[key][-1]))

    # plot generation: E/B profile
    setTeXPlot(sampling=2)
    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(111)
    makeEBProfile(ax, profilesr['EB'], plot_type)
    if coords == "physical":
        xlim = (1e-2, profilesr['EB']['mean_r'].max()*1.5)
        ylim = (1e3, 1e7)
    else:
        xlim = ax.get_xlim()
        ylim = (-0.15*pivot, 1.25*pivot)
    ax.set_xlim(xlim)
    #ax.set_ylim(ylim)
    makeAxisLabels(ax, coords, plot_type)
    fig.subplots_adjust(wspace=0, hspace=0, left=0.17, bottom=0.13, right=0.98, top=0.98)
    plotfile = outdir + name + '%s_EB.png' % coords
    fig.savefig(plotfile)

    # sliced profile plots
    for key in config1['splittings'].keys():
        print "  " + key
        fig = plt.figure(figsize=(5, 4))
        ax = fig.add_subplot(111)
        makeSlicedProfile(ax, key, profilesr[key], plot_type, config1['splittings'][key], xlim, ylim)
        makeAxisLabels(ax, coords, plot_type)
        fig.subplots_adjust(wspace=0, hspace=0, left=0.16, bottom=0.13, right=0.98, top=0.97)
        plotfile = outdir + name + "%s_%s.png" % (coords, key)
        fig.savefig(plotfile)

