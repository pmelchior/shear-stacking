#!/bin/env python

import os, errno, json
import matplotlib
matplotlib.use('agg')
import pylab as plt
import numpy as np
from sys import argv

# use actual LaTeX to render plot and fonts
from pylab import rcParams
def setTeXPlot(sampling=1):
    params = {
        'backend': 'ps',
        'ps.distiller.res': 6000,
        'axes.labelsize': sampling*9,
        'axes.linewidth' : sampling*0.25,
        'font.size': sampling*8,
        'text.fontsize': sampling*8,
        'legend.fontsize': sampling*8,
        'legend.markerscale' : sampling*0.5,
        'xtick.labelsize': sampling*8,
        'ytick.labelsize': sampling*8,
        'font.family': 'serif',
        'font.serif': 'Times',
        'font.weight': 'medium',
        'text.usetex': 'times',
        'figure.subplot.right' : 0.995,
        'figure.subplot.top' : 0.97,
        'figure.subplot.left' : 0.125,
        'figure.subplot.bottom' : 0.07,
    }
    rcParams.update(params)

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

def makeSlicedProfile(ax, key_name, profile, plot_type, limits, lw=1):
    if config['coords'] == "angular":
        ax.plot(xlim, [0,0], 'k:')
    if plt.matplotlib.rcParams['text.usetex']:
        title = r'\texttt{' + key_name.replace("_", "\_") + '}'
    else:
        title = key_name

    # make the profile for all
    ax.errorbar(profile['all']['mean_r'], profile['all']['mean_q'], yerr=profile['all']['std_q'], c='k', marker='.', label='all', lw=lw)

    # make the profile for each split
    colors = getColors(len(limits))
    for s in range(len(limits)-1):
        pname = key_name + "_%d" % s
        label = '$\in ['
        if isinstance(limits[s], (int, long)):
            label += '%d, ' % limits[s]
        else:
            label += '%.2f, ' % limits[s]
        if isinstance(limits[s+1], (int, long)):
            label += '%d)$' % limits[s+1]
        else:
            label += '%.2f)$' % limits[s+1]
        if plot_type == "shear" or plot_type == "scalar":
            ax.errorbar(profile[pname]['mean_r'], profile[pname]['mean_q'], yerr=profile[pname]['std_q'], c=colors[s], marker='.', label=label, lw=lw)
        else:
            ax.errorbar(profile[pname]['mean_r'], profile[pname]['sum_w'], yerr=None, c=colors[s], marker='.', label=label, lw=lw)

    # xlimits
    xmin = profile['all']['mean_r'].min() / 2
    xmax = profile['all']['mean_r'].max() * 2
    ax.set_xlim(xmin, xmax)

    # avoid negative values in shear plots
    if plot_type == "shear":
        ymin, ymax = (1e3, 1e7)
        ax.set_ylim(ymin, ymax)
    
    # decorations
    n_pair = profile['all']['n'].sum()
    ax.text(0.05, 0.95, r'$n_\mathrm{pair} = %.2f\cdot 10^9$' % (n_pair/1e9), ha='left', va='top', transform=ax.transAxes, fontsize='small')
    legend = ax.legend(loc='upper right', numpoints=1, title=title, frameon=False, fontsize='small')
    plt.setp(legend.get_title(),fontsize='small')

def makeAxisLabels(ax, plot_type, config):
    if plot_type == "shear":
        ax.set_ylabel(r'$\Delta\Sigma\ [10^{14}\ \mathrm{M}_\odot \mathrm{Mpc}^{-2}]$')
    if plot_type == "boost":
        ax.set_ylabel(r'$\sum_\mathrm{pairs}{\langle\Sigma_\mathrm{crit}^{-2}\rangle}_w$')
    if plot_type == "scalar":
        if plt.matplotlib.rcParams['text.usetex']:
            ax.set_ylabel(r'\texttt{' + config['shape_scalar_key'].replace("_", "\_") + '}')
        else:
            ax.set_ylabel(config['shape_scalar_key'])
        
    if config['coords'] == "physical":
        ax.set_xlabel('Radius [Mpc/$h$]')
        ax.set_xscale('symlog', linthreshx=1e-2)
        ax.xaxis.set_minor_locator(matplotlib.ticker.LogLocator(subs=np.arange(2, 10)))
        if plot_type == "shear":
            ax.set_yscale('symlog', linthreshy=1e3)
            ax.yaxis.set_minor_locator(matplotlib.ticker.LogLocator(subs=np.arange(2, 10)))
        if plot_type == "boost":
            ax.set_yscale('log')
            ax.yaxis.set_minor_locator(matplotlib.ticker.LogLocator(subs=np.arange(2, 10)))
    else:
        ax.set_xlabel('Radius [arcmin]')


if __name__ == '__main__':
    # parse inputs
    try:
        configfile = argv[1]
        plot_type = argv[2]
    except IndexError:
        print "usage: " + argv[0] + " <config file> <shear/boost/scalar>"
        raise SystemExit
    try:
        fp = open(configfile)
        print "opening configfile " + configfile
        config = json.load(fp)
        fp.close()
    except IOError:
        print "configfile " + configfile + " does not exist!"
        raise SystemExit
    
    if plot_type not in ['shear', 'boost', 'scalar']:
        print "specify plot_type from ['shear', 'boost', 'scalar']"
        raise SystemExit
    
    indir = os.path.dirname(configfile) + "/"
    outdir = indir

    # load profiles
    name = "shear_"
    if plot_type == "boost":
        name = "boost_"
    if plot_type == "scalar":
        name = "scalar_" + config['shape_scalar_key'] + "_"  

    pnames = ['all']
    for key, limit in config['splittings'].iteritems():
        for s in xrange(len(limit)-1):
            pnames.append("%s_%d" % (key, s))
    profiles = {}
    for pname in pnames:
        profiles[pname] = np.load(outdir + name + pname + ".npz")
        
    # plot generation
    setTeXPlot(sampling=2)
    
    # new plot for each slice
    for key in config['splittings'].keys():
        print "  " + key
        fig = plt.figure(figsize=(5, 4))
        ax = fig.add_subplot(111)
        makeSlicedProfile(ax, key, profiles, plot_type, config['splittings'][key])
        makeAxisLabels(ax, plot_type, config)
        fig.subplots_adjust(wspace=0, hspace=0, left=0.16, bottom=0.15, right=0.98, top=0.95)
        plotfile = outdir + name + "%s.png" % key
        fig.savefig(plotfile)

