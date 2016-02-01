#!/bin/env python

import os, errno, json, fitsio, copy
#import matplotlib
#matplotlib.use('agg')
import pylab as plt
import numpy as np
from sys import argv, path
path.append(os.path.abspath(os.path.join(os.getcwd(), os.pardir)))
from shear_stacking import *
from glob import glob
from matplotlib.ticker import NullFormatter, FuncFormatter


def makeWeightProfile(fig, p, p_balrog, p_r, config, lw=1):
    # comparison of the sum-of-weights
    ax = fig.add_subplot(211)
    ax.errorbar(p['mean_r'], p['sum_w'], yerr=None, c='k', marker='.', label=r'$w$', lw=lw, zorder=10)
    ax.errorbar(p_balrog['mean_r'], p_balrog['sum_w'] * p['sum_w'][-4:].mean() / p_balrog['sum_w'][-4:].mean(), yerr=None, c='m', marker='.', label=r'$w_b$', lw=lw, zorder=10)

    legend = ax.legend(loc='upper right', numpoints=1, frameon=False, fontsize='small')
    makeAxisLabels(ax, 'weight', config, stacked=True)

    # xlimits
    xmin = p['mean_r'].min() / 2
    xmax = p['mean_r'].max() * 2
    ax.set_xlim(xmin, xmax)

    ymin = p['sum_w'].min() / 1.2
    ymax = p['sum_w'].max() * 1.2
    mag = int(np.floor(np.log10(ymax)))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: ('%.1f')%(x/(10**mag))))
    ax.set_ylim(ymin, ymax)
              
    
    # boost factors
    ax = fig.add_subplot(212)
    ax.errorbar(p_r['mean_r'], p_r['mean_boost'], yerr=p_r['std_boost'], c='b', marker='.', label=r'$b$', lw=lw, zorder=10)
    makeAxisLabels(ax, 'boost', config)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(0.95, 1.55)
    ax.plot([xmin, xmax], [1,1], 'k:', zorder=1)

def makeSlicedProfile(ax, key_name, profile, plot_type, config, lw=1):
    if plt.matplotlib.rcParams['text.usetex']:
        title = r'\texttt{' + key_name.replace("_", "\_") + '}'
    else:
        title = key_name

    limits = config['splittings'][key]
        
    # make the profile for all
    if plot_type == "shear":
        ax.errorbar(profile['all']['mean_r'], profile['all']['mean_q'], yerr=profile['all']['std_q'], c='k', marker='.', label='all', lw=lw)
    if plot_type == "weight":
        ax.errorbar(profile['all']['mean_r'], profile['all']['sum_w'], yerr=None, c='k', marker='.', label='all', lw=lw)
    if plot_type == "boost":
        ax.errorbar(profile['all']['mean_r'], profile['all']['mean_boost'], yerr=profile['all']['std_boost'], c='k', marker='.', label='all', lw=lw)
        
    # make the profile for each split
    colors = getColors(len(limits))
    for s in range(len(limits)-1):
        pname = key_name + "_%d" % s
        if profile[pname]['n'].sum() > 0:
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
                ax.errorbar(profile[pname]['mean_r'] * (1 + 0.05*(s-0.5)), profile[pname]['mean_q'], yerr=profile[pname]['std_q'], c=colors[s], marker='.', label=label, lw=lw)
            if plot_type == "weight":
                ax.errorbar(profile[pname]['mean_r'] * (1 + 0.05*(s-0.5)), profile[pname]['sum_w'], yerr=None, c=colors[s], marker='.', label=label, lw=lw)
            if plot_type == "boost":
                ax.errorbar(profile[pname]['mean_r'] * (1 + 0.05*(s-0.5)), profile[pname]['mean_boost'], yerr=profile[pname]['std_boost'], c=colors[s], marker='.', label=label, lw=lw)
    legend = ax.legend(loc='upper right', numpoints=1, title=title, frameon=False, fontsize='small')
    plt.setp(legend.get_title(),fontsize='small')
    makeAxisLabels(ax, plot_type, config)

    # xlimits
    xmin = p['all']['mean_r'].min() / 2
    xmax = p['all']['mean_r'].max() * 2
    ax.set_xlim(xmin, xmax)

    if plot_type == "shear":
        ax.set_ylim(1e3,1e7)
    if plot_type == "boost":
        ax.set_ylim(0.95, 1.55)
        ax.plot([xmin, xmax], [1,1], 'k:', zorder=1)
    


def makeBoostedProfile(shapes, balrog, common_r=None):
    p = {}
    if common_r is None:
        r = shapes['all']['mean_r']
    else:
        r = common_r

    for pname in shapes.keys():
        p[pname] = {}
        w = extrap(r, shapes[pname]['mean_r'], shapes[pname]['sum_w'])
        w_b = extrap(r, balrog[pname]['mean_r'], balrog[pname]['sum_w'])

        p[pname]['mean_r'] = r
        p[pname]['n'] = shapes[pname]['n']
        p[pname]['sum_w'] = w_b * w[-4:].mean() / w_b[-4:].mean()
        p[pname]['mean_boost'] = w / p[pname]['sum_w']
        # we have no means to estimate error of the boost from the aggregate profiles
        p[pname]['std_boost'] = np.zeros_like(p[pname]['mean_boost'])
        try:
            p[pname]['mean_q'] = extrap(r, shapes[pname]['mean_r'], shapes[pname]['mean_q']) * p[pname]['mean_boost']
        except KeyError:
            p[pname]['mean_q'] = extrap(r, shapes[pname]['mean_r'], shapes[pname]['q']) * p[pname]['mean_boost']
        p[pname]['std_q'] = extrap(r, shapes[pname]['mean_r'], shapes[pname]['std_q']) * p[pname]['mean_boost']
    return p

def computeJackknife(profile, key):
    n_jack = len(profile)-1
    for pname in profile[0].keys():
        q = []
        missing = []
        for i in xrange(n_jack):
            missing.append(profile[i+1][pname]['n'] == 0)
            q.append(profile[i+1][pname]['mean_' + key])
        missing = np.array(missing)
        q = np.ma.masked_array(q, mask=missing)
        mean_q = q.mean(axis=0)

        # aggregated/non-jackknife profile
        mask = (profile[0][pname]['n'] > 0)
        mean0 = profile[0][pname]['mean_' + key]

        # variance and bias-corrected mean needs number of actual jackknifes:
        # to be corrected for available data in each radial bin
        n_avail = n_jack - missing.sum(axis=0)
        mean_q = n_avail*mean0 - (n_avail - 1)*mean_q
        std_q = ((n_avail - 1.)/n_avail * ((q - mean_q)**2).sum(axis=0))**0.5

        # update aggregated profile
        profile[0][pname]['mean_' + key] = mean_q.data
        profile[0][pname]['std_' + key] = std_q.data
    
def loadProfiles(indir1, indir2, config, n_jack=0):
    # define keys
    file_name = "shear_"
    pnames = ['all']
    for key, limit in config['splittings'].iteritems():
        for s in xrange(len(limit)-1):
            pnames.append("%s_%d" % (key, s))

    # apply boost factors to aggregated profiles
    if n_jack == 0:
        # load profiles
        p = {}
        p_balrog = {}
        for pname in pnames:
            p[pname] = np.load(indir1 + file_name + pname + ".npz")
            p_balrog[pname] = np.load(indir2 + file_name + pname + ".npz")

        p_r = makeBoostedProfile(p, p_balrog)
        return p, p_balrog, p_r
    else:
        p, p_balrog, p_r = loadProfiles(indir1, indir2, config, n_jack=0)
        common_r = p_r['all']['mean_r']
        p_r_i = []
        p_r_i.append(p_r)
        for i in xrange(n_jack):
            p = {}
            p_balrog = {}
            for pname in pnames:
                p[pname] = np.load(indir1 + "n_jack/" + file_name + pname + ("_%d" % i) +".npz")
                p_balrog[pname] = np.load(indir2 + "n_jack/" + file_name + pname + ("_%d" % i) + ".npz")

            p_r_i.append(makeBoostedProfile(p, p_balrog, common_r=common_r))
        computeJackknife(p_r_i, 'q') 
        computeJackknife(p_r_i, 'boost')
        return p, p_balrog, p_r_i[0]

    
if __name__ == '__main__':
    # parse inputs
    try:
        configfile1 = argv[1]  
        configfile2 = argv[2]
    except IndexError:
        print "usage: " + argv[0] + " <lensing config> <balrog config>"
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

    if config1['coords'] not in ['angular', 'physical']:
        print "specify either angular or physical coordinates"
        raise SystemExit
    
    indir1 = os.path.dirname(configfile1) + "/"
    indir2 = os.path.dirname(configfile2) + "/"
    outdir = indir1
    label = os.path.dirname(configfile1).split("/")[-1]

    # load profiles
    p, p_balrog, p_r = loadProfiles(indir1, indir2, config1, n_jack=config1['n_jack'])
    
    # plot generation: E/B profile
    name = {'shear': outdir + 'shear_boosted_' + label + '_%s.png',
            'boost': outdir + 'boost_factor_' + label + '_%s.png' }
    setTeXPlot(sampling=2)

    # weight / boost profile
    key = 'all'
    fig = plt.figure(figsize=(5, 6))
    makeWeightProfile(fig, p[key], p_balrog[key], p_r[key], config1)
    fig.subplots_adjust(wspace=0, hspace=0.14, left=0.17, bottom=0.10, right=0.98, top=0.97)
    plotfile = name['boost'] % key
    fig.show()
    # fig.savefig(plotfile)


    # sliced profile plots
    for key,limits in config1['splittings'].iteritems():
        print "  " + key
        for plot_type in ['boost', 'shear']:
            fig = plt.figure(figsize=(5, 4))
            ax = fig.add_subplot(111)
            makeSlicedProfile(ax, key, p_r, plot_type, config1)
            fig.subplots_adjust(wspace=0, hspace=0, left=0.16, bottom=0.15, right=0.98, top=0.97)
            plotfile = name[plot_type] % key
            fig.savefig(plotfile)
            fig.show()

        # save shear and boost factors to new npz
        for s in xrange(len(limits)-1):
            pname = key + "_%d" % s
            filename = (name['shear'] % pname).replace(".png", ".npz")
            np.savez(filename, **p_r[pname])


