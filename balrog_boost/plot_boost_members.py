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

def makeAxisLabels(ax, r, coords, plot_type, stacked=False):
    if coords == "physical":
        xlim = (r.min()/1.5, r.max()*1.5)
        if not stacked:
            ax.set_xlabel('Radius [Mpc/$h$]')
        ax.set_xscale('log')
    else:
        xlim = ax.get_xlim()
        if not stacked:
            ax.set_xlabel('Radius [arcmin]')
    ax.set_xlim(xlim)
    #ax.set_ylim(ylim)
    
    if plot_type == "shear":
        ax.set_ylabel(r'$\Delta\Sigma\ [10^{14}\ \mathrm{M}_\odot \mathrm{Mpc}^{-2}]$')
        ax.set_ylim(1e3,1e7)
        ax.set_yscale('symlog', linthreshx=1e2)
    if plot_type == "weight":
        # FIXME: make sure bins are given in Mpc
        # even when physical: they are currently in Mpc/h
        # otherwise they are in arcmin
        ax.set_ylabel(r'$\sum_\mathrm{pairs}{\langle\Sigma_\mathrm{crit}^{-2}\rangle}$')
        ax.xaxis.set_major_formatter(NullFormatter())
        ylim = ax.get_ylim()
        if max(ylim) < 1e-4:
            ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: ('%.1f')%(x*1e5)))
        else:
            ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: ('%d')%(x*1e5)))
    if plot_type == "boost":
        ax.set_ylabel('boost')
        ax.set_ylim(0.95, 1.55)
        ax.plot(xlim, [1,1], 'k:', zorder=1)

        
def makeShearProfile(ax, profile, coords, lw=1):
    colors = getColors(3)
    ax.errorbar(profile['mean_r'], profile['mean_e'], yerr=profile['std_e'], c='k', marker='.', label='E-mode', lw=lw, zorder=10)
    ax.errorbar(profile['mean_r'], profile['mean_b'], yerr=profile['std_b'], c=colors[1], marker='.', label='B-mode', lw=lw, zorder=9)
    xlim = ax.get_xlim()
    if coords == "angular":
        ax.plot(xlim, [0,0], 'k:')
        ax.set_xlim(xlim)
    n_pairs = profile['n'].sum()
    print n_pairs
    title = r'$n_\mathrm{pair} = %.2f\cdot 10^9$' % (n_pairs/1e9)
    legend = ax.legend(loc='upper right', numpoints=1, title=title, frameon=False, fontsize='small')
    makeAxisLabels(ax, profile['mean_r'], coords, 'shear')

def makeWeightProfile(fig, p_r, p_all, p_wo_member, p_balrog, coords, lw=1):
    # comparison of the sum-of-weights
    ax = fig.add_subplot(311)
    ax.errorbar(p_all['mean_r'], p_all['sum_w'], yerr=None, c='k', marker='.', label=r'$w$', lw=lw, zorder=10)
    ax.errorbar(p_wo_member['mean_r'], p_wo_member['sum_w'], yerr=None, c='b', marker='.', label=r'$w_{\backslash u}$', lw=lw, zorder=10)
    ax.errorbar(p_balrog['mean_r'], p_balrog['sum_w'] * p_r['sum_w'][-3:].mean() / p_balrog['sum_w'][-3:].mean(), yerr=None, c='m', marker='.', label=r'$w_b$', lw=lw, zorder=10)

    legend = ax.legend(loc='upper right', numpoints=1, frameon=False, fontsize='small')
    makeAxisLabels(ax, p_r['mean_r'], coords, 'weight', stacked=True)

    # member weights
    ax = fig.add_subplot(312)
    ax.errorbar(p_wo_member['mean_r'], p_all['sum_w'] - p_wo_member['sum_w'], yerr=None, c='r', marker='.', label=r'$w_u$', lw=lw)
    legend = ax.legend(loc='upper right', numpoints=1, frameon=False, fontsize='small')
    makeAxisLabels(ax, p_r['mean_r'], coords, 'weight', stacked=True)

    # boost factors
    ax = fig.add_subplot(313)
    ax.errorbar(p_r['mean_r'], p_r['boost_all'], yerr=p_r['std_boost_all'], c='k', marker='.', label=r'$b$', lw=lw, zorder=10)
    ax.errorbar(p_r['mean_r'], p_r['boost'], yerr=p_r['std_boost'], c='b', marker='.', label=r'$b_{\backslash u}$', lw=lw, zorder=10)
    legend = ax.legend(loc='upper right', numpoints=1, frameon=False, fontsize='small')
    makeAxisLabels(ax, p_r['mean_r'], coords, 'boost')

    
def makeSlicedProfile(ax, key_name, profile, plot_type, limits, lw=1):
    if plt.matplotlib.rcParams['text.usetex']:
        title = r'\texttt{' + key_name.replace("_", "\_") + '}'
    else:
        title = key_name

    # make the profile for each split
    colors = getColors(len(limits))
    for s in range(len(limits)-1):
        if profile[s]['n'].sum() > 0:
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
                try:
                    ax.errorbar(profile[s]['mean_r'], profile[s]['mean_e'], yerr=profile[s]['std_e'], c=colors[s], marker='.', label=label, lw=lw)
                except AssertionError:
                    print key_name, s
                    print profile[s]['mean_r']
                    print profile[s]['mean_e']
                    print profile[s]['std_e']
            if plot_type == "weight":
                ax.errorbar(profile[s]['mean_r'], profile[s]['sum_w'], yerr=None, c=colors[s], marker='.', label=label, lw=lw)
            if plot_type == "boost":
                ax.errorbar(profile[s]['mean_r'], profile[s]['boost'], yerr=profile[s]['std_boost'], c=colors[s], marker='.', label=label, lw=lw)
    legend = ax.legend(loc='upper right', numpoints=1, title=title, frameon=False, fontsize='small')
    plt.setp(legend.get_title(),fontsize='small')
    makeAxisLabels(ax, profile[0]['mean_r'], coords, plot_type)


def makeBoostedProfile(all, wo_member, balrog, common_r=None):
    p = {}
    if common_r is None:
        r = wo_member['mean_r']
    else:
        r = common_r
    if all['mean_r'].size:
        w = extrap(r, all['mean_r'], all['sum_w'])
        w_o = extrap(r, wo_member['mean_r'], wo_member['sum_w'])
        w_b = extrap(r, balrog['mean_r'], balrog['sum_w'])

        p['sum_w'] = w_b * w_o[-4:-1].mean() / w_b[-4:-1].mean()
        p['boost'] = w_o / (w_b * w_o[-4:-1].mean() / w_b[-4:-1].mean())
        p['boost_all'] = w / (w_b * w[-4:-1].mean() / w_b[-4:-1].mean())
        p['mean_e'] = extrap(r, wo_member['mean_r'], wo_member['mean_e']) * p['boost']
        p['std_e'] = extrap(r, wo_member['mean_r'], wo_member['std_e']) * p['boost']
        try:
            p['mean_b'] = extrap(r, wo_member['mean_r'], wo_member['mean_b']) * p['boost']
            p['std_b'] = extrap(r, wo_member['mean_r'], wo_member['std_b']) * p['boost']
        except KeyError:
            pass
        p['mean_r'] = r
        p['n'] = extrap(r, wo_member['mean_r'], wo_member['n'])
    else:
        p['mean_r'] = r
        p['n'] = np.zeros_like(r, dtype='i1')
        p['sum_w'] = np.zeros_like(r)
        p['mean_e'] = np.zeros_like(r)
        p['std_e'] = np.zeros_like(r)
    return p

def computeMeanStdForJackknifeProfile(profile, quantity, name, s=None, n_jack=0):
    q = []
    missing = []
    for i in xrange(1, n_jack + 1):
        if s is None:
            p = profile[i][name]
        else:
            p = profile[i][name][s]
        n_, q_ = p['n'], p[quantity]
        missing.append(n_ == 0)
        q.append(q_)

    missing = np.array(missing)
    q = np.ma.masked_array(q, mask=missing)
    mean_q = q.mean(axis=0)

    # result for normal/non-jackknife profile
    if s is None:
        p = profile[0][name]
    else:
        p = profile[0][name][s]
    n, mean0 =  p['n'], p[quantity]
        
    # variance and bias-corrected mean needs number of actual jackknifes:
    # to be corrected for available data in each radial bin
    n_avail = n_jack - missing.sum(axis=0)
    mean_q = n_avail*mean0 - (n_avail - 1)*mean_q
    std_q = ((n_avail - 1.)/n_avail * ((q - mean_q)**2).sum(axis=0))**0.5
    return mean_q.data, std_q.data

def mergeJackknifeProfiles(profile, key, s=None, n_jack=0):
    if n_jack:
        p = profile[0][key]
    else:
        p = profile[key]
    if s is None:
        p_ = p
    else:
        p_ = p[s]
    has_boost = 'boost' in p_.keys()
    p_['mean_e'], p_['std_e'] = computeMeanStdForJackknifeProfile(profile, 'mean_e', key, s=s, n_jack=n_jack)
    if has_boost:
        p_['boost'], p_['std_boost'] = computeMeanStdForJackknifeProfile(profile, 'boost', key, s=s, n_jack=n_jack)
        p_['boost_all'], p_['std_boost_all'] = computeMeanStdForJackknifeProfile(profile, 'boost_all', key, s=s, n_jack=n_jack)
    return p_

def mergeAllSubProfiles(profile, config, n_jack=0):
    p = {}
    p['EB'] = mergeJackknifeProfiles(profile, 'EB', n_jack=n_jack)
    for key, limit  in config['splittings'].iteritems():
        p[key] = []
        for s in xrange(len(limit)-1):
            p[key].append(mergeJackknifeProfiles(profile, key, s=s, n_jack=n_jack))
    return p

# all profiles should have 19 elements
# some have fewer because the bins were no filled
def makeNineteen(x):
    try:
        it = iter(x)
        if x.size < 19:
            x = np.concatenate((np.zeros(19-x.size, dtype=x.dtype), x))
    except TypeError:
        pass
    return x

def load_close(filename):
    f = np.load(filename)
    d = {}
    for k,v in f.iteritems():
        d[k] = makeNineteen(v)
    f.close()
    return d

def loadProfiles(indir1, indir2, indir3, config, coords, n_jack=0, common_r=None, set_r=False):
    if n_jack == 0:
        # load profiles
        p_r = {}
        p_all = {}
        p_wo_member = {}
        p_balrog = {}

        # E/B
        p_all['EB'] = load_close(indir1 + 'shear_profile_%s_EB.npz' % coords)
        p_wo_member['EB'] = load_close(indir2 + 'shear_profile_%s_EB.npz' % coords)
        if set_r:
            common_r = p_wo_member['EB']['mean_r']
        p_balrog['EB'] = load_close(indir3 + 'shear_profile_%s_EB.npz' % coords)
        p_r['EB'] = makeBoostedProfile(p_all['EB'], p_wo_member['EB'], p_balrog['EB'], common_r=common_r)

        # splits
        for key, limit  in config1['splittings'].iteritems():
            p_all[key] = []
            p_wo_member[key] = []
            p_balrog[key] = []
            p_r[key] = []
            for s in xrange(len(limit)-1):
                filename = indir1 + 'shear_profile_%s_%s_%d.npz' % (coords, key, s)
                p_all[key].append(load_close(filename))
                filename = indir2 + 'shear_profile_%s_%s_%d.npz' % (coords, key, s)
                p_wo_member[key].append(load_close(filename))
                filename = indir3 + 'shear_profile_%s_%s_%d.npz' % (coords, key, s)
                p_balrog[key].append(load_close(filename))
                p_r[key].append(makeBoostedProfile(p_all[key][-1], p_wo_member[key][-1], p_balrog[key][-1], common_r=common_r))
        return p_r, p_all, p_wo_member, p_balrog
    else:
        lp_r = []
        lp_all = []
        lp_wo_member = []
        lp_balrog = []
        subdir = "n_jack_%d/" % n_jack
        
        # load non-jackknife profile for bias correction
        p_r_, p_all_, p_wo_member_, p_balrog_ = loadProfiles(indir1, indir2, indir3, config, coords, n_jack=0, set_r=True)
        if set_r:
            common_r = p_r_['EB']['mean_r']
        lp_r.append(p_r_)
        lp_all.append(p_all_)
        lp_wo_member.append(p_wo_member_)
        lp_balrog.append(p_balrog_)

        for i in xrange(n_jack):
            jack_nr = "_jack_%d" % i
            p_r_, p_all_, p_wo_member_, p_balrog_ = loadProfiles(indir1 + subdir, indir2 + subdir, indir3 + subdir, config, coords + jack_nr, n_jack=0, common_r=common_r)
            lp_r.append(p_r_)
            lp_all.append(p_all_)
            lp_wo_member.append(p_wo_member_)
            lp_balrog.append(p_balrog_)

        p_r = mergeAllSubProfiles(lp_r, config, n_jack=n_jack)
        p_all = mergeAllSubProfiles(lp_all, config, n_jack=n_jack)
        p_wo_member = mergeAllSubProfiles(lp_wo_member, config, n_jack=n_jack)
        #p_balrog = mergeAllSubProfiles(lp_balrog, config, n_jack=n_jack)
        p_balrog = lp_balrog[0]
        return p_r, p_all, p_wo_member, p_balrog

    
if __name__ == '__main__':
    # parse inputs
    try:
        configfile1 = argv[1]  
        configfile2 = argv[2]
        configfile3 = argv[3]
        coords = argv[4]
    except IndexError:
        print "usage: " + argv[0] + " <default file> <w/o members file> <balrog file> <angular/physical>"
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
    try:
        fp = open(configfile3)
        print "opening configfile " + configfile3
        config3 = json.load(fp)
        fp.close()
    except IOError:
        print "configfile " + configfile3 + " does not exist!"
        raise SystemExit
    
    if coords not in ['angular', 'physical']:
        print "specify either angular or physical coordinates"
        raise SystemExit
    
    indir1 = os.path.dirname(configfile1) + "/"
    indir2 = os.path.dirname(configfile2) + "/"
    indir3 = os.path.dirname(configfile3) + "/"
    outdir = indir2
    label = os.path.dirname(configfile2).split("/")[-1]

    # load profiles
    n_jack = 40
    p_r, p_all, p_wo_member, p_balrog = loadProfiles(indir1, indir2, indir3, config1, coords, n_jack=n_jack, set_r=True)

    
    # plot generation: E/B profile
    name = {'shear': outdir + 'shear_profile_boosted_' + label + '_%s_%s.png',
            'boost': outdir + 'boost_factor_corrected_' + label + '_%s_%s.png' }
    setTeXPlot(sampling=2)

    """
    key = 'EB'
    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(111)
    makeShearProfile(ax, p_r[key], coords)
    fig.subplots_adjust(wspace=0, hspace=0, left=0.17, bottom=0.13, right=0.98, top=0.98)
    #plotfile = name['shear'] % (coords, key)
    #fig.show()
    #fig.savefig(plotfile)
    """
    
    # weight / boost profile
    key = 'EB'
    fig = plt.figure(figsize=(5, 7))
    makeWeightProfile(fig, p_r[key], p_all[key], p_wo_member[key], p_balrog[key], coords)
    fig.subplots_adjust(wspace=0, hspace=0.15, left=0.17, bottom=0.09, right=0.98, top=0.98)
    plotfile = name['boost'] % (coords, key)
    fig.show()
    fig.savefig(plotfile)
    
    # sliced profile plots
    for key in config1['splittings'].keys()[-1:]:
        print "  " + key
        for plot_type in ['shear', 'boost']:
            fig = plt.figure(figsize=(5, 4))
            ax = fig.add_subplot(111)
            makeSlicedProfile(ax, key, p_r[key], plot_type, config1['splittings'][key])
            fig.subplots_adjust(wspace=0, hspace=0, left=0.16, bottom=0.15, right=0.98, top=0.97)
            plotfile = name[plot_type] % (coords, key)
            fig.savefig(plotfile)

            if plot_type == 'boost':
                fig.show()
                
            for s in xrange(len(p_r[key])):
                np.savez(plotfile.replace(".png", "_%d.npz" % s), **p_r[key][s])

            """
            # Show weight profiles for each slice
            for s in xrange(len(p_r[key])):
                fig = plt.figure(figsize=(5, 7))
                makeWeightProfile(fig, p_r[key][s], p_all[key][s], p_wo_member[key][s], p_balrog[key][s], coords)
                fig.subplots_adjust(wspace=0, hspace=0.15, left=0.17, bottom=0.09, right=0.98, top=0.98)
                fig.show()
            """
