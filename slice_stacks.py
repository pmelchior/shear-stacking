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

import matplotlib.colors as colors
import matplotlib.cm as cmx
def makeSlicedProfilePlot(ax, radius, DeltaSigma, weight, slices, splittings, key_name, lw=1):
    if plt.matplotlib.rcParams['text.usetex']:
        label = r'\texttt{' + key_name.replace("_", "\_") + '} all'
    else:
        label = key_name + " all"
    mean_r, n, mean_q, std_q = scalarProfile(bins, radius, DeltaSigma, weight)
    ax.errorbar(mean_r, mean_q, yerr=std_q, c='k', marker='.', label=label, lw=lw)

    # use a diverging colormap for the colors of lines
    my_map = cm = plt.get_cmap('RdYlGn') 
    cNorm  = colors.Normalize(vmin=0, vmax=(len(splittings)-2))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=my_map)

    # make the profile for each split
    for s in range(len(splittings)-1):
        label = '$\in[%.2f, %.2f)$' % (splittings[s], splittings[s+1])
        if s == 0:
            label = '$<%.2f$' % splittings[s+1]
        if s == len(splittings)-2:
            label = '$\geq %.2f$' % splittings[s]
        if slices[s].size:
            mean_r, n, mean_q, std_q = scalarProfile(bins, radius[slices[s]], DeltaSigma[slices[s]], weight[slices[s]])
            ax.errorbar(mean_r, mean_q, yerr=std_q, c=scalarMap.to_rgba(s), marker='.', label=label, lw=lw)

def sizeRatio(data):
    return 2*data['RADIUS'] / (data['PSF_FWHM'] * 0.263)

def SNRadius(data):
    return data['SNR'] * data['RADIUS']

from struct import unpack
class HTMFile:
    """Class to read in HTM match files sequentially
    
    Provides two convenient iterators:
      htmf = HTMFile(filename)
      for m1, m2, d12 in htmf:
          # do somthing with a single matched m1, m2
      for m1, m2s, d12s in htmf.matches():
          # do something with the list of matches m2s of a single m1
    """
    def __init__(self, filename):
        self.fp = open(filename, 'rb')
        self.n_matches = unpack('q', self.fp.read(8))[0]
        self.m1_current = -1
    def __iter__(self):
        return self
    def next(self):
        """Line iterator.

        Returns one match of m1 and m2 with the relative distance d12 (in deg).
        """
        line = self.fp.read(24)
        if line != '':
            return unpack('qqd', line)
        else:
            raise StopIteration
    def matches(self):
        """Match iterator.
        
        Returns the current match index m1, the list of matches m2 and their
        respective distances (in deg).
        """
        while self.fp.tell() < self.n_matches * 24:
            m1, m2, d12 = self.next()
            self.m1_current = m1
            m2s = [m2]
            d12s = [d12]
            while True:
                try:
                    m1, m2, d12 = self.next()
                    if m1 == self.m1_current:
                        m2s.append(m2)
                        d12s.append(d12)
                    else: # if next m1: rewind to previous line
                        self.fp.seek(-24, 1)
                        break
                except StopIteration: # at end of file, return current set
                    break
            yield self.m1_current, m2s, d12s
    def __del__(self):
        self.fp.close()


if __name__ == '__main__':
    if len(argv) < 2:
        print "usage: " + argv[0] + " <band>"
        exit(1)

    band = argv[1]
    lensfile = 'data/sva1_gold_1.0_run_redmapper_v5.9_lgt20_catalog.fit'
    lens_z_key = 'Z_LAMBDA'
    shapefile = 'data/011-im3shape-3_' + band.lower() + '_photoz_gold_Feb12.fits'
    matchfile = '/tmp/matches_' + band.lower() + '.bin'
    stackfile = '/tmp/shear_stack_' + band.lower() + '.npz'
    plotfile = 'shear_stack_' + band.lower() + '.pdf'
    keys = ['SNR', 'RADIUS', 'PSF_FWHM', sizeRatio, 'RGPP_RP', SNRadius, 'ZP_2', 'FLAGS_' + band.upper()]
    # slices are either equal-volumne or pre-determined
    split = 2 # equal-volume slices
    splittings = [split, split, split, split, split, split, split, [0,1,10]]

    # open lens catalog, apply selection if desired
    hdu = pyfits.open(lensfile)
    lenses = hdu[1].data
    #good_cl = (lenses[lens_z_key] < 0.8) & (lenses[lens_z_key] > 0.1)
    #lenses = lenses[good_cl]
    print "lens sample: %d" % lenses.size

    # open shapes, apply post-run selections
    shdu = pyfits.open(shapefile)
    good_sh = ((shdu[1].data['FLAG'] & (1+4+8)) == 0) & invert(isnan(shdu[1].data['SNR'])) & (shdu[1].data['SNR'] > 5) & (shdu[1].data['FLAGS_I'] == 0) & (shdu[1].data['FLAGS_R'] == 0)
    shapes = shdu[1].data[good_sh]

    # find all galaxies in shape catalog within maxranfe arcmin 
    # of each lens center
    # NOTE: maxrange can be changed for each cluster when working in Mpc instead of arcmin
    print "matching lens and source catalog..."
    maxrange = 30. # arcmin
    if exists(matchfile) is False:
        # CAVEAT: make sure to have enough space where you put the match file
        # it has 24 byte per match, which quickly becomes Gb's of data 
        h = eu.htm.HTM(8)
        h.match(lenses['RA'], lenses['DEC'], shapes['RA'], shapes['DEC'], maxrange/60, maxmatch=-1, file=matchfile)
        del h
    else:
        print " re-using existing matchfile", matchfile
 
    htmf = HTMFile(matchfile)
    print " found ", htmf.n_matches, "matches"

    if exists(stackfile) is False:
        # set up the container for the shears
        Ngal = htmf.n_matches
        DeltaSigma = np.empty(Ngal, dtype='float32')
        radius = np.empty(Ngal, dtype='float32')
        weight = np.empty(Ngal, dtype='float32')

        # define the slices and set up their index containers
        # can be any key in the shape catalog or a function thereof
        slices = {}
        print "determining slice ranges..."
        for i in range(len(keys)):
            if callable(keys[i]):
                key_name = keys[i].__name__
            else:
                key_name = keys[i]
            slices[key_name] = {}
            # FIXME: there is no selection on the shapes here, e.g. no BG cut
            values = getValues(shapes, keys[i])
            # if the slices are set already: find out how many source will fall in each
            if iterable(splittings[i]):
                for s in range(len(splittings[i])-1):
                    Ngal_slice = getSliceMask(values, splittings[i][s], splittings[i][s+1], return_num=True)*(1.*Ngal)/len(shapes)
                    slices[key_name][s] = np.zeros(Ngal_slice, dtype='int64')
            # for percentile splits: create equal-size containers and determine limits
            else:
                for s in range(splittings[i]):
                    # CAVEAT: this is not exactly correct, but deviations
                    # should be minimal for large stacks
                    # -> Check for the size of the contained before filling it!
                    slices[key_name][s] = np.zeros(Ngal/splittings[i], dtype='int64')
                delta = 100./splittings[i]
                ranges = [k*delta for k in range(splittings[i]+1)]
                splittings[i] = percentile(values, ranges)
                del ranges
            del values

        # remember the last element (of the whole array DeltaSigma) for each
        # of the slices
        last_element = {'all': 0}
        for k,v in slices.iteritems():
            last_element[k] = {}
            for vv in v:
                last_element[k][vv] = 0

        # iterate over all lenses
        print "stacking lenses..."
        specz_calib = getSpecZCalibration()
        counter = 1
        for m1, m2, d12 in htmf.matches():
            lens = lenses[m1]
            shapes_lens = shapes[m2]

            # compute effective Sigma_crit
            z_phot, cz = getSigmaCritCorrection(specz_calib, lens[lens_z_key])
            sigma_crit = getSigmaCritEffective(z_phot, cz, shapes_lens['ZP_2'])
            # determine extent in DeltaSigma array
            lower, upper = last_element['all'], last_element['all'] + len(m2)
            elements = np.arange(lower, upper, dtype='int64')
            DeltaSigma[lower:upper] = sigma_crit * tangentialShear(shapes_lens['RA'], shapes_lens['DEC'], shapes_lens['S1'], -shapes_lens['S2'], lens['RA'], lens['DEC'], computeB=False)
            radius[lower:upper] = d12
            weight[lower:upper] = 0.2/(0.2**2 + (0.1*20/shapes_lens['SNR'])**2)**0.5/sigma_crit**2
            last_element['all'] += len(m2)

            # get indices for all sources in each slice
            for i in xrange(len(keys)):
                if callable(keys[i]):
                    key_name = keys[i].__name__
                else:
                    key_name = keys[i]
                values = getValues(shapes_lens, keys[i])
                for s in xrange(len(splittings[i])-1):
                    mask = getSliceMask(values, splittings[i][s], splittings[i][s+1])
                    sum_mask = sum(mask)
                    # check if slice array can keep all indices
                    # if not, double its size
                    if last_element[key_name][s] + sum_mask > slices[key_name][s].size:
                        slices[key_name][s].resize(max(slices[key_name][s].size * 2,last_element[key_name][s] + sum_mask))
                        print " (extending slice index list for %s[%d])" % (key_name, s)
                    slices[key_name][s][last_element[key_name][s] : (last_element[key_name][s] + sum_mask)] = elements[mask]
                    last_element[key_name][s] += sum_mask
                    del mask
                del values
            del shapes_lens, elements, z_phot, cz, sigma_crit
            if counter % 100 == 0:
                print counter
            counter += 1
        print counter

        # reduce the size of the containers to last_element
        DeltaSigma.resize((last_element['all']), refcheck=False)
        weight.resize((last_element['all']), refcheck=False)
        radius.resize((last_element['all']), refcheck=False)
        radius *= 60 # distances now in arcmin
        for k,v in slices.iteritems():
            for vv in v.keys():
                slices[k][vv].resize((last_element[k][vv]), refcheck=False)
        # save the entire shebang
        kwargs = {"DeltaSigma": DeltaSigma, "weight": weight, "radius": radius}
        for k,v in slices.iteritems():
            for vv in v.keys():
                argname = k + "_%d" % vv
                kwargs[argname] = slices[k][vv]
        kwargs['splittings'] = splittings
        np.savez(stackfile, **kwargs)
    else:
        # load from previously saved file
        print " loading from file", stackfile
        X = np.load(stackfile)
        DeltaSigma = X['DeltaSigma']
        weight = X['weight']
        radius = X['radius']
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
    print "generating plots..."
    bins = np.arange(1, maxrange, 2)
    # FIXME: needs to be automatically adjusted???
    ylim = (-0.25, 1)
    xlim = (1, maxrange + 1)
    maxcol = 4 # how many columns per row
    rows = (len(keys)-1)/maxcol + 1
    fig = plt.figure(figsize=(12, 4*rows))

    for i in range(len(keys)):
        if callable(keys[i]):
            key_name = keys[i].__name__
        else:
            key_name = keys[i]
        ax = fig.add_subplot(rows, min(maxcol, len(keys)), i+1)
        makeSlicedProfilePlot(ax, radius, DeltaSigma, weight, slices[key_name], splittings[i], key_name)
        ax.plot(xlim, [0,0], 'k:')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.legend(loc='upper right', numpoints=1, frameon=False)
        if i/maxcol == rows-1:
            ax.set_xlabel('Radius [arcmin]')
        if i%maxcol == 0:
            ax.set_ylabel(r'$\Delta\Sigma\ [10^{14} M_\odot \mathrm{Mpc}^{-2}]$')
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

