#!/bin/env python

import fitsio
from matplotlib.colors import LogNorm
from matplotlib.ticker import NullFormatter
import pylab as plt
import weighting
from esutil import numpy_util

def modestify(data):
    galcut = (data['flags'] <= 3) & -( ((data['class_star'] > 0.3) & (data['mag_auto'] < 18.0)) | ((data['spread_model'] + 3*data['spreaderr_model']) < 0.003) | ((data['mag_psf'] > 30.0) & (data['mag_auto'] < 21.0)))
    return galcut

def getUniqueIndex(data):
    unq, unq_idx, unq_cnt = np.unique(data, return_index=True, return_counts=True)
    return unq_idx

def fields_view(arr, fields):
    dtype2 = np.dtype({name:arr.dtype.fields[name] for name in fields})
    return np.ndarray(arr.shape, dtype2, arr, 0, arr.strides)

def getResampledIndex(data, weights, size=1):
    idx = np.empty(size, dtype='i8')
    done = 0
    while done < size:
        draw = np.nonzero(np.random.random(data.shape[0]) < weights/weights.max())[0]
        draw_len = min(size-done, draw.size)
        idx[done:done+draw_len] = draw[:draw_len]
        done += draw_len
    return idx

def plotHistograms(balrog_data, weights, des_data, band):
    # plot results
    bins = [np.linspace(20, 25, 51), np.linspace(0.5, 6, 45)]

    fig = plt.figure(figsize=(12, 5))
    ax = fig.add_subplot(131)
    ax.hist2d(balrog_data['mag_auto'], balrog_data['flux_radius'], bins=bins, normed=LogNorm, cmap=plt.cm.GnBu_r)
    ax.set_ylabel('FLUX_RADIUS ' + band)
    ax.set_xlabel('MAG_AUTO ' + band)
    ax.text(0.05, 0.95, 'Balrog', color='w', ha='left', va='top', transform=ax.transAxes)

    ax = fig.add_subplot(132)
    ax.hist2d(balrog_data['mag_auto'], balrog_data['flux_radius'], bins=bins, normed=LogNorm, cmap=plt.cm.GnBu_r, weights=weights)
    ax.set_xlabel('MAG_AUTO ' + band)
    ax.yaxis.set_major_formatter(NullFormatter())
    ax.text(0.05, 0.95, 'Balrog re-weighted', color='w', ha='left', va='top', transform=ax.transAxes)

    ax = fig.add_subplot(133)
    ax.hist2d(des_data[:,0], des_data[:,1], bins=bins, normed=LogNorm, cmap=plt.cm.GnBu_r)
    ax.set_xlabel('MAG_AUTO ' + band)
    ax.yaxis.set_major_formatter(NullFormatter())
    ax.text(0.05, 0.95, 'DES shapes', color='w', ha='left', va='top', transform=ax.transAxes)
    plt.subplots_adjust(left=0.06, right=0.97, top=0.95, wspace=0)
    plt.show()


band = 'i'

# get balrog data and eliminate duplicates
balrog = fitsio.FITS('/n/des/pmelchior/des/SV/catalogs/Balrog/BalrogObsFile-' + band + '.fits')
w = balrog[1].where('flags == 0')
balrog_ids = balrog[1][w]['balrog_index'].astype('i8')
balrog_unq_idx = getUniqueIndex(balrog_ids)
w = w[balrog_unq_idx]

# load data for matching and to make modest cut
balrog_data = balrog[1]['ra', 'dec', 'balrog_index', 'mag_auto', 'flux_radius', 'flags', 'class_star', 'spread_model', 'spreaderr_model', 'mag_psf'][w]
modest = modestify(balrog_data)
balrog_data = balrog_data[modest]
del balrog_ids, balrog_unq_idx, modest, w

# get that DES shape catalog
des = fitsio.FITS('/n/des/pmelchior/des/SV/catalogs/v18/des_sv_wl_info.fits')
wd = des[1].where('sva1_flags == 0 && sva1_gold_flags == 0 && sva1_spte_flags == 0 && sva1_gold_mag_flags == 0 && photoz_bin >= 0 && ngmix_flags == 0 && dec < -40 && 60 < ra && ra< 95')
des_mag = des[1][wd]['mag_auto_' + band]
des_ids = des[1][wd]['coadd_objects_id']

# get DES flux_radius from secondary table
des2 = fitsio.FITS('/n/des/pmelchior/des/SV/catalogs/sva1_coadd_object_flux_radius.fits')
des2_ids = des2[1]['coadd_objects_id'][:]
des2_index, idx_ = numpy_util.match(des2_ids, des_ids)
des_size = des2[1]['flux_radius_' + band][:][des2_index]
des_data = np.dstack((des_mag, des_size))[0]
del des_mag, des_size, des_ids, des2_ids, idx_

# match 2D distributions according to nearest neighbors
# separate the 3 photoz-bins
weights = []
n_near = 20
balrog_data_ = np.dstack((balrog_data['mag_auto'], balrog_data['flux_radius']))[0]
newcat = np.empty(des_data.shape[0], dtype=[('ra', 'f4'), ('dec', 'f4'), ('photoz_bin', 'i1')])
done = 0
for pzb in [0,1,2]:
    mask = des[1][wd]['photoz_bin'][:] == pzb
    wn = weighting.weight_match(balrog_data_, des_data[mask], n_near)
    weights.append(wn.get_weights())
    # create a matched resmpled balrog catalog with ra,dec,photoz_bin
    newidx = getResampledIndex(balrog_data, weights[pzb], size=mask.sum())
    newcat[done:done+newidx.size]['ra'] = balrog_data[newidx]['ra']
    newcat[done:done+newidx.size]['dec'] = balrog_data[newidx]['dec']
    newcat[done:done+newidx.size]['photoz_bin'] = pzb
    done += newidx.size
newfits = fitsio.FITS('/n/des/pmelchior/des/SV/catalogs/v18/balrog_matched_ngmix_i_photoz-bin.fits', 'rw')
newfits.write(newcat)
newfits.close()


"""
# cross-check: i-band weights applied to other bands
for band in ['g', 'r', 'z']:
    des_mag = des[1][wd]['mag_auto_' + band]
    des_size = des2[1]['flux_radius_' + band][:][des2_index]

    balrog2 = fitsio.FITS('/n/des/pmelchior/des/SV/catalogs/Balrog/BalrogObsFile-' + band + '.fits')
    w2 = balrog2[1].where('flags == 0')
    balrog2_ids = balrog2[1][w2]['balrog_index'].astype('i8')
    balrog2_unq_idx = getUniqueIndex(balrog2_ids)
    w2 = w2[balrog2_unq_idx]

    balrog2_data = balrog2[1]['balrog_index', 'mag_auto', 'flux_radius', 'flags', 'class_star', 'spread_model', 'spreaderr_model', 'mag_psf'][w2]
    modest = modestify(balrog2_data)
    balrog2_data = balrog2_data[modest]
    del balrog2_ids, balrog2_unq_idx, modest, w2

    idx, idx2 = numpy_util.match(balrog_data['balrog_index'], balrog2_data['balrog_index'])
    
    fig = plt.figure(figsize=(12, 5))
    ax = fig.add_subplot(131)
    ax.hist2d(balrog2_data['mag_auto'][idx2], balrog2_data['flux_radius'][idx2], bins=bins, normed=LogNorm, cmap=plt.cm.GnBu_r)
    ax.set_ylabel('FLUX_RADIUS ' + band)
    ax.set_xlabel('MAG_AUTO ' + band)
    ax.text(0.05, 0.95, 'Balrog', color='w', ha='left', va='top', transform=ax.transAxes)
                                 
    ax = fig.add_subplot(132)
    ax.hist2d(balrog2_data['mag_auto'][idx2], balrog2_data['flux_radius'][idx2], bins=bins, normed=LogNorm, cmap=plt.cm.GnBu_r, weights=weights[idx])
    ax.set_xlabel('MAG_AUTO ' + band)
    ax.yaxis.set_major_formatter(NullFormatter())
    ax.text(0.05, 0.95, 'Balrog re-weighted', color='w', ha='left', va='top', transform=ax.transAxes)
                                 
    ax = fig.add_subplot(133)
    ax.hist2d(des_mag, des_size, bins=bins, normed=LogNorm, cmap=plt.cm.GnBu_r)
    ax.set_xlabel('MAG_AUTO ' + band)
    ax.yaxis.set_major_formatter(NullFormatter())
    ax.text(0.05, 0.95, 'DES shapes', color='w', ha='left', va='top', transform=ax.transAxes)
    plt.subplots_adjust(left=0.06, right=0.97, top=0.95, wspace=0)
    plt.show()
"""
