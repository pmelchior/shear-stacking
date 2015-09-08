#!/bin/env python

import fitsio, numpy as np, weighting
from esutil import numpy_util
import desdb

def getUniqueIndex(data):
    unq, unq_idx, unq_cnt = np.unique(data, return_index=True, return_counts=True)
    return unq_idx

def fields_view(arr, fields):
    dtype2 = np.dtype({name:arr.dtype.fields[name] for name in fields})
    return np.ndarray(arr.shape, dtype2, arr, 0, arr.strides)

def getResampledIndex(data, weights, size=1):
    """idx = np.empty(size, dtype='i8')
    done = 0
    while done < size:
        draw = np.nonzero(np.random.random(data.shape[0]) < weights/weights.max())[0]
        draw_len = min(size-done, draw.size)
        idx[done:done+draw_len] = draw[:draw_len]
        done += draw_len
    """
    idx = np.random.choice(np.arange(data.shape[0], dtype='i8'), size=size, replace=False, p=weights)
    return idx

def saveWeightFile(filename, data, weights):
    newdata = np.empty(data.shape[0], dtype=[('ra', 'f4'), ('dec', 'f4'), ('mag_auto', 'f4'), ('flux_radius', 'f4'), ('match_weight', 'f8')])
    newdata['ra'] = data['ra']
    newdata['dec'] = data['dec']
    newdata['mag_auto'] = data['mag_auto']
    newdata['flux_radius'] = data['mag_auto']
    newdata['match_weight'] = weights
    fits = fitsio.FITS(filename, 'rw')
    fits.write(newdata)
    fits.close()

band = 'i'

balrog = fitsio.FITS('/n/des/suchyta.1/DataFiles/BalrogDB/cluster-lensing-fix-cut/sim-clus-lens-like.fits')
w = balrog[1].where('modest_i == 0')
balrog_ids = balrog[1][w]['balrog_index'].astype('i8')
balrog_version = balrog[1][w]['version']
# create almost unique superindex and select unique ones
balrog_unq_idx = getUniqueIndex(balrog_ids + 100000000 * balrog_version)
w = w[balrog_unq_idx]
balrog_data = balrog[1]['alphawin_j2000_' + band, 'deltawin_j2000_' + band, 'mag_auto_' + band, 'flux_radius_' + band][w]

# get that DES shape catalog
des = fitsio.FITS('/n/des/pmelchior/des/SV/catalogs/v18/des_sv_wl_info.fits')
wd = des[1].where('sva1_flags == 0 && sva1_gold_flags == 0 && sva1_spte_flags == 0 && sva1_gold_mag_flags == 0 && photoz_bin >= 0 && ngmix_flags == 0 && dec < -40 && 60 < ra && ra < 95')
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
n_near = 10
balrog_data_ = np.dstack((balrog_data['mag_auto_' + band], balrog_data['flux_radius_' + band]))[0]
newcat = np.empty(des_data.shape[0], dtype=[('ra', 'f4'), ('dec', 'f4'), ('photoz_bin', 'i1')])

done = 0
for pzb in [0,1,2]:
    mask = (des_pzb == pzb)
    wn = weighting.weight_match(balrog_data_, des_data[mask], n_near)
    weights = wn.get_weights()
    
    # save a balrog catalog with weights
    #filename = '/n/des/pmelchior/des/SV/catalogs/v18/balrog_matched_ngmix_i_weights_%d.fits' % pzb
    #saveWeightFile(filename, balrog_data, weights)

    # create a matched resmpled balrog catalog with ra,dec,photoz_bin
    newidx = getResampledIndex(balrog_data, weights, size=mask.sum())
    newcat[done:done+newidx.size]['ra'] = balrog_data[newidx]['alphawin_j2000_' + band]
    newcat[done:done+newidx.size]['dec'] = balrog_data[newidx]['deltawin_j2000_' + band]
    newcat[done:done+newidx.size]['photoz_bin'] = pzb
    done += newidx.size
    
newfits = fitsio.FITS('/n/des/pmelchior/des/SV/catalogs/v18/balrog_matched_ngmix_i_photoz-bin.fits', 'rw')
newfits.write(newcat)
newfits.close()
