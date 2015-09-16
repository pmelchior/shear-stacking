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

def skyDistance(ra, dec, ra_ref, dec_ref):
    # CAUTION: this needs to be a pseudo-Cartesian coordinate frame
    # (not pure RA/DEC), otherwise distances are skewed
    return (((ra-ra_ref)*cos(dec*pi/180))**2 + (dec-dec_ref)**2)**0.5

band = 'i'

balrog = fitsio.FITS('/n/des/suchyta.1/DataFiles/BalrogDB/cluster-lensing-3/sim-clus-lens-like.fits')
w = balrog[1].where('modest_i == 0')
balrog_ids = balrog[1][w]['balrog_index'].astype('i8')
balrog_version = balrog[1][w]['version']
# create almost unique superindex and select unique ones
balrog_unq_idx = getUniqueIndex(balrog_ids + 100000000 * balrog_version)
w = w[balrog_unq_idx]
balrog_info = balrog[1]['alphawin_j2000_' + band, 'deltawin_j2000_' + band, 'tilename_' + band, 'mag_auto_' + band, 'flux_radius_' + band][w]
balrog_data = np.dstack((balrog[1]['mag_auto_' + band][w], balrog[1]['flux_radius_' + band][w], balrog[1]['alphawin_j2000_' + band][w], balrog[1]['deltawin_j2000_' + band][w]))[0]
del balrog_ids, balrog_version, balrog_unq_idx, w

# get that DES shape catalog
des = fitsio.FITS('/n/des/pmelchior/des/SV/catalogs/v18/des_sv_wl_info.fits')
w = des[1].where('sva1_flags == 0 && sva1_gold_flags == 0 && sva1_spte_flags == 0 && sva1_gold_mag_flags == 0 && photoz_bin >= 0 && ngmix_flags == 0 && dec < -40 && 60 < ra && ra < 95')
des_mag = des[1][w]['mag_auto_' + band]
des_ids = des[1][w]['coadd_objects_id']
des_pzb = des[1][w]['photoz_bin']

# get DES flux_radius and tilename from secondary table
des2 = fitsio.FITS('/n/des/pmelchior/des/SV/catalogs/sva1_coadd_object_tile_flux_radius.fits')
des2_ids = des2[1]['coadd_objects_id'][:]
des2_index, idx_ = numpy_util.match(des2_ids, des_ids)
des_size = des2[1]['flux_radius_' + band][:][des2_index]
des_tile = des2[1]['tilename'][:][des2_index]
des_data = np.dstack((des_mag, des_size, des[1][w]['ra'], des[1][w]['dec']))[0]
tiles = np.unique(des_tile)

# normalize mag/size of data vector to have ~ mean zero and unit variance
# use global instead of tile means/stds
des_means = np.median(des_data, axis=0)
des_stds = np.std(des_data, axis=0)
des_data[:,0] = (des_data[:,0]-des_means[0])/des_stds[0]
des_data[:,1] = (des_data[:,1]-des_means[1])/des_stds[1]
balrog_data[:,0] = (balrog_data[:,0]-des_means[0])/des_stds[0]
balrog_data[:,1] = (balrog_data[:,1]-des_means[1])/des_stds[1]

# check if matching is complete and unique
if not (des_ids == des2_ids[des2_index]).all():
    raise RuntimeError("matching all DES objects from secondary table failed!")
del des_ids, w, des2_ids, des2_index, idx_, des_mag, des_size 

# match both data sets distributions according to nearest neighbors
# separate the 3 photoz-bins
n_near = 10
oversampling = 1
newcat = np.empty(int(des_data.shape[0]*oversampling), dtype=[('ra', 'f4'), ('dec', 'f4'), ('mag_auto_' + band, 'f4'), ('flux_radius_' + band, 'f4'), ('photoz_bin', 'i1')])
done = 0

for tile in tiles:
    des_tile_mask = (des_tile == tile)
    balrog_tile_mask = (balrog_info['tilename_' + band] == tile)
    print "Tile "+ tile + ", matching ", des_tile_mask.sum() , "vs", balrog_tile_mask.sum()

    if des_tile_mask.sum() > balrog_tile_mask.sum():
        print "Balrog data missing"
        raise SystemExit
    
    if balrog_tile_mask.sum() > n_near:
        balrog_info_ = balrog_info[balrog_tile_mask]
        balrog_data_ = balrog_data[balrog_tile_mask]
        des_data_ = des_data[des_tile_mask]
        des_pzb_ = des_pzb[des_tile_mask]

        # normalize the ra/dec coordinates of data vector
        # project onto tangent plane
        # use balrog sample because of higher spatial density
        ra_mean, dec_mean = np.median(balrog_data_[:,2]), np.median(balrog_data_[:,3])
        balrog_data_[:,2] = (balrog_data_[:,2] - ra_mean)*np.cos(dec_mean*np.pi/180)
        balrog_data_[:,3] = (balrog_data_[:,3] - dec_mean)
        des_data_[:,2] = (des_data_[:,2] - ra_mean)*np.cos(dec_mean*np.pi/180)
        des_data_[:,3] = (des_data_[:,3] - dec_mean)
        ra_std, dec_std = balrog_data_[:,2].std(), balrog_data_[:,3].std()
        # give higher weight to mag/size
        ra_std *= 10.
        dec_std *= 10.
        balrog_data_[:,2] /= ra_std
        balrog_data_[:,3] /= dec_std
        des_data_[:,2] /= ra_std
        des_data_[:,3] /= dec_std

        # cycle through photo-z bins at random
        # because we remove selected Balrog gaalxies from sample,
        # the last bin would otherwise suffer most from shrinking sampling
        for pzb in np.random.permutation(3):
            des_pzb_mask = (des_pzb_ == pzb)
            if des_pzb_mask.sum() > n_near:
                wn = weighting.weight_match(balrog_data_, des_data_[des_pzb_mask], n_near)
                weights = wn.get_weights()
                weights /= weights.sum()

                # create a matched resampled balrog catalog with ra,dec,photoz_bin
                newidx = getResampledIndex(balrog_data_, weights, size=int(des_pzb_mask.sum()*oversampling))
                newcat[done:done+newidx.size]['ra'] = balrog_info_[newidx]['alphawin_j2000_' + band]
                newcat[done:done+newidx.size]['dec'] = balrog_info_[newidx]['deltawin_j2000_' + band]
                newcat[done:done+newidx.size]['mag_auto_' + band] = balrog_info_[newidx]['mag_auto_' + band]
                newcat[done:done+newidx.size]['flux_radius_' + band] = balrog_info_[newidx]['flux_radius_' + band]
                newcat[done:done+newidx.size]['photoz_bin'] = pzb
                done += newidx.size

                # remove newidx from balrog_data_ to prevent multiple selection
                balrog_data_ = np.delete(balrog_data_, newidx, axis=0)
                balrog_info_ = np.delete(balrog_info_, newidx, axis=0)

newfits = fitsio.FITS('/n/des/pmelchior/des/SV/catalogs/v18/balrog_matched_ngmix_tiles_radec10_near10_i_photoz-bin.fits', 'rw')
newfits.write(newcat[:done])
newfits.close()






