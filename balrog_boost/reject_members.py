#!/bin/env python

import fitsio, numpy as np
import esutil as eu
from sys import argv

# load cluster member and DES shape catalog
um = fitsio.FITS('/n/des/pmelchior/des/SV/catalogs/redmapper/sva1_gold_1.0.2_run_redmapper_v6.3.3_ubermem_lgt5_catalog_members.fit')
infocat = '/n/des/pmelchior/des/SV/catalogs/v18/des_sv_wl_info.fits'
sc = fitsio.FITS(infocat)

# match both
h = eu.htm.HTM(8)
maxrange = 1./3600 # 1 arcsec
print "matching (%d vs %d) ..." % (sc[1].get_nrows(), um[1].get_nrows())
m1, m2, d12 = h.match(sc[1]['ra'][:], sc[1]['dec'][:], um[1]['RA'][:], um[1]['DEC'][:], maxrange, maxmatch=1)
del d12
print "found %d matches" % m1.size

# reject members based on their membership probability
print "rejecting ..."
probs = um[1]['P'][:][m2]
rejected = np.random.random(size=m2.size) < probs

# matched and rejected members
r1 = m1[rejected]
print "rejected %d (avg = %.3f) " % (r1.size, probs.mean())

# non-rejected galaxies
print "finding non-rejected indices ..."
n1 = np.arange(sc[1].get_nrows())
n1 = n1[np.in1d(n1, r1, assume_unique=True, invert=True)]

# save the two shape catalogs
print "writing output files ..."
shapecat = '/n/des/pmelchior/des/SV/catalogs/v18/des_sv_wl_ngmix.fits'
wl = fitsio.FITS(shapecat)
wlr = fitsio.FITS(shapecat.replace('.fits', '_ubermem.fits'), 'rw')
wlr.write(wl[1][r1])
wlr.close()
wlr = fitsio.FITS(infocat.replace('.fits', '_ubermem.fits'), 'rw')
wlr.write(sc[1][r1])
wlr.close()

wln = fitsio.FITS(shapecat.replace('.fits', '_wo_ubermem.fits'), 'rw')
wln.write(wl[1][n1])
wln.close()
wln = fitsio.FITS(infocat.replace('.fits', '_wo_ubermem.fits'), 'rw')
wln.write(sc[1][n1])
wln.close()




