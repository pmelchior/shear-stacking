#!/bin/env python

import json, errno
import healpy as hp
import healpix_util as hu
import numpy as np
from sys import argv
from shear_stacking import *

def makeDensityMap(outfile, config, shapes, nside=512):
    ipix = hp.ang2pix(nside, (90-shapes[config['shape_dec_key']])/180*np.pi, shapes[config['shape_ra_key']]/180*np.pi, nest=False)
    bc = np.bincount(ipix, minlength=hp.nside2npix(nside))
    hp.write_map(outfile, bc)

if __name__ == '__main__':
    # parse inputs
    try:
        configfile = argv[1]
    except IndexError:
        print "usage: " + argv[0] + " <config file>"
        raise SystemExit
    try:
        fp = open(configfile)
        print "opening configfile " + configfile
        config = json.load(fp)
        fp.close()
    except IOError:
        print "configfile " + configfile + " does not exist!"
        raise SystemExit

    outdir = os.path.dirname(configfile)
    if outdir[-1] != '/':
        outdir += '/'
    try:
        os.makedirs(outdir)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(outdir):
            pass
        else: raise
    
    # open shape catalog
    shapefile = config['shape_file']
    # since all selection are normally in the extra file (if present)
    # we speed up the process by changing extra and shape and dropping shape
    try:
        extrafile = config['shape_file_extra']
        config['shape_file'] = extrafile
        del config['shape_file_extra']
    except KeyError:
        pass
    shapes = getShapeCatalog(config, verbose=True)

    # Troxel's photo-z bin edges
    zbin_edges = [0.3, 0.57, 0.85, 1.3]
    
    if shapes.size:
        basename = os.path.basename(shapefile)
        basename = ".".join(basename.split(".")[:-1])
        densityfile = outdir + basename + '_density.fits'

        """
        # make healpix map of density of all shapes
        makeDensityMap(densityfile, config, shapes, nside=1024)
        print "created healpix density map %s" % densityfile
        dmap=hu.readDensityMap(densityfile)
        """
        
        # define shapes selection that are at higher than cluster
        # redshift, i.e. no high-z cutoff
        dmaps = []
        for b in xrange(3):
            mask = shapes[config['shape_z_key']] >= b
            print "bin %d: %d shapes" % (b, mask.sum())
            densityfile = outdir + basename + '_density_bin%d.fits' %b
            makeDensityMap(densityfile, config, shapes[mask], nside=1024)
            print "created healpix density map %s" % densityfile
            dmaps.append(hu.readDensityMap(densityfile))
        
        # open lens catalog for quadrant check
        lenses = getLensCatalog(config, verbose=True)

        # check quadrants around the input points
        # make sure weighted position ellipticity in adjacent quadrants
        # less than 0.05
        ellip_max=0.05
        data = np.zeros(lenses.size, dtype=[('quad_flags', 'i1'), ('quad_ellip', '4f4')])
        # 30 Mpc/h outer radius
        radius_degrees = Dist2Ang(30., lenses[config['lens_z_key']])
        for i in xrange(lenses.size):
            lens = lenses[i]
            zl = lens[config['lens_z_key']]
            # use the dmap the is has shapes from bin above the cluster redshift
            # FIXME:we lose all clusters above 0.85, but that's a small fraction
            dmap_ = None
            for b in xrange(3):
                zmin = zbin_edges[b]
                zmax = zbin_edges[b+1]
                if zl < zmin:
                    dmap_ = dmaps[b]
                    break
            # get quadrant ellipticty and mask flags
            if dmap_ is not None:
                data['quad_flags'][i] = dmap_.check_quad(lens['RA'], lens['DEC'], radius_degrees[i], ellip_max)
                data['quad_ellip'][i] = dmap_.get_quad_ellip(lens['RA'], lens['DEC'], radius_degrees[i])
        
        # save result as table
        lensfile = config['lens_file']
        basename = os.path.basename(lensfile)
        basename = ".".join(basename.split(".")[:-1])
        quadfile = outdir + basename + '_quadrant-check.fits'
        fits = fitsio.FITS(quadfile, 'rw', clobber=True)
        fits.write(data)
        fits.close()
        print "created quadrant check file %s" % quadfile
        
