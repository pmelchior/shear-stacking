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

    # see if we need to do anything
    append_to_extra = False
    try:
        hdu = fitsio.FITS(config['lens_extra_file'])
        columns = hdu[1].get_colnames()
        hdu.close()
        if 'quad_flags' in columns:
            print "Quadrant check flags already in " + config['lens_extra_file']
            print "Delete file if you want to regenerate them."
            raise SystemExit
        else:
            append_to_extra = True
    except (KeyError, IOError) as exc: # not in config or file doesn't exist
        pass
    
    
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
        data = np.zeros(lenses.size, dtype=[('quad_flags', 'i1')])
        # 30 Mpc/h outer radius
        radius_degrees = Dist2Ang(30., lenses[config['lens_z_key']])
        for i in xrange(lenses.size):
            lens = lenses[i]
            zl = lens[config['lens_z_key']]
            # use the dmap the is has shapes from bin above the cluster redshift
            # FIXME: we lose all clusters above 0.85, but that's only
            # a small fraction
            dmap_ = None
            for b in xrange(3):
                zmin = zbin_edges[b]
                zmax = zbin_edges[b+1]
                if zl < zmin:
                    dmap_ = dmaps[b]
                    break
            # get quadrant ellipticty and mask flags
            if dmap_ is not None:
                data['quad_flags'][i] = dmap_.check_quad(lens[config['lens_ra_key']], lens[config['lens_dec_key']], radius_degrees[i], ellip_max)
        
        # save result as table
        if append_to_extra == False:
            lensfile = config['lens_file']
            basename = os.path.basename(lensfile)
            basename = ".".join(basename.split(".")[:-1])
            quadfile = outdir + basename + '_quadrant-check.fits'
            fits = fitsio.FITS(quadfile, 'rw', clobber=True)
            fits.write(data)
            fits.close()
            print "created quadrant check file %s" % quadfile

            if 'lens_extra_file' not in config.keys():
                print "\nBefore proceeding: add"
                print "    \"lens_extra_file\": \"%s\"" % quadfile
                print "to your config file!"
        else:
            fits = fitsio.FITS(config['lens_extra_file'], 'rw')
            fits[1].insert_column('quad_flags', data['quad_flags'])
            fits.close()    
