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
    return bc

""" for plotting only"""
def lon2RA(lon):
    lon = 360 - lon
    hours = int(lon)/15
    minutes = int(float(lon - hours*15)/15 * 60)
    minutes = '{:>02}'.format(minutes)
    return "%d:%sh" % (hours, minutes)

def getCountLocation(config, shapes, nside=512):
    ipix = hp.ang2pix(nside, (90-shapes[config['shape_dec_key']])/180*np.pi, shapes[config['shape_ra_key']]/180*np.pi, nest=False)
    bc = np.bincount(ipix)
    pixels = np.nonzero(bc)[0]
    bc = bc[bc>0] / hp.nside2resol(nside, arcmin=True)**2 # in arcmin^-2
    theta, phi = hp.pix2ang(nside, pixels, nest=False)
    lat = 90 - theta*180/np.pi
    lon = phi*180/np.pi
    return bc, lat, lon

from mpl_toolkits.basemap import Basemap
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def plotDensityMap(config, shapes, nside=512):
    # set up figure
    setTeXPlot(2*nside/512)
    fig = plt.figure(figsize=(6.5*nside/512,6*nside/512))
    ax = fig.add_axes([0.07,0.07,0.84,0.9], aspect='equal')
    # equal-area map straight above the footprint center
    m = Basemap(projection='aea',width=2000000,height=2200000,
                lat_0=-52.5, lat_1=-61, lat_2=-42., lon_0=-75.)

    # after cuts
    vmin,vmax = 0,10
    bc, lat, lon = getCountLocation(config, shapes, nside=nside)
    x,y  = m(-lon, lat)
    sc = m.scatter(x,y,c=bc, linewidths=0, s=10, marker='s', cmap=cm.YlOrRd, vmin=vmin, vmax=vmax, rasterized=True, ax=ax)
    #sc = m.scatter(x,y,c=bc, linewidths=0, s=8, marker='h', cmap=cm.jet, vmin=vmin, vmax=vmax, rasterized=True)#, norm=matplotlib.colors.LogNorm())

    # draw parallels and meridians.
    # label on left and bottom of map.
    parallels = np.arange(-75.,0.,5.)
    m.drawparallels(parallels,labels=[1,0,0,0], labelstyle="+/-", linewidth=0.5)
    meridians = np.arange(0.,360.,5.)
    m.drawmeridians(meridians,labels=[0,0,0,1], fmt=lon2RA, linewidth=0.5)

    # add colorbar
    cb = m.colorbar(sc,"right", size="3%", pad='0%')
    cb.set_label('$n_g\ [\mathrm{arcmin}^{-2}]$')
    cb.solids.set_edgecolor("face")
    #plt.show()
    plt.savefig('depth_map_quadrant_check.pdf', transparent=True)
    plt.savefig('depth_map_quadrant_check.png')
""" end plotting """
    
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

    if config['coords'] not in ['angular', 'physical']:
        print "config: specify either 'angular' or 'physical' coordinates"
        raise SystemExit
    
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
    outdir = os.path.dirname(configfile) + "/"
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
    zbin_edges = [0.3, 0.55, 0.83, 1.3]
    
    if shapes.size:
        basename = os.path.basename(shapefile)
        basename = ".".join(basename.split(".")[:-1])
        densityfile = outdir + basename + '_density.fits'

        # make healpix map of density of all shapes
        makeDensityMap(densityfile, config, shapes, nside=1024)
        print "created healpix density map %s" % densityfile
        dmap=hu.readDensityMap(densityfile)
        #plotDensityMap(config, shapes, nside=1024)

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
        """
        
        # open lens catalog for quadrant check
        # we need to remove any lens cuts since we want the check for all
        # lenses in the lens_file
        config['lens_cuts'] = []
        lenses = getLensCatalog(config, verbose=True)

        # check quadrants around the input points
        # make sure weighted position ellipticity in adjacent quadrants
        # less than 0.05
        ellip_max=0.05
        data = np.zeros(lenses.size, dtype=[('quad_flags', 'i1')])

        # match the outer radius to the range asked for stacking
        if config['coords'] == "physical":
            radius_degrees = Dist2Ang(config['maxrange'], lenses[config['lens_z_key']])
        else:
            radius_degrees = config['maxrange'] * np.ones(lenses.size)
            
        for i in xrange(lenses.size):
            lens = lenses[i]
            """
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
            """
            data['quad_flags'][i] = dmap.check_quad(lens[config['lens_ra_key']], lens[config['lens_dec_key']], radius_degrees[i], ellip_max)
        
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
