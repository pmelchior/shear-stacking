import numpy as np
from math import pi, sqrt
import os, fitsio

def skyAngle(ra, dec, ra_ref, dec_ref):
    # CAUTION: this needs to be a pseudo-Cartesian coordinate frame
    # (not pure RA/DEC), otherwise angles are skewed
    return np.arctan2(dec-dec_ref, (ra-ra_ref)*np.cos(dec*pi/180))

def skyDistance(ra, dec, ra_ref, dec_ref):
    # CAUTION: this needs to be a pseudo-Cartesian coordinate frame
    # (not pure RA/DEC), otherwise distances are skewed
    return (((ra-ra_ref)*np.cos(dec*pi/180))**2 + (dec-dec_ref)**2)**0.5

def tangentialShear(ra, dec, e1, e2, ra_ref, dec_ref, computeB=False):
    phi = skyAngle(ra, dec, ra_ref, dec_ref)
    if computeB is False:
        return -e1*np.cos(2*phi) + e2*np.sin(2*phi)
    else:
        return -e1*np.cos(2*phi) + e2*np.sin(2*phi), e1*np.sin(2*phi) + e2*np.cos(2*phi)

# CAUTION: assumes Gaussian errors and large samples
# replace with Jackknife/Bootstrap estimate for more accurate errors
class WeightedMeanVar:
    def __init__(self):
        self.N = 0.
        self.Wi = 0.
        self.WiXi = 0.
        self.WiXi2 = 0.
        self.WiSi = 0.
    def getMean(self):
        if self.Wi > 0:
            if self.WiSi > 0:
                return self.WiXi / self.WiSi
            else:
                return self.WiXi / self.Wi
        else:
            return 0
    def getStd(self):
        if self.Wi > 0:
            if self.WiSi > 0:
                # this is not entirely correct since we ignore the extra variance 
                # in the sensitivity itself
                # again: use bootstraps of the mean for more accurate errors
                return ((self.WiXi2 - (self.WiXi**2)/self.Wi) / ((self.N - 1) * self.WiSi))**0.5
            else:
                return ((self.WiXi2 - (self.WiXi**2)/self.Wi) / ((self.N - 1) * self.Wi))**0.5
        else:
            return 0
    def insert(self, X, W, S=None):
        if X.size:
            self.N += X.size
            self.Wi += W.sum()
            self.WiXi += (W*X).sum()
            self.WiXi2 += (W*X**2).sum()
            if S is not None:
                self.WiSi += (W*S).sum()
    def __iadd__(self, other):
        self.N += other.N
        self.Wi += other.Wi
        self.WiXi += other.WiXi
        self.WiXi2 += other.WiXi2
        self.WiSi += other.WiSi
        return self

class BinnedScalarProfile:
    def __init__(self, bins):
        self.bins = bins
        self.Q = [] # binned quantity
        self.R = [] # center of radial bins
        for i in xrange(len(self.bins)-1):
            self.Q.append(WeightedMeanVar())
            self.R.append(0.)
    def __iadd__(self, other):
        if len(self.R) == len(other.R):
            for i in xrange(len(self.bins)-1):
                self.Q[i] += other.Q[i]
                self.R[i] += other.R[i]
            return self
        else:
            raise AssertionError("Profiles do not have the same length.")
    def insert(self, R, Q, W, S=None):
        for i in xrange(len(self.bins)-1):
            mask = (R >= self.bins[i]) & (R < self.bins[i+1])
            if S is None:
                self.Q[i].insert(Q[mask], W[mask])
            else:
                self.Q[i].insert(Q[mask], W[mask], S[mask])
            self.R[i] += R[mask].sum()
            del mask
    def getProfile(self):
        mean_q = np.empty(len(self.bins)-1)
        std_q = np.empty(len(self.bins)-1)
        n = np.empty(len(self.bins)-1)
        r = np.empty(len(self.bins)-1)
        for i in xrange(len(self.bins)-1):
            n[i] = self.Q[i].N
            if n[i] > 0:
                r[i] = self.R[i] / n[i]
                mean_q[i] = self.Q[i].getMean()
                std_q[i] = self.Q[i].getStd()
        return r, n, mean_q, std_q
    """
    def getSummmedProfile(self):
        mean_q = np.empty(len(self.bins)-1)
        std_q = np.empty(len(self.bins)-1)
        n = np.empty(len(self.bins)-1)
        r = np.empty(len(self.bins)-1)
        for i in xrange(len(self.bins)-1):
            n[i] = self.Q[i].N
            if n[i] > 0:
                r[i] = self.R[i] / n[i]
                mean_q[i] = self.Q[i].WiXi / (np.pi*(self.bins[i+1]**2 - self.bins[i]**2))
                std_q[i] = ((self.Q[i].Wi*self.Q[i].WiXi2 - self.Q[i].WiXi**2) / (self.Q[i].N - 1) * self.Q[i].N)**0.5
            else:
                r[i] = 0
                mean_q[i] = 0
                std_q[i] = 0
        return r, n, mean_q, std_q
    """

# extrapolation function from
# http://stackoverflow.com/questions/2745329/how-to-make-scipy-interpolate-give-an-extrapolated-result-beyond-the-input-range
def extrap(x, xp, yp):
    """np.interp function with linear extrapolation"""
    x_ = np.array(x)
    y = np.array(np.interp(x_, xp, yp))
    y[x_ < xp[0]] = yp[0] + (x_[x_ < xp[0]] -xp[0]) * (yp[0] - yp[1]) / (xp[0] - xp[1])
    y[x_ > xp[-1]] = yp[-1] + (x_[x_ > xp[-1]] -xp[-1])*(yp[-1] - yp[-2])/(xp[-1] - xp[-2])
    return y  

from galsim import Cosmology
cosmo = Cosmology()

# get separation in deg for distance L in Mpc/h at redshift z
# uses c/H0 = 3000 Mpc/h
def Dist2Ang(L, z):
    global cosmo
    return L / cosmo.Da(z) / 3000. * 180./np.pi

def Ang2Dist(theta, z):
    global cosmo
    return theta * cosmo.Da(z) * 3000. / 180. * np.pi

def getBeta(z_c, z):
    if z_c >= z:
        return 0
    else:
        return cosmo.Da(z, z_c)/cosmo.Da(z)  

def getSigmaCrit(z_c, z):
    c2_4piG = 3.882 # in 1e14 M_solar / Mpc^2 (since cosmo.Da comes in units of c/H0)
    return c2_4piG / getBeta(z_c, z) / cosmo.Da(z_c)

# From Troxel: <Sigma_crit ^-power w> / <w> for each photo-z bin
# calculated for flat LCDM model with Omega_m = 0.27, h=0.7 and distances in Mpc
# FIXME: need to adjust to Omega_m = 0.3 of our reference cosmology
def getWZ(power=1):
    thisdir = os.path.dirname(os.path.realpath(__file__))
    if power != 1 and power != 2:
        raise RuntimeError("Must be integer power 1 or 2")
    filename = 'invsigcrit-skynetsmooth6-false_z_mean.txt'
    if power == 2:
        filename = 'invsigcrit2-skynetsmooth6-false_z_mean.txt'
    data = np.genfromtxt(thisdir + '/data/' + filename, dtype=[('z', 'float32'), ('bin0', 'float32'), ('bin1', 'float32'), ('bin2', 'float32')])

    c2_4piG = 1.661e4 # in 1e14 M_solar / Mpc, for distances in Mpc
    for b in xrange(3):
        data['bin%d' % b] /= c2_4piG**power
    return data


# join to names arrays of with different columns but same number of rows
# if column names are identical, those of data1 will prevail
# FIXME: does not work if data2 has non-trivial dtypes (e.g. multi-dim array)
def joinDataSets(data1, data2):
    from numpy.lib import recfunctions
    columns = data1.dtype.names    
    columns_ = []
    dtypes_ = []
    for i in xrange(len(data2.dtype.names)):
        col = data2.dtype.names[i]
        if col not in columns:
            columns_.append(col)
            dtypes_.append(data2.dtype[col])
    return recfunctions.rec_append_fields(data1, columns_, [data2[c] for c in columns_], dtypes=dtypes_)

def getShapeCatalog(config, verbose=False, chunk_index=None):
    # open shapes file(s)
    shapefile = config['shape_file']
    chunk_size = config['shape_chunk_size']
    shdu = fitsio.FITS(shapefile)
    extra = None
    if verbose:
        print "opening shapefile %s (%d entries)" % (shapefile, shdu[1].get_nrows())

    if len(config['shape_cuts']) == 0:
        if chunk_index is None:
            shapes = shdu[1][:]
        else:
            shapes = shdu[1][chunk_index*chunk_size : (chunk_index+1)*chunk_size]
        try:
            ehdu = fitsio.FITS(config['shape_file_extra'])
            if verbose:
                print "opening extra shapefile " + config['shape_file_extra']
            if chunk_index is None:
                extra = ehdu[1][:]
            else:
                extra = ehdu[1][chunkindex*chunk_size : (chunkindex+1)*chunk_size]
            ehdu.close()
        except KeyError:
            pass
    else:
    # apply shape cuts: either on the file itself of on the extra file
    # since we're working with FITS type selections, we can't apply it
    # directly to the shapes array, but need to go back to the catalogs.
    # that's not really elegant since the .where runs on entire table
        cuts = " && ".join(config['shape_cuts'])
        try:
            ehdu = fitsio.FITS(config['shape_file_extra'])
            mask = ehdu[1].where(cuts)
            if verbose:
                print "opening extra shapefile " + config['shape_file_extra']
                print "selecting %d shapes" % mask.size
            if chunk_index is not None:
                mask = mask[chunk_index*chunk_size : (chunk_index+1)*chunk_size]
            shapes = shdu[1][mask]
            extra = ehdu[1][mask]
            ehdu.close()
        except KeyError:
            mask = shdu[1].where(cuts)
            if verbose:
                print "selecting %d shapes" % mask.size
            if chunk_index is not None:
                mask = mask[chunk_index*chunk_size : (chunk_index+1)*chunk_size]
            shapes = shdu[1][mask]
        del mask
    if verbose:
        print "shape sample (this chunk): %d" % shapes.size
    shdu.close()

    # if there's an extra file: join data with shapes
    if extra is not None:
        shapes = joinDataSets(shapes, extra)
    return shapes

def getLensCatalog(config, verbose=False):
    lensfile = config['lens_file']
    hdu = fitsio.FITS(lensfile)
    if verbose:
        print "opening lensfile %s (%d entries)" % (lensfile, hdu[1].get_nrows())
    mask = None
    if len(config['lens_cuts']) == 0:
        lenses = hdu[1][:]
    else:
        cuts = " && ".join(config['lens_cuts'])
        mask = hdu[1].where(cuts)
        if verbose:
            print "selecting %d lenses" % mask.size
        lenses = hdu[1][mask]
    hdu.close()
    if verbose:
        print "lens sample: %d" % lenses.size

    # see if there's an extra file
    try:
        hdu = fitsio.FITS(config['lens_extra_file'])
        if verbose:
            print "opening extra lensfile %s (%d entries)" % (config['lens_extra_file'], hdu[1].get_nrows())
        if mask is None:
            extra = hdu[1][:]
        else:
            extra = hdu[1][mask]
        hdu.close()
        lenses = joinDataSets(lenses, extra)
    except (KeyError, IOError) as exc: # not in config or file doesn't exist
        pass
    
    return lenses

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

# use actual LaTeX to render plot and fonts
from pylab import rcParams
def setTeXPlot(sampling=1):
    params = {
        'backend': 'ps',
        'ps.distiller.res': 6000,
        'axes.labelsize': sampling*9,
        'axes.linewidth' : sampling*0.25,
        'font.size': sampling*8,
        'text.fontsize': sampling*8,
        'legend.fontsize': sampling*8,
        'legend.markerscale' : sampling*0.5,
        'xtick.labelsize': sampling*8,
        'ytick.labelsize': sampling*8,
        'font.family': 'serif',
        'font.serif': 'Times',
        'font.weight': 'medium',
        'text.usetex': 'times',
        'figure.subplot.right' : 0.995,
        'figure.subplot.top' : 0.97,
        'figure.subplot.left' : 0.125,
        'figure.subplot.bottom' : 0.07,
    }
    rcParams.update(params)
