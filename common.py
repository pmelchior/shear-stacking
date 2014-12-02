import numpy as np
from math import pi, sqrt
import os

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
            return None
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
            return None
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
        for i in xrange(len(self.bins)-1):
            n[i] = self.Q[i].N
            mean_q[i] = self.Q[i].getMean()
            std_q[i] = self.Q[i].getStd()
        return np.array(self.R) / n, n, mean_q, std_q

# extrapolation function from
# http://stackoverflow.com/questions/2745329/how-to-make-scipy-interpolate-give-an-extrapolated-result-beyond-the-input-range
def extrap(x, xp, yp):
    """np.interp function with linear extrapolation"""
    y = np.interp(x, xp, yp)
    y[x < xp[0]] = yp[0] + (x[x<xp[0]]-xp[0]) * (yp[0]-yp[1]) / (xp[0]-xp[1])
    y[x > xp[-1]]= yp[-1] + (x[x>xp[-1]]-xp[-1])*(yp[-1]-yp[-2])/(xp[-1]-xp[-2])
    return y  


from galsim import Cosmology
cosmo = Cosmology()
def getBeta(z_c, z):
    if z_c >= z:
        return 0
    else:
        return cosmo.Da(z, z_c)/cosmo.Da(z)  

def getSigmaCrit(z_c, z):
    c2_4piG = 4. # in 1e14 M_solar / Mpc^2 (since cosmo.Da comes in units of c/H0)
    return c2_4piG / getBeta(z_c, z) / cosmo.Da(z_c)

def getSpecZCalibration():
    thisdir = os.path.dirname(os.path.realpath(__file__))
    return np.loadtxt(thisdir + '/data/checkphotoz_sv_deep_i24_psf1.2_sva1.dat')

def getSigmaCritCorrection(specz_calib, z_c):
    z_phot = 0.05 + 0.1*np.arange(20)
    above = z_phot > z_c
    z_phot = z_phot[above]
    cz = np.zeros_like(z_phot)
    for i in range(len(z_phot)):
        z_s = z_phot[i]
        mask = np.abs(specz_calib[:,0] - z_s) < 0.01
        SigmaCrit_ = getSigmaCrit(z_c, z_s)
        z_spec = specz_calib[:,1][mask]
        prob = specz_calib[:,2][mask]
        for j in xrange(len(z_spec)):
            if z_spec[j] > z_c:
                cz[i] += prob[j]*getSigmaCrit(z_c, z_spec[j])**-1
        cz[i] = cz[i]**-1
    return z_phot, cz

def getSigmaCritEffective(z_phot, cz, z):
    return extrap(z, z_phot, cz)

# radius of galaxy compared to PSF
def sizeRatio(data):
    return 2*data['RADIUS'] / (data['PSF_FWHM'] * 0.263)

# approximation of SNR of the size (instead of the flux)
def SNRadius(data):
    return data['SNR'] * data['RADIUS']

# Eli's modest S/G classification
def ModestSG(data):
    return invert(((data['CLASS_STAR_I'] > 0.3) & (data['MAG_AUTO_I'] < 18.0)) | ((data['SPREAD_MODEL_I'] + 3*data['SPREADERR_MODEL_I'] < 0.003) | ((data['MAG_PSF_I'] > 30.0) & (data['MAG_AUTO_I'] < 21.0)))) & (abs(data['SPREAD_MODEL_I']) < 0.1)  &  (data['FLAGS_I'] <= 3)

# bulge/disk flux ratio
def B_D(data, band='r'):
    return data['im3shape_' + band.lower() + '_bulge_flux'] / data['im3shape_' + band.lower() + '_disc_flux']

# SNR-dependent weights for im3shape
def SNR_weight(data, band='r'):
    return 0.2/(0.2**2 + (0.1*20/data['im3shape_' + band.lower() + '_snr'])**2)**0.5

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
