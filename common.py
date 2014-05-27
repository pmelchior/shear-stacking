from numpy import *
from math import pi, sqrt

def skyAngle(ra, dec, ra_ref, dec_ref):
    # CAUTION: this needs to be a pseudo-Cartesian coordinate frame
    # (not pure RA/DEC), otherwise angles are skewed
    return arctan2(dec-dec_ref, (ra-ra_ref)*cos(dec*pi/180))

def skyDistance(ra, dec, ra_ref, dec_ref):
    # CAUTION: this needs to be a pseudo-Cartesian coordinate frame
    # (not pure RA/DEC), otherwise distances are skewed
    return (((ra-ra_ref)*cos(dec*pi/180))**2 + (dec-dec_ref)**2)**0.5

def tangentialShear(ra, dec, e1, e2, ra_ref, dec_ref, computeB=False):
    # RA/Dec is a left-handed coordinate frame
    # but it needs to be right-handed, hence the minus sign
    phi = skyAngle(ra, dec, ra_ref, dec_ref)
    if computeB is False:
        return -e1*cos(2*phi) - e2*sin(2*phi)
    else:
        return -e1*cos(2*phi) - e2*sin(2*phi), e1*sin(2*phi) - e2*cos(2*phi)

def scalarProfile(bins, r, q, weight):
    mean_q = zeros(len(bins)-1)
    std_q = zeros(len(bins)-1)
    mean_r = zeros(len(bins)-1)
    n = zeros(len(bins)-1, dtype='int64')
    for i in range(len(bins)-1):
        mask = (r >= bins[i]) & (r < bins[i+1])
        n[i] = sum(mask)
        mean_r[i] = r[mask].mean()
        V1 = weight[mask].sum()
        mean_q[i] = (q[mask]*weight[mask]).sum()/V1
        # CAUTION: assumes Gaussian errors and large samples
        # replace with Jackknife/Bootstrap estimate!
        std_q[i] = ((weight[mask]*(q[mask]-mean_q[i])**2).sum()/V1)**0.5/sqrt(sum(mask))
    return mean_r, n, mean_q, std_q

# extrapolation function from
# http://stackoverflow.com/questions/2745329/how-to-make-scipy-interpolate-give-an-extrapolated-result-beyond-the-input-range
def extrap(x, xp, yp):
    """np.interp function with linear extrapolation"""
    y = interp(x, xp, yp)
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
    return loadtxt('data/checkphotoz_sv_deep_i24_psf1.2_sva1.dat')

def getSigmaCritCorrection(specz_calib, z_c):
    z_phot = 0.05 + 0.1*arange(20)
    above = z_phot > z_c
    z_phot = z_phot[above]
    cz = zeros_like(z_phot)
    for i in range(len(z_phot)):
        z_s = z_phot[i]
        mask = abs(specz_calib[:,0] - z_s) < 0.01
        SigmaCrit_ = getSigmaCrit(z_c, z_s)
        z_spec = specz_calib[:,1][mask]
        prob = specz_calib[:,2][mask]
        for j in range(len(z_spec)):
            if z_spec[j] > z_c:
                cz[i] += prob[j]*getSigmaCrit(z_c, z_spec[j])**-1
        cz[i] = cz[i]**-1
    return z_phot, cz

def getSigmaCritEffective(z_phot, cz, z):
    return extrap(z, z_phot, cz)

