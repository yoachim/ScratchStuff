import numpy as np
import healpy as hp
from scipy.optimize import minimize
from healpyUtils import _healbin
import matplotlib.pylab as plt
from iminuit import Minuit, describe
from iminuit.util import make_func_code


def haversine(long1, lat1, long2, lat2):
    """
    Return the angular distance between two points in radians
    @param [in] long1 is the longitude of point 1 in radians
    @param [in] lat1 is the latitude of point 1 in radians
    @param [in] long2 is the longitude of point 2 in radians
    @param [in] lat2 is the latitude of point 2 in radians
    @param [out] the angular separation between points 1 and 2 in radians
    From http://en.wikipedia.org/wiki/Haversine_formula
    """
    t1 = np.sin(lat2 / 2.0 - lat1 / 2.0)**2
    t2 = np.cos(lat1) * np.cos(lat2) * np.sin(long2 / 2.0 - long1 / 2.0)**2
    return 2 * np.arcsin(np.sqrt(t1 + t2))


def make_ang_ij(ra, dec):
    """
    given array of xx, that has ra, dec
    make arrays that are ra_i, ra_j
    """

    # ra = xx[:xx.size/2]
    # dec = xx[xx.size/2:]
    ra = np.array(ra)
    dec = np.array(dec)

    z = np.zeros((ra.size, ra.size), dtype=float)
    ra_i = ra + z
    ra_j = ra[:, np.newaxis] + z

    dec_i = dec + z
    dec_j = dec[:, np.newaxis] + z

    return ra_i, ra_j, dec_i, dec_j


def electron_sphere_potential(ra, dec):
    """
    find the potential of electrons on a sphere
    xx = [ra, dec] in radians
    """

    ra_i, ra_j, dec_i, dec_j = make_ang_ij(ra, dec)

    # Only need the upper triabgle of the matrix
    tri_indices = np.triu_indices(ra_i.shape[0], k=1)

    # Find the angle between the ij points
    distance = haversine(ra_i[tri_indices], dec_i[tri_indices], ra_j[tri_indices], dec_j[tri_indices])

    # Convert to physical distance
    distance = 2.*np.sin(distance/2.)

    # potential
    potential = np.sum(1./distance)
    return potential


# Need to class it up so minuit can call it
class potential_functor:
    def __init__(self, ra, dec):
        self.ra = ra
        self.dec = dec
        f_sig = describe(electron_sphere_potential)
        self.func_code = make_func_code(f_sig)  # docking off independent variable
        self.func_defaults = None  # this keeps np.vectorize happy

    def __call__(self, *arg):
        # ra = np.array(arg[0:len(arg)/2], dtype=float)
        # dec = np.array(arg[len(arg/2):], dtype=float)
        import pdb ; pdb.set_trace()
        return electron_sphere_potential(arg[0:len(arg)/2], arg[len(arg)/2:])



def bin_xx(ra, dec, nside=16):
    # ra = xx[:np.size(xx)/2]
    # dec = xx[np.size(xx)/2:]

    # Make sure ra wraps.
    ra = ra % 2.*np.pi
    neg = np.where(ra < 0)
    ra[neg] += 2.*np.pi
    dec = np.arcsin(np.sin(dec))

    binned_map = _healbin(ra, dec, dec*0+1, nside=nside)
    return binned_map


if __name__ == "__main__":

    ra = np.array([0., np.pi/2.])
    dec = np.array([0.,0.])
    print electron_sphere_potential(ra, dec)

    npts = 100
    ra = np.random.rand(npts)*2.*np.pi
    dec = np.random.rand(npts)*np.pi - np.pi/2.
    xx = np.hstack((ra, dec))

    pf = potential_functor(ra, dec)
    pf()

    # fit = minimize(electron_sphere_potential, xx)
    m = Minuit(pf, ra=ra, dec=dec, print_level=1)
    m.migrad()



