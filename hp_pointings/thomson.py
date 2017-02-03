import numpy as np
import healpy as hp
from scipy.optimize import minimize
from fib_sphere_grid import fib_sphere_grid
from healpyUtils import *
import matplotlib.pylab as plt

# Would like to be able to get N points evenly distributed on a sphere.

# Note from here: http://www.mcs.anl.gov/~zippy/publications/cgoprl/node3.html
# This may not actually solve what I want it to do.


def ra_wrap(ra_in):

    ra = ra_in % 2.*np.pi
    ra[np.where(ra < 0)] += 2.*np.pi
    return ra


def radec2xyz(ra, dec):
    x = np.cos(dec)*np.sin(ra)
    y = np.cos(dec)*np.cos(ra)
    z = np.sin(dec)
    return x, y, z


def xyz2radec(xx, x=None, y=None, z=None):
    """
    Given (x,y,z), return the ra and dec of the point on the unit sphere
    """

    if xx is not None:
        if np.ndim(xx) == 1:
            xx = xx.reshape(3, xx.size/3)


        x = xx[0]
        y = xx[1]
        z = xx[2]

    # normalize by the radius.
    norm = (x**2+y**2+z**2)**0.5
    x = x/norm
    y = y/norm
    z = z/norm

    dec = np.arcsin(z)
    ra = np.arctan2(y, x) + np.pi
    return ra, dec


def limit_wrap(ra_in, dec_in):
    dec = dec_in.copy()
    ra = ra_in.copy()
    while dec.max() > np.pi/2.:
        over = np.where(dec > np.pi/2)
        dec[over] -= np.pi/2
        ra[over] += np.pi
    while dec.min() < -np.pi/2.:
        over = np.where(dec < -np.pi/2)
        dec[over] += np.pi/2
        ra[over] += np.pi

    ra = ra % 2.*np.pi
    while ra.min() < 0:
        ra[np.where(ra < 0)] += 2.*np.pi

    return ra, dec


def points_on_sphere_potential(xx):
    if xx.ndim == 1:
        # x = x.reshape(x.size/2, 2)
        xx = xx.reshape(xx.size/3, 3)
    x = xx[:, 0]
    y = xx[:, 1]
    z = xx[:, 2]

    norm = (x**2+y**2+z**2)**0.5
    x = x/norm
    y = y/norm
    z = z/norm

    di = np.diag_indices(x.size)
    r_sq = (x[:,np.newaxis] - x)**2 + (y[:,np.newaxis] - y)**2 +(z[:,np.newaxis] - z)**2 

    r_sq[di] = 1.

    potential = 1./r_sq
    potential[di] = 0.

    if True in np.isinf(potential):
        import pdb ; pdb.set_trace()

    return np.sum(potential)



def blah_points_on_sphere_potential(x):
    """
    return the potential energy for n points on a sphere
    x should be x,y,z. Projects them down radius vector to unit sphere.

    Got equation from https://www.cs.purdue.edu/homes/pengh/reports/590OP.pdf
    """

    if x.ndim == 1:
        # x = x.reshape(x.size/2, 2)
        x = x.reshape(x.size/3, 3)

    ra, dec = xyz2radec(None, x=x[:, 0], y=x[:, 1], z=x[:, 2])

    #ra = x[:, 0]
    #dec = x[:, 1]

    cos_phi = np.cos(dec)
    sin_phi = np.sin(dec)

    di = np.diag_indices(cos_phi.size)

    cos_theta_diff = np.cos(ra[:, np.newaxis] - ra)
    denom = cos_phi[:, np.newaxis]*cos_phi
    denom *= cos_theta_diff
    denom += sin_phi[:, np.newaxis]*sin_phi
    denom -= 1.
    denom[di] = 1.

    func = 1./denom

    # Make sure the diagonal is zero
    func[di] = 0

    if np.inf in func:
        import pdb ; pdb.set_trace()

    return np.abs(np.sum(func))


def temp_points_on_sphere_potential(x):
    """
    x should be ra, dec pairs
    """
    di = np.diag_indices(x[:, 0].size)
    phi = np.pi/2. - x[:, 1]
    cos_dec_diff = np.cos(phi[:, np.newaxis]-phi)
    sin_prod = np.sin(x[:, 0, np.newaxis])*np.sin(x[:, 0])
    cos_prod = np.cos(x[:, 0, np.newaxis]) * np.cos(x[:, 0])
    energy = 1./(1.-2.*(sin_prod*cos_dec_diff+cos_prod))
    energy[di] = 0

    return np.sum(energy)


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


def haver_potential(xx):
    """
    Find the potential energy if force is along radial distance
    """
    if xx.ndim == 1:
        # x = x.reshape(x.size/2, 2)
        xx = xx.reshape(xx.size/3, 3)
    x = xx[:, 0]
    y = xx[:, 1]
    z = xx[:, 2]

    #ra = x[:, 0]
    #dec = x[:, 1]
    #ra, dec = limit_wrap(ra, dec)

    ra, dec = xyz2radec(None, x=x, y=y, z=z)

    di = np.diag_indices(ra.size)
    # There's a factor of 2 savings to be had since i,j and j,i are symetric.
    distances = haversine(ra[:, np.newaxis], dec[:, np.newaxis], ra, dec)
    distances[di] = 1.
    potential = 1./distances
    potential[di] = 0
    return np.abs(np.sum(potential))


if __name__ == "__main__":

    pts = np.array([[0, 0, 1], [0, -1, 0]])

    e1 = points_on_sphere_potential(pts)

    he1 = haver_potential(pts)

    pts = np.array([[0, 0, 1], [0., 0, -1.]])
    e2 = points_on_sphere_potential(pts)
    he2 = haver_potential(pts)

    print 'e1 = %f, e2 = %f' % (e1, e2)
    print 'he1 = %f, he2 = %f' % (he1, he2)

    lons = np.radians(np.arange(45., 405, 90.))
    lats = np.radians(np.array([45., -45.]))

    cube_points = np.zeros((2, 8), dtype=float)
    i = 0
    for lon in lons:
        for lat in lats:
            cube_points[0, i] = lon
            cube_points[1, i] = lat
            i += 1

    cx, cy, cz = radec2xyz(cube_points[0,:], cube_points[1,:])
    cube_points = np.vstack((cx, cy, cz)).T
    print 'cube potential = %f' % points_on_sphere_potential(cube_points)
    print 'cube potential_angular = %f' % haver_potential(cube_points)

    # I'd like to make a solver to show that minimizing the regular potential makes something strange, and 
    # minimzing the angular distance potential makes a cube when given 8 points.

    haver_opt = minimize(haver_potential, cube_points)
    reg_opt = minimize(points_on_sphere_potential, cube_points)

    import pdb ; pdb.set_trace()

    ra, dec = xyz2radec(cube_points.T)
    hp.mollview(_healbin(ra, dec, dec*0+1, nside=16, reduceFunc=np.sum), title='input')

    ra, dec = xyz2radec(haver_opt.x)
    hp.mollview(_healbin(ra, dec, dec*0+1, nside=16), title='haver_opt')

    ra, dec = xyz2radec(reg_opt.x)
    hp.mollview(_healbin(ra, dec, dec*0+1, nside=16), title='reg_opt')

    plt.show()



    # Let's try solving a 200 point grid
    ra, dec = fib_sphere_grid(20)
    x, y, z = radec2xyz(ra, dec)
    # fib_points = np.hstack((ra,dec))
    fib_points = np.hstack((x, y, z))

    haver_opt = minimize(haver_potential, fib_points)
    reg_opt = minimize(points_on_sphere_potential, fib_points)
    

    print 'haver potential = %f' % haver_potential(haver_opt.x)
    print 'reg potential = %f' % points_on_sphere_potential(reg_opt.x)

    start = _healbin(ra, dec, ra*0+1, nside=16)
    ra, dec = xyz2radec(haver_opt.x)
    end = _healbin(ra, dec, ra*0+1, nside=16)


    # Ugh, minimizing on a sphere is a pain in the ass. maybe I should mimize where I vary x,y,z? Then just use those to compute theta and phi
    # That way I don't have to worry about angle wrapping things, points near the poles will be well behaved. 


    # OK, this shit is failing miserably. I need to work out given an array of ra,dec find the physical distances and angular distances
    # between all the pairs of points. 
    # Generate angular distances between list of points. Use numpy.triu_indices to get the i,j distances.
    # Make arrays of ra_i, ra_j, dec_i, dec_j. 
    # Ang dist to physical dist. 
    # Potential given distance matrix. 



