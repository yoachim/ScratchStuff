import numpy as np
import healpy as hp

# Would like to be able to get N points evenly distributed on a sphere.

# Note from here: http://www.mcs.anl.gov/~zippy/publications/cgoprl/node3.html
# This may not actually solve what I want it to do.

def points_on_sphere_potential(x):
    """
    return the potential energy for n points on a sphere
    x should be theta, phi pairs in radians

    Got equation from https://www.cs.purdue.edu/homes/pengh/reports/590OP.pdf

    XXX--need to get this in dec rather than phi I think.
    """

    # Arg, need to make nxn arrays, set diagonal to zero, and sum to avoid loop

    cos_phi = np.cos(x[:, 1])
    sin_phi = np.sin(x[:, 1])

    di = np.diag_indices(cos_phi.size)

    cos_theta_diff = np.cos(x[:, 0, np.newaxis] - x[:, 0])
    denom = cos_phi[:, np.newaxis]*cos_phi
    denom *= cos_theta_diff
    denom += sin_phi[:, np.newaxis]*sin_phi
    denom -= 1.
    denom[di] = 1.

    func = 1./denom

    # Make sure the diagonal is zero
    func[di] = 0


    return np.abs(np.sum(func))

# Should see what happens if the force is based on angular distance rather than 3-D distance.


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


def haver_potential(x):
    """
    Find the potential energy if force is along radial distance
    """
    di = np.diag_indices(x[:, 0].size)
    # There's a factor of 2 savings to be had since i,j and j,i are symetric.
    distances = haversine(x[:, 0, np.newaxis], x[:, 1, np.newaxis], x[:, 0], x[:, 1])

    potential = 1./distances
    potential[di] = 0
    return np.abs(np.sum(potential))


if __name__ == "__main__":


    pts = np.array([[0,0],[0,np.pi]])
    #e1 = points_on_sphere_potential(pts)

    #he1 = haver_potential(pts)

    pts = np.array([[0,0],[np.pi/2,0]])
    #e2 = points_on_sphere_potential(pts)
    #he2 = haver_potential(pts)

    #print 'e1 = %f, e2 = %f' %(e1, e2)
    #print 'he1 = %f, he2 = %f' %(he1, he2)

    lons = np.radians(np.arange(45., 405, 90.))
    lats = np.radians(np.array([45., -45.]))

    cube_points = np.zeros((8,2), dtype=float)
    i = 0
    for lon in lons:
        for lat in lats:
            cube_points[i,0] = lon
            cube_points[i,1] = lat
            i += 1

    print 'cube potential = %f' % points_on_sphere_potential(cube_points)
    print 'cube potential_angular = %f' % haver_potential(cube_points)

