import numpy as np
from scipy.optimize import minimize
# One more time trying to solve the Thomson problem with python since I'm a glutton for punishment


def potential(x0):
    """
    Assume x0 is theta and phi values
    """

    theta = x0[0:x0.size/2]
    phi = x0[x0.size/2:]

    x = np.sin(phi)*np.cos(theta)
    y = np.sin(phi)*np.sin(theta)
    z = np.cos(phi)

    dsq = 0.

    indices = np.triu_indices(x.size, k=1)

    for coord in [x, y, z]:
        coord_i = np.tile(coord, (coord.size, 1))
        coord_j = coord_i.T
        d = (coord_i[indices]-coord_j[indices])**2
        dsq += d

    U = np.sum(1./np.sqrt(dsq))
    return U


# let's generate some random points on a sphere
npts = 100
theta = np.random.rand(npts)*np.pi*2.
phi = np.arccos(2.*np.random.rand(npts)-1.)

x = np.concatenate((theta, phi))

ack = minimize(potential, x, method='CG')
