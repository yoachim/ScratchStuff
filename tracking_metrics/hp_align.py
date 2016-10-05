import numpy as np
from utils import wrapRADec


def rotate_ra_dec(ra, dec, ra_target, dec_target, init_rotate=0.):
    """
    Rotate ra and dec coordinates to be centered on a new dec.

    Inputs
    ------
    ra : float or np.array
        RA coordinate(s) to be rotated in radians
    dec : float or np.array
        Dec coordinate(s) to be rotated in radians
    ra_rotation : float
        RA distance to rotate in radians
    dec_target : float
        Dec distance to rotate in radians
    init_rotate : float (0.)
        The amount to rotate the points around the x-axis first (radians).
    """
    # point (ra,dec) = (0,0) is at x,y,z = 1,0,0

    x, y, z = ra_dec_2_xyz(ra, dec)

    # Rotate around the x axis to start
    c_i = np.cos(init_rotate)
    s_i = np.sin(init_rotate)
    xp = x
    yp = c_i*y - s_i*z
    zp = s_i*y + c_i*z

    theta_y = dec_target
    c_ty = np.cos(theta_y)
    s_ty = np.sin(theta_y)

    # Rotate about y
    xp2 = c_ty*xp + s_ty*zp
    zp2 = -s_ty*xp + c_ty*zp

    # Convert back to RA, Dec
    ra_p = np.arctan2(yp, xp2)
    dec_p = np.arcsin(zp2)

    # Rotate to the correct RA
    ra_p += ra_target

    ra_p, dec_p = wrapRADec(ra_p, dec_p)

    return ra_p, dec_p


def ra_dec_2_xyz(ra, dec):
        """Calculate x/y/z values for ra/dec points, ra/dec in radians."""
        # Note ra/dec can be arrays.
        x = np.cos(dec) * np.cos(ra)
        y = np.cos(dec) * np.sin(ra)
        z = np.sin(dec)
        return x, y, z
