import numpy as np
import healpy as hp


# Thinking about how one could dynamically place LSST pointings to minimize the power spectrum.

# Want to avoid the psudo-deep drilling fields as much as possible.

# given some selected region of healpixels, calc the min and max ra and dec.
# Select fields that have the same range but centered on

# Ah, maybe just do an evenish tesselation on the range dec +/- 45 degrees and RA +/- 15 degrees?
# Then rotate that tesselation to the center of the blob I want to do?

# rotate pointings to be centered on region center. Compute dict and reverse dict of healpixel to pointings?


# Maybe generate a grid of starting tesselations, and then rotate all of them to the new coordinates?

field_data = np.loadtxt('fieldID.dat', delimiter='|', skiprows=1,
                        dtype=zip(['id', 'ra', 'dec'], [int, float, float]))

ra_range = 15.  # Degrees
dec_range = 45.
good = np.where(((field_data['ra'] <= ra_range) | (field_data['ra'] >= 360.-ra_range)) &
                ((field_data['dec'] >= -dec_range) & (field_data['dec'] <= dec_range)))

field_data = field_data[good]

def rotate_around_ra_dec(ra, dec, ra_rotation, dec_rotaion):
    """
    Rotate Ra,dec points around a specified axis.

    handy sheet here:
    http://inside.mines.edu/fs_home/gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.pdf
    """





def rotate_ra_dec(ra, dec, ra_rotation, dec_rotation):
    """
    Rotate ra and dec coordinates to be centered on a new ra and dec.
    Coords are rotated first in dec, then in RA.

    Inputs
    ------
    ra : float or np.array
        RA coordinate(s) to be rotated in radians
    dec : float or np.array
        Dec coordinate(s) to be rotated in radians
    ra_rotation : float
        RA distance to rotate in radians
    dec_rotation : float
        Dec distance to rotate in radians
    """

    # OK, rotating to the correct ra is a rotation around the z axis
    # then rotating to the correct dec is a y rotation)
    # want to do the y-rotation first

    # proint (ra,dec) = (0,0) is at x,y,z = 1,0,0


    # XXXXX--wait, can't I just rotate the points around y-axis to get them all
    # on the correct decs, then add the ra-offset?

    x, y, z = ra_dec_2_xyx(ra, dec)

    theta_y = dec_rotation
    theta_z = ra_rotation

    c_ty = np.cos(theta_y)
    s_ty = np.sin(theta_y)
    # c_tz = np.cos(theta_z)
    # s_tz = np.sin(theta_z)

    # Rotate about y
    xp = c_ty*x + s_ty*z
    yp = y
    zp = -s_ty*x + c_ty*z

    # rotate about z
    # xp = c_tz*xp - s_tz*yp
    # yp = s_tz*xp + c_tz*yp

    #xp = x*c_ty*c_tz - y*s_tz + z*c_tz*s_ty
    #yp = s_tz*c_ty*x + y*c_tz + z*s_ty*s_tz
    #zp = -x*s_ty + z*c_ty

    ra_p = np.arctan2(yp, xp)
    dec_p = np.arcsin(zp)

    ra_p += ra_rotation


    return ra_p, dec_p


def ra_dec_2_xyx(ra, dec):
        """Calculate x/y/z values for ra/dec points, ra/dec in radians."""
        # Note ra/dec can be arrays.
        x = np.cos(dec) * np.cos(ra)
        y = np.cos(dec) * np.sin(ra)
        z = np.sin(dec)
        return x, y, z

# genenate a bunch of fields. make arrays that shift them by d_ra, d_dec, d_theta.  
# Fit the peak to find the best d_ra, d_dec, d_theta.
