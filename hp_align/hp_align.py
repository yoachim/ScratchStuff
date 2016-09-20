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


def rotate_ra_dec(ra, dec, ra_rotation, dec_rotation):
    """
    rotate some x,y,z positions
    """

    # OK, rotating to the correct ra is a rotation around the z axis
    # then rotating to the correct dec is a y rotation)
    # want to do the y-rotation first

    # proint (ra,dec) = (0,0) is at x,y,z = 1,0,0

    x, y, z = _treexyz(ra, dec)

    theta_y = dec_rotation
    theta_z = ra_rotation

    c_ty = np.cos(theta_y)
    s_ty = np.sin(theta_y)
    c_tz = np.cos(theta_z)
    s_tz = np.sin(theta_z)

    xp = x*c_ty*c_tz - y*s_tz + z*c_tz*s_ty
    yp = s_tz*c_ty*x + y*c_tz + z*s_ty*s_tz
    zp = -x*s_ty + z*c_ty

    ra_p = np.arctan2(yp, xp)
    dec_p = np.arcsin(zp)

    return ra_p, dec_p


def _treexyz(ra, dec):
        """Calculate x/y/z values for ra/dec points, ra/dec in radians."""
        # Note ra/dec can be arrays.
        x = np.cos(dec) * np.cos(ra)
        y = np.cos(dec) * np.sin(ra)
        z = np.sin(dec)
        return x, y, z
