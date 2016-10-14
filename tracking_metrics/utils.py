import numpy as np
from surveyStatus import HealpixLookup
import healpy as hp
import lsst.sims.utils as utils
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy import units as u


def wrapRADec(ra, dec):
    # XXX--from MAF, should put in general utils
    """
    Wrap RA into 0-2pi and Dec into +/0 pi/2.

    Parameters
    ----------
    ra : numpy.ndarray
        RA in radians
    dec : numpy.ndarray
        Dec in radians

    Returns
    -------
    numpy.ndarray, numpy.ndarray
        Wrapped RA/Dec values, in radians.
    """
    # Wrap dec.
    low = np.where(dec < -np.pi / 2.0)[0]
    dec[low] = -1 * (np.pi + dec[low])
    ra[low] = ra[low] - np.pi
    high = np.where(dec > np.pi / 2.0)[0]
    dec[high] = np.pi - dec[high]
    ra[high] = ra[high] - np.pi
    # Wrap RA.
    ra = ra % (2.0 * np.pi)
    return ra, dec


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


def generate_taget_maps(nside=128):
    """
    Generate a suite of target depths for each filter.
    """
    npix = hp.nside2npix(nside)
    ra, dec = utils._hpid2RaDec(nside, np.arange(npix))
    coord = SkyCoord(ra=ra*u.rad, dec=dec*u.rad)
    g_long, g_lat = coord.galactic.l.radian, coord.galactic.b.radian
    eclip_long, eclip_lat = coord.barycentrictrueecliptic.lon.radian, coord.barycentrictrueecliptic.lat.radian

    # Set borders of different regions
    main_dec = 0

    south_dec = np.radians(-60.)

    nes_dec = 0.
    nes_eclip = np.radians(15.)

    gal_lat = np.radians(10.)
    gal_long_end1 = np.radians(70.)
    gal_long_end2 = np.radians(270.)

    wfd_region = np.where((dec <= main_dec) &
                          ((np.abs(g_lat) >= gal_lat) | (g_long > gal_long_end1)) &
                          (dec >= south_dec))
    south_pole_region = np.where(dec < south_dec)
    nes_region = np.where((dec > nes_dec) & (eclip_lat < nes_eclip))

    # Need to add in the galactic plane taper.
    galatic_plane_region = np.where((np.abs(g_lat) < gal_lat) &
                                    ((g_long < gal_long_end1) | (g_long > gal_long_end2)))

    wfd_depth = {'u': 26.1, 'g': 27.4, 'r': 27.5, 'i': 26.8, 'z': 26.1, 'y': 24.9}
    nes_depth = {'g': 26.9, 'r': 27.0, 'i': 26.3, 'z': 23.6}
    south_pole_depth = {'u': 25.6, 'g': 26.9, 'r': 27.0, 'i': 26.3, 'z': 23.6, 'y': 24.4}
    galactic_plane_depth = {'u': 25.6, 'g': 26.9, 'r': 27.0, 'i': 26.3, 'z': 23.6, 'y': 24.4}

    all_filters = []
    all_filters.extend(wfd_depth.keys())
    all_filters.extend(nes_depth.keys())
    all_filters.extend(south_pole_depth.keys())
    all_filters.extend(galactic_plane_depth.keys())
    all_filters = list(set(all_filters))

    all_depth = {}
    for filtername in all_filters:
        all_depth[filtername] = np.empty(npix, dtype=float)
        all_depth[filtername].fill(hp.UNSEEN)

    # Should probably throw an error if any of the regions overlap.

    for filtername in wfd_depth:
        all_depth[filtername][wfd_region] = wfd_depth[filtername]

    for filtername in south_pole_depth:
        all_depth[filtername][south_pole_region] = south_pole_depth[filtername]

    for filtername in nes_depth:
        all_depth[filtername][nes_region] = nes_depth[filtername]

    for filtername in galactic_plane_depth:
        all_depth[filtername][galatic_plane_region] = galactic_plane_depth[filtername]

    return all_depth


class min_coadd_power_tesselate(object):
    """
    Take a block of tesselated pointings and rotate them to a new spot such that
    the resulting power in the coadded depth map is minimized.
    """
    def __init__(self, nside=128):
        field_data = np.loadtxt('fieldID.dat', delimiter='|', skiprows=1,
                                dtype=zip(['id', 'ra', 'dec'], [int, float, float]))

        ra_range = 15.  # Degrees
        dec_range = 45.
        good = np.where(((field_data['ra'] <= ra_range) | (field_data['ra'] >= 360.-ra_range)) &
                        ((field_data['dec'] >= -dec_range) & (field_data['dec'] <= dec_range)))

        self.field_data = field_data[good]
        self.field_data['ra'] = np.radians(self.field_data['ra'])
        self.field_data['dec'] = np.radians(self.field_data['dec'])
        self.nside = nside
        # put a kdtree in here
        self.kdtree = HealpixLookup(nside=nside)

    def find_pointings(self, reward_map, m5_map, single_visit_depth=22.):
        """
        Paramters
        ---------
        reward_map : np.array
            A healpix map where unmasked pixels should be tesselated with pointings.
        m5_map : np.array
            A healpix map of the co-added 5-sigma limiting depth at each point for the
            current status of the survey.
        single_visit_depth : float
            The rough depth of a single visit.
        """

        # Check that target map is correct nside

        # Find the peak in the target_map
        # Maybe just use the peak pixel and it's 8 neighbors and go with that?
        max_index = np.where(reward_map == np.max(reward_map))[0].max()
        peak_pixels = hp.get_all_neighbours(max_index)
        pass
        #target_ra, target_dec = 