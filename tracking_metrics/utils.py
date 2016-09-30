import numpy as np
from surveyStatus import HealpixLookup
import healpy as hp

def generate_taget_maps(nside=128):
    """
    Generate a suite of target depths for each filter
    """
    target_median_depths = {'u': 26.1, 'g': 27.4, 'r': 27.5, 'i': 26.8, 'z': 26.1, 'y': 24.9}
    nes_depth = {}
    south_pole_depth = {}
    galactic_plane_depth = {}

    all_filters = []
    all_filters.extend(target_map.keys())
    all_fitlers.extend(nes_depth.keys())
    all_fitlers.extend(south_pole_depth.keys())
    all_filters.extend(galactic_plane_depth.keys())
    all_filters = list(set(all_filters))


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
            A healpix map of the co-added 5-sigma limiting depth at each point.
        single_visit_depth : float
            The rough depth of a single visit.
        """

        # Check that target map is correct nside

        # Find the peak in the target_map
        # Maybe just use the peak pixel and it's 8 neighbors and go with that?
        max_index = np.where(reward_map == np.max(reward_map))[0].max()
        peak_pixels = hp.get_all_neighbours(max_index)
        
        target_ra, target_dec = 