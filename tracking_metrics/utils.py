import numpy as np
from surveyStatus import HealpixLookup
import healpy as hp

def generate_taget_maps():
    """
    Generate a suite of target depths for each filter
    """

    target_median_depths = {'u':}


class simple_tesselate(object):
    """
    Take a block of tesselated pointings and rotate them to a new spot
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


    def find_pointings(self, reward_map, m5_map):
        """
        Paramters
        ---------
        target_map : np.array
            A healpix map where unmasked pixels should be tesselated with pointings
        """

        # Check that target map is correct nside

        # Find the peak in the target_map
        # Maybe just use the peak pixel and it's 8 neighbors and go with that?
        max_index = np.where(reward_map == np.max(reward_map))[0].max()
        peak_pixels = hp.get_all_neighbours(max_index)
        
        target_ra, target_dec = 