import numpy as np
import healpy as hp
import matplotlib.pylab as plt
from utils import rotate_ra_dec
from surveyStatus import countFilterStatus, HealpixLookup


field_data = np.loadtxt('fieldID.dat', delimiter='|', skiprows=1,
                        dtype=list(zip(['id', 'ra', 'dec'], [int, float, float])))

ra_range = 15.  # Degrees
dec_range = 25.
good = np.where(((field_data['ra'] <= ra_range) | (field_data['ra'] >= 360.-ra_range)) &
                ((field_data['dec'] >= -dec_range) & (field_data['dec'] <= dec_range)))

field_data = field_data[good]

new_ra, new_dec = rotate_ra_dec(np.radians(field_data['ra']), np.radians(field_data['dec']),
                                np.radians(0.), np.radians(-45.), init_rotate=np.radians(0))

nside = 128
hpl = HealpixLookup(nside=nside)
ssl = []
ssl.append(countFilterStatus(filter_name=['r'], nside=nside))


class obs(object):
    """
    Simple observation object
    """
    def __init__(self):
        self.filter = 'r'
        self.ra = 0.
        self.dec = 0

observation = obs()
for ra, dec in zip(new_ra, new_dec):
    observation.filter = 'r'
    observation.ra = ra
    observation.dec = dec
    pix = hpl.lookup(ra, dec)
    ssl[0].add_visit(observation, pix)


hp.mollview(ssl[0].survey_map)
plt.show()
