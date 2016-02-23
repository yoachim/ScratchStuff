import numpy as np
from medBD import medDB, single_frame
import healpy as hp
import matplotlib.pylab as plt
from lsst.sims.utils import Site
import ephem


nside = 32
median_map = np.zeros(hp.nside2npix(nside), dtype=float) + hp.UNSEEN

umjd = medDB(full_select='select DISTINCT(mjd) from medskybrightness;', dtypes=float)
site = Site('LSST')
sun = ephem.Sun()
moon = ephem.Moon()

obs = ephem.Observer()
obs.lat, obs.lon, obs.elevation = site.latitude_rad, site.longitude_rad,site.height
doff = ephem.Date(0)-ephem.Date('1858/11/17')
udjd = umjd - doff
moonAlts = umjd*0.
sunAlts = umjd*0.

for i,mjd in enumerate(umjd):
    obs.date = udjd[i]
    moon.compute(obs)
    sun.compute(obs)
    moonAlts[i] = moon.alt
    sunAlts[i] = sun.alt

goodDates = umjd[np.where((moonAlts < 0) & (sunAlts < np.radians(-18.)))]

for i in np.arange(hp.nside2npix(nside)):
    data = medDB(where_clause='hpindex = %i' % i)
    data = data[np.in1d(data['mjd'], goodDates)]
    if np.size(data) > 0:
        median_map[i] = np.median(data['R'])

np.savez('median_map.npz', median_map=median_map)
hp.mollview(median_map)
plt.show()
