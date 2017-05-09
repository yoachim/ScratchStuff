from __future__ import print_function
from lsst.sims.maf.utils import obs2sqlite
import numpy as np


filename = 'sim_results.csv'
names = ['mjd', 'ra', 'dec', 'filter']
types = [float, float, float, '|S1']
data = np.loadtxt(filename, skiprows=1, delimiter=',', dtype=zip(names, types))

# make a rough unix to MJD conversion
data['mjd'] = data['mjd']/86400. + 40587
print('Converting %i obserations to sqlite' %np.size(data))

obs2sqlite(data)

