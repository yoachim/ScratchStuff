from __future__ import print_function
from lsst.sims.maf.utils.obs2sqlite import obs2sqlite
import numpy as np
import pandas as pd
import sqlite3


input_db = 'MAFFBDE.db'
conn = sqlite3.connect(input_db)
data = pd.read_sql('select * from Schedule', conn, index_col='obsID')


# Do a quick remap to my favorite column names
names = ['filter', 'ra', 'dec', 'mjd', 'alt', 'az', 'night', 'sunAlt', 'moonAlt']
types = ['|S1']
types.extend([float]*(len(names)-1))

observations = np.zeros(data.shape[0], dtype=list(zip(names, types)))

convert_dict = {'filter': 'filter', 'ra': 'FieldRA', 'dec': 'FieldDec',
'mjd': 'MJD', 'alt': 'altitude', 'az': 'azimuth',
'night':'nightID',  'sunAlt':'sunAlt', 'moonAlt': 'moonAlt'}

for key in convert_dict:
    observations[key] = data[convert_dict[key]]

print('Converting %i obserations to sqlite' %np.size(data))

obs2sqlite(observations, radians=False, full_sky=True)
