import numpy as np
import lsst.sims.maf.db as db
import lsst.sims.skybrightness as sb
import matplotlib.pylab as plt

# Let's see how much my sky values differ from the current opsim values

stride = 100

opsdb = db.OpsimDatabase('ewok_1004_sqlite.db')
sql = ''
data = opsdb.fetchMetricData(['expMJD','filter','fieldRA','fieldDec','filtSkyBrightness'], sql)
data = data[::stride]


mysky = np.zeros(data.size, dtype=float)

sm = sb.SkyModel(mags=True)
ind2filt = {'u':0, 'g':1,'r':2,'i':3,'z':4,'y':5}

for i,dat in enumerate(data):
    sm.setRaDecMjd([dat['fieldRA']], [dat['fieldDec']], dat['expMJD'])
    mysky[i] = sm.returnMags()[0][ind2filt[dat['filter']]]


fig,axes = plt.subplots(2,3)
axes = axes.ravel()
keys = ['u','g','r','i','z','y']

print 'filter, median OpSim, median ESO'
for i,key in enumerate(keys):
    good = np.where(data['filter'] == key)
    axes[i].plot(data['filtSkyBrightness'][good],data['filtSkyBrightness'][good]-mysky[good], 'ko', alpha=.1)
    axes[i].set_title(r'$'+key+'$')
    axes[i].set_ylabel(r'$\mu_{OpSim}-\mu_{ESO}$')
    axes[i].set_xlabel(r'$\mu_{OpSim}$')
    print key, np.median(data['filtSkyBrightness'][good]), np.median(mysky[good])

fig.tight_layout()
fig.savefig('skycompare.pdf')
