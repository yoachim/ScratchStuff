import numpy as np
import lsst.sims.maf.db as db
import lsst.sims.skybrightness as sb
import matplotlib.pylab as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter

# Let's see how much my sky values differ from the current opsim values

stride = 100

opsdb = db.OpsimDatabase('ewok_1004_sqlite.db')
sql = ''
data = opsdb.fetchMetricData(['expMJD','filter','fieldRA','fieldDec',
                              'filtSkyBrightness', 'airmass', 'dist2Moon', 'moonAlt'], sql)
data = data[::stride]


mysky = np.zeros(data.size, dtype=float)

sm = sb.SkyModel(mags=True)
ind2filt = {'u':0, 'g':1,'r':2,'i':3,'z':4,'y':5}

for i,dat in enumerate(data):
    sm.setRaDecMjd([dat['fieldRA']], [dat['fieldDec']], dat['expMJD'])
    mysky[i] = sm.returnMags()[0][ind2filt[dat['filter']]]

print '%i pointings outside airmass limits of model' % np.size(mysky[ np.isinf(mysky)])

fig,axes = plt.subplots(2,3)
axes = axes.ravel()
keys = ['u','g','r','i','z','y']

print 'filter, median OpSim, median ESO'
for i,key in enumerate(keys):
    good = np.where((data['filter'] == key) & (mysky != 0))
    axes[i].plot(mysky[good],data['filtSkyBrightness'][good]-mysky[good], 'ko', alpha=.1)
    axes[i].set_title(r'$'+key+'$')
    axes[i].set_ylabel(r'$\mu_{OpSim}-\mu_{ESO}$')
    axes[i].set_xlabel(r'$\mu_{ESO}$')
    axes[i].invert_xaxis()
    axes[i].tick_params(axis='x', labelsize=7)
    axes[i].xaxis.set_major_formatter(FormatStrFormatter('%0.1f'))
    xloc = plt.MaxNLocator(7)
    axes[i].xaxis.set_major_locator(xloc)
    print key+', %.2f, %.2f' % (np.median(data['filtSkyBrightness'][good]),
                                np.median(mysky[good]))

fig.tight_layout()
fig.savefig('skycompare.png')


key = 'z'
good = np.where((data['filter'] == key) & (mysky != 0))
fig,ax = plt.subplots()

ack = ax.scatter(data['airmass'][good],mysky[good],  c=np.degrees(data['dist2Moon'][good]), alpha=.1 )
cb = plt.colorbar(ack)
ax.invert_yaxis()
ax.set_title('ESO $z$')
ax.set_ylabel('$\mu_{ESO}$')
ax.set_xlabel('Airmass')
cb.set_label('Degrees to Moon')

fig.savefig('zband_eso.png')

fig,ax = plt.subplots()

ack = ax.scatter(data['airmass'][good],  data['filtSkyBrightness'][good],
                    c=np.degrees(data['dist2Moon'][good]), alpha=.1 )
cb = plt.colorbar(ack)
ax.invert_yaxis()
ax.set_title('OpSim $z$')
ax.set_ylabel('$\mu_{OpSim}$')
ax.set_xlabel('Airmass')
cb.set_label('Degrees to Moon')

fig.savefig('zband_opsim.png')

plt.close('all')

fig,ax = plt.subplots()

ax.hist(np.degrees(data['dist2Moon'][good]), bins=50)
ax.set_xlabel('Distance to Moon')
ax.set_title('$z$')
ax.set_ylabel('# of Visits')
fig.savefig('zd2moon.png')
plt.close(fig)

fig,ax = plt.subplots()
ax.plot(np.degrees(data['moonAlt'][good]),np.degrees(data['dist2Moon'][good]),'ko', alpha=.1)
ax.set_xlabel('Moon Altitude (degrees)')
ax.set_title('$z$')
ax.set_ylabel('Distance to the Moon (degrees)')
fig.savefig('zmoonaltmoond.png')
plt.close(fig)
