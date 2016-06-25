import numpy as np
import lsst.sims.skybrightness as sb
import lsst.sims.maf.db as db
import sys

# Let's generate a new column of sky brightnesses and new m5 values

def calcM5(skybright, filters, airmass, seeing, expTime=30.):
    kAtm = {'u':0.451,'g':0.163,'r':0.087,'i':0.065,'z':0.043,'y':0.138}
    Cm = {'u':22.92,'g':24.29,'r':24.33,'i':24.20,'z':24.07,'y':23.69}
    dCm_infinity = {'u':0.67,'g':0.21,'r':0.11,'i':0.08,'z':0.05,'y':0.04}

    m5 = np.zeros(skybright.size, dtype=float)
    ufilters = np.unique(filters)
    for filterName in ufilters:
        good = np.where(filters == filterName)
        m5[good] = Cm[filterName] + 0.5*(skybright[good] - 21.) + 2.5*np.log10(0.7/seeing[good]) + \
            1.25*np.log10(expTime/30.) - kAtm[filterName]*(airmass[good]-1.)

        return m5

if __name__ == '__main__':
    sm = sb.SkyModel(observatory='LSST', mags=True)

    runName = 'minion_1016'
    opsdb = db.OpsimDatabase(runName+'_sqlite.db')
    data = opsdb.fetchMetricData(['fieldRA', 'fieldDec', 'expMJD', 'airmass', 
                                 'filter', 'filtSkyBrightness'], '', distinctExpMJD=False)

    names = ['newsky']
    types = [float]
    newSky = np.zeros(data.size, dtype=zip(names, types))

    maxI = data.size
    for i, d in enumerate(data):
        sm.setRaDecMjd([d['fieldRA']], [d['fieldDec']], d['expMJD'], degrees=False)
        newSky[i] = sm.returnMags()[d['filter']][0]

        progress = i/float(maxI)*100
        text = "\rprogress = %.1f%%"%progress
        sys.stdout.write(text)
        sys.stdout.flush()

    newm5 = calcM5(newSky, data['filter'], data['airmass'], data['FWHMeff'])
    # np.savez('newData.npz', data=data, newSky=newSky)
