# Make a text file dump for Z

import numpy as np
import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metricBundles as metricBundles
import ephem

from lsst.sims.maf.utils.telescopeInfo import TelescopeInfo
from lsst.sims.maf.stackers import mjd2djd


telescope = TelescopeInfo('LSST')
observatory = ephem.Observer()
observatory.lat = telescope.lat
observatory.lon = telescope.lon
observatory.elevation = telescope.elev



# Set up the database connection
opsdb = db.OpsimDatabase('sqlite:///enigma_1189_sqlite.db')
outDir = 'Dump'
resultsDb = db.ResultsDb(outDir=outDir)


metricList =[]
metricList.append(metrics.PassMetric('expMJD'))
metricList.append(metrics.CountMetric('filter'))
metricList.append(metrics.CountMetric('solarElong'))
metricList.append(metrics.CountMetric('fieldRA'))
metricList.append(metrics.CountMetric('fieldDec'))
metricList.append(metrics.CountMetric('filtSkyBrightness'))
metricList.append(metrics.CountMetric('finSeeing'))
metricList.append(metrics.CountMetric('fiveSigmaDepth'))
metricList.append(metrics.CountMetric('NEOGeoDist'))
metricList.append(metrics.CountMetric('eclipLat'))
metricList.append(metrics.CountMetric('eclipLon'))
metricList.append(metrics.CountMetric('altitude'))
metricList.append(metrics.CountMetric('azimuth'))

slicer = slicers.UniSlicer()
bDict={}

for i,metric in enumerate(metricList):
    bDict[i] = metricBundles.MetricBundle(metric, slicer,'')
bgroup = metricBundles.MetricBundleGroup(bDict, opsdb, outDir=outDir, resultsDb=resultsDb)
bgroup.runAll()

data = np.load('Dump/opsim_Pass_expMJD_UNIS.npz')
ack = data['metricValues'][0].copy()

f = open('Dump/enigma_dump.dat','w')
print >>f, 'obsHistID, MJD, fieldRA, fieldDec, eclipLat, eclipLon, filter, filtskybrightness, finSeeing, solarElong, fiveSigmaDepth,  NEOGeoDist, NEOHelioX, NEOHelioY, sunAlt, sunAz'

sun = ephem.Sun()

for blah in ack:
    observatory.date = mjd2djd( blah['expMJD'])
    sun.compute(observatory)
    sunAlt = sun.alt
    sunAz = sun.az
    print >>f, blah['obsHistID'], blah['expMJD'], np.degrees(blah['fieldRA']), np.degrees(blah['fieldDec']), np.degrees(blah['altitude']),np.degrees(blah['azimuth']),np.degrees(blah['eclipLat']), np.degrees(blah['eclipLon']),blah['filter'], blah['filtSkyBrightness'], blah['finSeeing'], blah['solarElong'], blah['fiveSigmaDepth'],  blah['NEOGeoDist'], blah['NEOHelioX'], blah['NEOHelioY'], np.degrees(sunAlt), np.degrees(sunAz)

f.close()
