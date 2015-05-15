# Make a text file dump for Z

import numpy as np
import matplotlib.pyplot as plt
import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.stackers as stackers
import lsst.sims.maf.metricBundles as metricBundles
from lsst.sims.maf.plots import NeoDetectPlotter
import lsst.sims.maf.plots as plotters

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

slicer = slicers.UniSlicer()
bDict={}

for i,metric in enumerate(metricList):
    bDict[i] = metricBundles.MetricBundle(metric, slicer,'')
bgroup = metricBundles.MetricBundleGroup(bDict, opsdb, outDir=outDir, resultsDb=resultsDb)
bgroup.runAll()

data = np.load('Dump/opsim_Pass_expMJD_UNIS.npz')
ack = data['metricValues'][0].copy()

f = open('Dump/enigma_dump.dat','w')
print >>f, 'obsHistID, fieldRA, fieldDec, eclipLat, eclipLon, filter, filtskybrightness, finSeeing, solarElong, fiveSigmaDepth,  NEOGeoDist, NEOHelioX, NEOHelioY'

for blah in ack:
    print >>f, blah['obsHistID'], blah['fieldRA'], blah['fieldDec'], blah['eclipLat'], blah['eclipLon'],blah['filter'], blah['filtSkyBrightness'], blah['finSeeing'], blah['solarElong'], blah['fiveSigmaDepth'],  blah['NEOGeoDist'], blah['NEOHelioX'], blah['NEOHelioY']

f.close()
