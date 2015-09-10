# Make count maps of only some of the chips

import numpy as np
import matplotlib.pyplot as plt
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.db as db
from lsst.sims.maf.plots import PlotHandler
import healpy as hp
from lsst.sims.maf.stackers import BaseStacker


class RotPairStacker(BaseStacker):
    """
    Add a new column that rotates the second (third, 4rth, etc) visit in a night by
    the given rotAmount (radians).
    """
    def __init__(self, rotAmount=np.pi , fieldIDCol='fieldID', nightCol='night', rotCol='rotSkyPos'):

        self.units = ['radians']
        self.colsAdded = ['rotatedRotSkyPos']
        self.nightCol = nightCol
        self.rotCol = rotCol
        self.fieldIDCol = fieldIDCol
        self.rotAmount=np.pi

        self.colsReq = [nightCol, rotCol, fieldIDCol]

    def run(self, simData):
        simData=self._addStackers(simData)

        # Fill in the old rotation angles as default
        simData['rotatedRotSkyPos'] = simData['rotatedRotSkyPos']*0+simData[self.rotCol]

        simData.sort(order=[self.nightCol,self.fieldIDCol])

        unights = np.unique(simData[self.nightCol])
        left = np.searchsorted(simData[self.nightCol], unights, side='left')
        right = np.searchsorted(simData[self.nightCol], unights, side='right')

        for l,r in zip(left,right):
            if r-l > 1:
                ufid = np.unique(simData[self.fieldIDCol][l:r])
                fLeft = np.searchsorted(simData[self.fieldIDCol][l:r], ufid, side='left')
                fRight = np.searchsorted(simData[self.fieldIDCol][l:r],ufid, side='right')
                for fL,fR in zip(fLeft,fRight):
                    if fR-fL > 1:
                        angleStart = simData[self.rotCol][l:r][fL:fR][0]
                        newAngles = self.rotAmount*np.arange(fR-fL)+angleStart
                        newAngles = newAngles % (2.*np.pi)
                        simData['rotatedRotSkyPos'][l:r][fL:fR] = newAngles
        return simData


database = 'enigma_1189_sqlite.db'
opsdb = db.OpsimDatabase(database)
outDir = 'Camera'
resultsDb = db.ResultsDb(outDir=outDir)

rafts = [         'R:0,1', 'R:0,2', 'R:0,3',
         'R:1,0', 'R:1,1', 'R:1,2', 'R:1,3', 'R:1,4',
         'R:2,0', 'R:2,1', 'R:2,2', 'R:2,3', 'R:2,4',
         'R:3,0', 'R:3,1', 'R:3,2', 'R:3,3', 'R:3,4',
                  'R:4,1', 'R:4,2', 'R:4,3'
        ]

rafts2 = [         'R:0,1', 'R:0,2', 'R:0,3',
         'R:1,0', 'R:1,1', 'R:1,2', 'R:1,3', 'R:1,4',
         'R:2,0', 'R:2,1'
        ]

rafts1 = [                 'R:2,2', 'R:2,3', 'R:2,4',
         'R:3,0', 'R:3,1', 'R:3,2', 'R:3,3', 'R:3,4',
                  'R:4,1', 'R:4,2', 'R:4,3'
        ]

sensors = ['S:0,0', 'S:0,1', 'S:0,2',
           'S:1,0', 'S:1,1', 'S:1,2',
           'S:2,0', 'S:2,1', 'S:2,2',]

chips1 = []
for raft in rafts1:
    for sensor in sensors:
        chips1.append(raft+' '+sensor)

chips2 = []
for raft in rafts2:
    for sensor in sensors:
        chips2.append(raft+' '+sensor)


nside = 64

years = [1]#[1,2,5,10] #[0.1]
bundleList = []
for year in years:
    sql='night < %i and filter="r"' % (356.25*year)

    metric = metrics.CountMetric(col='expMJD', metricName='Vendor1')
    slicer = slicers.HealpixSlicer(nside=nside, useCamera=True, chipNames=chips1,
                                   lonCol='ditheredRA', latCol='ditheredDec')
    bundle = metricBundles.MetricBundle(metric,slicer,sql)
    bundleList.append(bundle)

    metric = metrics.CountMetric(col='expMJD', metricName='Vendor2')
    slicer = slicers.HealpixSlicer(nside=nside, useCamera=True, chipNames=chips2,
                                   lonCol='ditheredRA', latCol='ditheredDec')
    bundle = metricBundles.MetricBundle(metric,slicer,sql)
    bundleList.append(bundle)

    # Now let's try it with the new stacker
    metric = metrics.CountMetric(col='expMJD', metricName='Vendor1, rotAdded')
    slicer = slicers.HealpixSlicer(nside=nside, rotSkyPosColName='rotatedRotSkyPos',
                                   useCamera=True, chipNames=chips1,
                                   lonCol='ditheredRA', latCol='ditheredDec')
    bundle = metricBundles.MetricBundle(metric,slicer,sql)
    bundleList.append(bundle)

    metric = metrics.CountMetric(col='expMJD', metricName='Vendor2, rotAdded')
    slicer = slicers.HealpixSlicer(nside=nside, rotSkyPosColName='rotatedRotSkyPos',
                                   useCamera=True, chipNames=chips2)
    bundle = metricBundles.MetricBundle(metric,slicer,sql)
    bundleList.append(bundle)

bd = metricBundles.makeBundlesDictFromList(bundleList)
bg = metricBundles.MetricBundleGroup(bd, opsdb,
                                     outDir=outDir, resultsDb=resultsDb)



bg.runAll()
bg.plotAll()
