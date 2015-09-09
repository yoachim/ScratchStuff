# Make count maps of only some of the chips

import numpy as np
import matplotlib.pyplot as plt
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.db as db
from lsst.sims.maf.plots import PlotHandler
import healpy as hp


database = 'enigma_1189_sqlite.db'
opsdb = db.OpsimDatabase(database)
outDir = 'Camera'
resultsDb = db.ResultsDb(outDir=outDir)


rafts = ['R:0,1', 'R:0,2', 'R:0,3',
         'R:1,0', 'R:1,1', 'R:1,2', 'R:1,3', 'R:1,4',
         'R:2,0', 'R:2,1', 'R:2,2', 'R:2,3', 'R:2,4',
         'R:3,0', 'R:3,1', 'R:3,2', 'R:3,3', 'R:3,4',
         'R:4,1', 'R:4,2', 'R:4,3',
        ]
rafts2 = ['R:0,1',  'R:0,3',
         'R:1,1', 'R:1,3',
         'R:2,0', 'R:2,2', 'R:2,4',
         'R:3,1', 'R:3,3',
         'R:4,1',  'R:4,3',
        ]
rafts1 = ['R:2,2', 'R:2,3', 'R:2,4',
         'R:3,0', 'R:3,1', 'R:3,2', 'R:3,3', 'R:3,4',
         'R:4,0', 'R:4,1', 'R:4,2', 'R:4,3',
        ]

sensors = ['S:0,0', 'S:0,1', 'S:0,2',
           'S:1,0', 'S:1,1', 'S:1,2',
           'S:2,0', 'S:2,1', 'S:2,2',]

chips1 = []
for raft in rafts1:
    for sensor in sensors:
        chips1.append(raft+','+sensor)

chips2 = []
for raft in rafts2:
    for sensor in sensors:
        chips2.append(raft+' '+sensor)


nside = 8

metric = metrics.CountMetric(col='expMJD')
slicer = slicers.HealpixSlicer(nside=nside, useCamera=True, chipNames=chips1)
sql='night < 100 and filter="r"'
bundle = metricBundles.MetricBundle(metric,slicer,sql)
bg = metricBundles.MetricBundleGroup({0:bundle}, opsdb,
                                     outDir=outDir, resultsDb=resultsDb)
bg.runAll()
bg.plotAll()
