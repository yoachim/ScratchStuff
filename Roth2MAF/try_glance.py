
import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.stackers as stackers
import lsst.sims.maf.plots as plots
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.utils as utils

from lsst.sims.maf.bundles import glanceBundle

conn = db.SimpleDatabase('observations.sqlite')
outDir='Glance'
colmap = {'ra': 'ra', 'dec': 'dec', 'mjd': 'mjd',
              'exptime': 'exptime', 'visittime': 'exptime', 'alt': 'alt',
              'az': 'az', 'filter': 'filter', 'fiveSigmaDepth': 'fivesigmadepth',
              'night': 'night', 'slewtime': 'slewtime', 'seeingGeom': 'seeing'}

gb = glanceBundle(colmap_dict=colmap)
outDir='Glance'
resultsDb = db.ResultsDb(outDir=outDir)

group = metricBundles.MetricBundleGroup(gb, conn, outDir=outDir, resultsDb=resultsDb, 
                                        runName='roth_1')

group.runAll()
group.plotAll()
