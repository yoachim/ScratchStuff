from lsst.sims.maf.bundles import glanceBundle
import lsst.sims.maf.utils as utils
import lsst.sims.maf.db as db
import lsst.sims.maf.metricBundles as metricBundles

# Let's run the Glance summary on several different surveys!

outDir = 'pontus_2026'
resultsDb = db.ResultsDb(outDir=outDir)
conn = utils.connectOpsimDb('pontus_2026.db')

gb = glanceBundle(utils.opsimColMapDict())
group = metricBundles.MetricBundleGroup(gb, conn, outDir=outDir, resultsDb=resultsDb, dbTable='SummaryAllProps')
group.runAll()
group.plotAll()
