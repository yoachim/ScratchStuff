from glanceBundle import glanceBundle
import lsst.sims.maf.utils as utils
import lsst.sims.maf.db as db
import lsst.sims.maf.metricBundles as metricBundles

# Let's run the Glance summary on several different surveys!

outDir = 'pontus_1168'
resultsDb = db.ResultsDb(outDir=outDir)
conn = utils.connectOpsimDb('pontus_1168.db')
gb = glanceBundle()
group = metricBundles.MetricBundleGroup(gb, conn, outDir=outDir, resultsDb=resultsDb)
group.runAll()
group.plotAll()
