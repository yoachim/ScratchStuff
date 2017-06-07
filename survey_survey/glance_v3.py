from lsst.sims.maf.bundles import glanceBundle
import lsst.sims.maf.utils as utils
import lsst.sims.maf.db as db
import lsst.sims.maf.metricBundles as metricBundles

# Let's run the Glance summary on several different surveys!

def v3Colmap():
    """Return a dict that maps column names.
    """
    result = {'ra': 'fieldRA', 'dec': 'fieldDec', 'mjd': 'expMJD',
              'exptime': 'visitExpTime', 'visittime': 'visitTime', 'alt': 'altitude',
              'az': 'azimuth', 'filter': 'filter', 'fiveSigmaDepth': 'fiveSigmaDepth',
              'night': 'night', 'slewtime': 'slewTime', 'seeingGeom': 'FWHMgeom'}

    return result

outDir = 'minion_1016'
resultsDb = db.ResultsDb(outDir=outDir)
conn = utils.connectOpsimDb('minion_1016_sqlite.db')

gb = glanceBundle(v3Colmap())
group = metricBundles.MetricBundleGroup(gb, conn, outDir=outDir, resultsDb=resultsDb, dbTable='Summary')
group.runAll()
group.plotAll()
