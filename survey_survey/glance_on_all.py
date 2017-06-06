import lsst.sims.maf.utils as utils
import lsst.sims.maf.db as db
import lsst.sims.maf.metricBundles as metricBundles
from lsst.sims.maf.bundles import glanceBundle

# Let's run the Glance summary on several different surveys!


def RothColMapDict():
    """Return a dict that maps column names.
    """
    result = {'ra': 'ra', 'dec': 'dec', 'mjd': 'mjd',
              'exptime': 'exptime', 'visittime': 'exptime', 'alt': 'alt',
              'az': 'az', 'filter': 'filter', 'fiveSigmaDepth': 'fivesigmadepth',
              'night': 'night', 'slewtime': 'slewtime', 'seeingGeom': 'seeing'}

    return result


outDir='Roth'

resultsDb = db.ResultsDb(outDir=outDir)
conn = db.SimpleDatabase('Roth.sqlite')
gb = glanceBundle(colmap_dict=RothColMapDict())
group = metricBundles.MetricBundleGroup(gb, conn, outDir=outDir, resultsDb=resultsDb)
group.runAll()
group.plotAll()



#opsdb = utils.connectOpsimDb('pontus_1074.db')