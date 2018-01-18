import lsst.sims.maf.batches as batches
import lsst.sims.maf.db as db
from lsst.sims.maf.metricBundles import MetricBundleGroup

colmap = batches.ColMapDict('barebones')
colmap['ra'] = 'RA'
colmap['seeingEff'] = 'FWHMeff'
colmap['seeingGeom'] = 'FWHM_geometric'
bd = batches.glanceBatch(colmap=colmap)


def run_glance(outDir, dbname):


    conn = db.Database(dbname, defaultTable='observations')
    resultsDb = db.ResultsDb(outDir=outDir)
    mbg = MetricBundleGroup(bd, conn, outDir=outDir, resultsDb=resultsDb)
    mbg.runAll()
    mbg.plotAll()
    conn.close()  


#run_glance('temp_2', 'feature_baseline_1yrs.db')
#run_glance('roll_check_third', 'feature_rolling_third_0yrs.db')
#run_glance('roll_check_half', 'feature_rolling_half_0yrs.db')
#run_glance('baseline_10', 'feature_baseline_10yrs.db')
#run_glance('half_roll_10', 'feature_rolling_half_10yrs.db')
#run_glance('third_roll_10', 'feature_rolling_third_10yrs.db')
#run_glance('half_roll_mask_10', 'feature_rolling_half_mask_10yrs.db')
#run_glance('baseline_0', 'feature_baseline_0yrs.db')
#run_glance('baseline_1', 'feature_baseline_1yrs.db')
run_glance('roll_twothirds', 'feature_rolling_twoThird_10yrs.db')