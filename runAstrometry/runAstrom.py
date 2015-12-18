import matplotlib.pyplot as plt
import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.plots as plots
import healpy as hp

# Run the Astrometry metics on a number of runs

runNames = ['enigma_1189', 'ops2_1093']
nside = 64

slicer = slicers.HealpixSlicer(nside=nside)

# Make it so we don't bother with the silly power spectra
plotFuncs = [plots.HealpixSkyMap(), plots.HealpixHistogram()]

# The old runs have the seeing in finSeeing
seeingCol = 'finSeeing'

# Try it out for a 20th mag star with a flat SED (can change mag or to OBAFGKM)
rmag = 20.
SedTemplate='flat'

# Use all the observations. Can change if you want a different time span
sqlconstraint = ''

# run some summary stats on everything
summaryMetrics = [metrics.MedianMetric()]

for run in runNames:
    # Open the OpSim database
    opsdb = db.OpsimDatabase(run+'_sqlite.db')
    # Set where the output should go
    outDir = run
    resultsDb = db.ResultsDb(outDir=outDir)

    bundleList = []

    # Configure the metrics
    metric = metrics.ParallaxMetric(rmag=rmag, seeingCol=seeingCol, SedTemplate=SedTemplate)
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, runName=run,
                                        plotFuncs=plotFuncs, summaryMetrics=summaryMetrics)
    bundleList.append(bundle)

    metric=metrics.ProperMotionMetric(rmag=rmag, seeingCol=seeingCol, SedTemplate=SedTemplate)
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, runName=run,
                                        plotFuncs=plotFuncs, summaryMetrics=summaryMetrics)
    bundleList.append(bundle)

    metric = metrics.ParallaxCoverageMetric(rmag=rmag, seeingCol=seeingCol, SedTemplate=SedTemplate)
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, runName=run,
                                        plotFuncs=plotFuncs, summaryMetrics=summaryMetrics)
    bundleList.append(bundle)

    metric = metrics.ParallaxHADegenMetric(rmag=rmag, seeingCol=seeingCol, SedTemplate=SedTemplate)
    bundle = metricBundles.MetricBundle(metric, slicer, sqlconstraint, runName=run,
                                        plotFuncs=plotFuncs, summaryMetrics=summaryMetrics)
    bundleList.append(bundle)

    # Run everything and make plots
    bundleDict = metricBundles.makeBundlesDictFromList(bundleList)
    bgroup = metricBundles.MetricBundleGroup(bundleDict, opsdb, outDir=outDir, resultsDb=resultsDb)
    bgroup.runAll()
    bgroup.plotAll()
