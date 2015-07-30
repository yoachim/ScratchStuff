#! /usr/bin/env python
import os, sys, argparse
import numpy as np
# Set matplotlib backend (to create plots where DISPLAY is not set).
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
import healpy as hp
import warnings

import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.stackers as stackers
import lsst.sims.maf.plots as plots
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.utils as utils
import lsst.sims.maf.plots as plotters

def makeBundleList(dbFile, runName=None, nside=128, benchmark='design',
                   lonCol='fieldRA', latCol='fieldDec'):
    """
    make a list of metricBundle objects to look at the scientific performance
    of an opsim run.
    """

    # List to hold everything we're going to make
    bundleList = []
    filtorder = {'u':1,'g':2,'r':3,'i':4,'z':5,'y':6}
    plotList = [plotters.HealpixSkyMap(), plotters.HealpixHistogram()]

    slicer=slicers.HealpixSlicer(nside=nside)
    plotDict = {'cmap':'jet'}

    metricList = []
    metricList.append(metrics.IntraNightGapsMetric())
    metricList.append(metrics.InterNightGapsMetric())
    metricList.append(metrics.IntraNightGapsMetric(reduceFunc=np.size,
                                                   metricName='Number of in-night pairs'))

    filters = ['u','g','r','i','z','y']

    summaryStats = [metrics.MedianMetric()]

    ddList = [ {'group':'Median gap, in-night'}, {'group':'Median gap, between nights'},
               {'group': 'In-night Pairs'} ]

    for f in filters:
        sql = "filter='%s'"%f
        for metric,dd in zip(metricList,ddList):
            displayDict = dd
            displayDict['order'] = filtorder[f]
            bundleList.append(metricBundles.MetricBundle(metric,slicer,sql, plotFuncs=plotList,
                                                         plotDict=plotDict, summaryMetrics=summaryStats,
                                                         displayDict=displayDict))
    for metric,dd in zip(metricList,ddList):
        sql = ''
        displayDict = dd
        displayDict['order'] = 7
        bundleList.append(metricBundles.MetricBundle(metric,slicer,sql, plotFuncs=plotList,
                                                     plotDict=plotDict, summaryMetrics=summaryStats,
                                                     displayDict=displayDict))



    slicer = slicers.OneDSlicer(sliceColName='night', binsize=1)
    metric = metrics.UniqueRatioMetric(col='fieldID')
    sql='filter="g" or filter="r" or filter="i" or filter="z"'
    displayDict = {'group':'Fields Per Night'}
    bundleList.append(metricBundles.MetricBundle(metric,slicer,sql,
                                                 summaryMetrics=summaryStats,
                                                 displayDict=displayDict))
    metric = metrics.CountUniqueMetric(col='fieldID')
    bundleList.append(metricBundles.MetricBundle(metric,slicer,sql,
                                                 summaryMetrics=summaryStats,
                                                 displayDict=displayDict))
    return metricBundles.makeBundlesDictFromList(bundleList)


if __name__=="__main__":

    parser = argparse.ArgumentParser(description='Python script to run MAF with the science performance metrics')
    parser.add_argument('dbFile', type=str, default=None,help="full file path to the opsim sqlite file")

    parser.add_argument("--outDir",type=str, default='./Out', help='Output directory for MAF outputs. Default "Out"')
    parser.add_argument("--nside", type=int, default=128,
                        help="Resolution to run Healpix grid at (must be 2^x). Default 128.")
    parser.add_argument("--lonCol", type=str, default='fieldRA',
                        help="Column to use for RA values (can be a stacker dither column). Default=fieldRA.")
    parser.add_argument("--latCol", type=str, default='fieldDec',
                        help="Column to use for Dec values (can be a stacker dither column). Default=fieldDec.")
    parser.add_argument('--benchmark', type=str, default='design',
                        help="Can be 'design' or 'requested'")
    parser.add_argument('--plotOnly', dest='plotOnly', action='store_true',
                        default=False, help="Reload the metric values from disk and re-plot them.")

    parser.set_defaults()
    args, extras = parser.parse_known_args()

    # Build metric bundles.
    bundleDict = makeBundleList(args.dbFile, nside=args.nside,
                                lonCol=args.lonCol, latCol=args.latCol,
                                benchmark=args.benchmark)

    # Set up / connect to resultsDb.
    resultsDb = db.ResultsDb(outDir=args.outDir)
    # Connect to opsimdb.
    opsdb = utils.connectOpsimDb(args.dbFile)

    # Set up metricBundleGroup.
    group = metricBundles.MetricBundleGroup(bundleDict, opsdb,
                                            outDir=args.outDir, resultsDb=resultsDb)
    # Read or run to get metric values.
    if args.plotOnly:
        group.readAll()
    else:
        group.runAll()
    # Make plots.
    group.plotAll()

    # Get config info and write to disk.
    utils.writeConfigs(opsdb, args.outDir)
