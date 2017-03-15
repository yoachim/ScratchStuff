#!/usr/bin/env python

from __future__ import print_function
from builtins import zip
from builtins import str
from builtins import range
import os
import argparse
import copy
import numpy as np
import warnings
import matplotlib
# Set matplotlib backend (to create plots where DISPLAY is not set).
matplotlib.use('Agg')
import matplotlib.cm as cm

import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.plots as plots
import lsst.sims.maf.utils as utils


def makeBundleList(dbFile, night=1, nside=64, latCol='ditheredDec', lonCol='ditheredRA'):
    """
    Make a bundleList of things to run
    """

    # Construct sql queries for each filter and all filters
    filters = ['u', 'g', 'r', 'i', 'z', 'y']
    sqls = ['filter="%s"' % f for f in filters]
    sqls.append('')

    bundleList = []
    plotFuncs = [plots.LambertSkyMap()]

    reg_slicer = slicers.HealpixSlicer(nside=nside, lonCol=lonCol, latCol=latCol)
    altaz_slicer = slicers.HealpixSlicer(nside=nside, latCol='zenithDistance',
                                         lonCol='azimuth', useCache=False)
    for sql in sqls:

        # Number of exposures
        metric = metrics.CountMetric('expMJD', metricName='Nvisits')
        bundle = metricBundles.MetricBundle(metric, reg_slicer, sql)
        bundleList.append(bundle)
        bundle = metricBundles.MetricBundle(metric, altaz_slicer, sql, plotFuncs=plotFuncs)
        bundleList.append(bundle)

    return metricBundles.makeBundlesDictFromList(bundleList)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Python script to run the science performance metrics.')
    parser.add_argument('dbFile', type=str, default=None, help="full file path to the opsim sqlite file")
    parser.add_argument("--outDir", type=str, default='./Out', help='Output directory for MAF outputs.' +
                        ' Default "Out"')
    parser.add_argument("--nside", type=int, default=64,
                        help="Resolution to run Healpix grid at (must be 2^x). Default 64.")
    parser.add_argument("--lonCol", type=str, default='fieldRA',
                        help="Column to use for RA values (can be a stacker dither column)." +
                        " Default=fieldRA.")
    parser.add_argument("--latCol", type=str, default='fieldDec',
                        help="Column to use for Dec values (can be a stacker dither column)." +
                        " Default=fieldDec.")
    parser.add_argument('--night', type=int, default=1)

    parser.set_defaults()
    args, extras = parser.parse_known_args()

    bundleDict = makeBundleList(args.dbFile, nside=args.nside, lonCol=args.lonCol, latCol=args.latCol,
                                night=args.night)

    # Set up / connect to resultsDb.
    resultsDb = db.ResultsDb(outDir=args.outDir)
    # Connect to opsimdb.
    opsdb = utils.connectOpsimDb(args.dbFile)


    # Set up metricBundleGroup.
    group = metricBundles.MetricBundleGroup(bundleDict, opsdb, uotDir=args.outDir, resultsDb=resultsDb)
    group.runAll()
    group.plotAll()

