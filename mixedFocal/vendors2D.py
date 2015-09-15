import numpy as np
import matplotlib.pyplot as plt
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.metricBundles as metricBundles
import lsst.sims.maf.db as db
from lsst.sims.maf.plots import PlotHandler
from lsst.sims.maf.stackers import BaseStacker
import healpy as hp


class RotPairStacker(BaseStacker):
    """
    Add a new column that rotates the second (third, 4rth, etc) visit in a night by
    the given rotAmount (radians).
    """
    def __init__(self, rotAmount=np.pi , fieldIDCol='fieldID', nightCol='night', rotCol='rotSkyPos'):

        self.units = ['radians']
        self.colsAdded = ['rotatedRotSkyPos']
        self.nightCol = nightCol
        self.rotCol = rotCol
        self.fieldIDCol = fieldIDCol
        self.rotAmount=np.pi

        self.colsReq = [nightCol, rotCol, fieldIDCol]

    def run(self, simData):
        simData=self._addStackers(simData)

        # Fill in the old rotation angles as default
        simData['rotatedRotSkyPos'] = simData['rotatedRotSkyPos']*0+simData[self.rotCol]

        simData.sort(order=[self.nightCol,self.fieldIDCol])

        unights = np.unique(simData[self.nightCol])
        left = np.searchsorted(simData[self.nightCol], unights, side='left')
        right = np.searchsorted(simData[self.nightCol], unights, side='right')

        for l,r in zip(left,right):
            if r-l > 1:
                ufid = np.unique(simData[self.fieldIDCol][l:r])
                fLeft = np.searchsorted(simData[self.fieldIDCol][l:r], ufid, side='left')
                fRight = np.searchsorted(simData[self.fieldIDCol][l:r],ufid, side='right')
                for fL,fR in zip(fLeft,fRight):
                    if fR-fL > 1:
                        angleStart = simData[self.rotCol][l:r][fL:fR][0]
                        newAngles = self.rotAmount*np.arange(fR-fL)+angleStart
                        newAngles = newAngles % (2.*np.pi)
                        simData['rotatedRotSkyPos'][l:r][fL:fR] = newAngles
        return simData




database = 'enigma_1189_sqlite.db'
opsdb = db.OpsimDatabase(database)
outDir = '2DCamera'
resultsDb = db.ResultsDb(outDir=outDir)


rafts = [         'R:0,1', 'R:0,2', 'R:0,3',
         'R:1,0', 'R:1,1', 'R:1,2', 'R:1,3', 'R:1,4',
         'R:2,0', 'R:2,1', 'R:2,2', 'R:2,3', 'R:2,4',
         'R:3,0', 'R:3,1', 'R:3,2', 'R:3,3', 'R:3,4',
                  'R:4,1', 'R:4,2', 'R:4,3'
        ]

rafts2 = [         'R:0,1', 'R:0,2', 'R:0,3',
         'R:1,0', 'R:1,1', 'R:1,2', 'R:1,3', 'R:1,4',
         'R:2,0', 'R:2,1'
        ]

rafts1 = [                 'R:2,2', 'R:2,3', 'R:2,4',
         'R:3,0', 'R:3,1', 'R:3,2', 'R:3,3', 'R:3,4',
                  'R:4,1', 'R:4,2', 'R:4,3'
        ]

sensors = ['S:0,0', 'S:0,1', 'S:0,2',
           'S:1,0', 'S:1,1', 'S:1,2',
           'S:2,0', 'S:2,1', 'S:2,2',]

chips1 = []
for raft in rafts1:
    for sensor in sensors:
        chips1.append(raft+' '+sensor)

chips2 = []
for raft in rafts2:
    for sensor in sensors:
        chips2.append(raft+' '+sensor)



filters = ['u','g','r','i','z','y']
year =  10

read = True

nside = 16

extras = ['None', 'dither', 'rotation']

for extra in extras:
    md = ''
    if extra == 'dither':
        latCol = 'ditheredDec'
        lonCol = 'ditheredRA'
        md += ' Dithered,'
    else:
        latCol = 'fieldDec'
        lonCol = 'fieldRA'

    if extra == 'rotation':
        rotSkyPosColName = 'rotatedRotSkyPos'
        md += ' Rotated 2nd visit,'
    else:
        rotSkyPosColName = 'rotSkyPos'


    for filterName in filters:
        bundleList=[]
        mdf = md+' %s,' % filterName

        sql = 'filter="%s" and night < %i' % (filterName,year*365.25)
        metric = metrics.AccumulateCountMetric()

        bins = np.arange(0,np.ceil(year*365.25)+1,1)
        slicer = slicers.Healpix2dSlicer(nside=nside, bins=bins, useCamera=True,
                                         latCol=latCol, lonCol=lonCol,
                                         rotSkyPosColName=rotSkyPosColName)
        bundle = metricBundles.MetricBundle(metric,slicer,sql, metadata=mdf+'Single Vendor')
        bundle.Single=True
        bundleList.append(bundle)

        slicer = slicers.Healpix2dSlicer(nside=nside, bins=bins, useCamera=True, chipNames=chips1,
                                         latCol=latCol, lonCol=lonCol,
                                         rotSkyPosColName=rotSkyPosColName)
        bundle = metricBundles.MetricBundle(metric,slicer,sql, metadata=mdf+'Vendor 1')
        bundleList.append(bundle)

        slicer = slicers.Healpix2dSlicer(nside=nside, bins=bins, useCamera=True, chipNames=chips2,
                                         latCol=latCol, lonCol=lonCol,
                                         rotSkyPosColName=rotSkyPosColName)
        bundle = metricBundles.MetricBundle(metric,slicer,sql, metadata=mdf+'Vendor 2')
        bundleList.append(bundle)


        bd = metricBundles.makeBundlesDictFromList(bundleList)
        bg = metricBundles.MetricBundleGroup(bd, opsdb,
                                             outDir=outDir, resultsDb=resultsDb)
        if read:
            bg.readAll()
        else:
            bg.runAll()
        bg.plotAll(closefigs=True)

        nLimits = [2,4,8,16,32]
        pix2area = hp.nside2pixarea(nside, degrees=True)
        for bundle in bundleList:
            if hasattr(bundle,'Single'):
                refBundle = bundle
        for bundle in bundleList:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            for limit in nLimits:
                good = np.where(bundle.metricValues > limit)
                goodRef = np.where(refBundle.metricValues > limit)
                nRef = np.zeros(bundle.metricValues.shape)
                nRef[goodRef] = 1.
                nRef = np.sum(nRef, axis=0)
                nHp = np.zeros(bundle.metricValues.shape)
                nHp[good] = 1.
                nHp = np.sum(nHp, axis=0)
                ax.plot(bins,nHp/nRef, label='%i' % limit)
            ax.set_xlabel('Night')
            ax.set_ylabel('Area/Single Vendor Area')
            #ax.set_ylim([0,35000])
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles, labels, loc='lower right')
            ax.set_title(bundle.metadata)
            filename = outDir+'/%s' % 'timeEvo'+'_'+bundle.metadata.replace(' ','').replace(',','_')+'.png'
            fig.savefig(filename)
            print 'Made file %s' % filename
            plt.close(fig)
