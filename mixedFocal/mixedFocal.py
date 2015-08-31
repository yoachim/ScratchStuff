# Do some experiments to check how well we perform with a mixed vendor chips
import numpy as np
import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.stackers as stackers
import lsst.sims.maf.plots as plots
import lsst.sims.maf.metricBundles as metricBundles


# Make a new stacker to modify single visit depths
class V2m5Stacker(stackers.BaseStacker):
    """
    Make a stacker to create a new coaddM5 column for a
    different sensitivity chip.
    """
    def __init__(self, filterCol='filter', m5Col='fiveSigmaDepth',
                 m5Deltas={'u':-0.3,'g':0.05,'r':0.03}):
        """
        m5Deltas = key for each filter, value is added to the single visit
        m5 values.
        """
        self.colsReq = [filterCol,m5Col]
        self.colsAdded = ['v2fiveSigmaDepth']
        self.filterCol = filterCol
        self.m5Col = m5Col
        self.m5Deltas = m5Deltas

    def run(self, simData):
        simData=self._addStackers(simData)

        for filterName in self.m5Deltas.keys():
            good = np.where(simData[self.filterCol] == filterName)[0]
            simData['v2fiveSigmaDepth'][good] = simData[self.m5Col][good] + self.m5Deltas[filterName]
        return simData


class MixedM5Metric(mertics.BaseMetric):
    """
    Modify the m5 values on some of the rafts then compute the final co-added m5.
    """

    def __init__(self, m5v1Col = 'fiveSigmaDepth', m5v2Col = 'v2fiveSigmaDepth', units='mag',
                 rafts1 = [], rafts2= []):
        super(MixedM5Metric, self).__init__(col=[m5v1Col,m5v2Col],units=units, **kwargs)
        self.m5v1Col = m5v1Col
        self.m5v2Col = m5v2Col
