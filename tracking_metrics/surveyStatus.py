import numpy as np
import healpy as hp
from lsst.sims.utils import _hpid2RaDec
from scipy.spatial import cKDTree as kdtree
from builtins import zip

"""
Want to make some maps that can track the progress of a survey

How should one go about this?

what are the things that are easy to track?
These are things where we can have the map and easily add a new single observation.
Just select the relevant healpixels and use some algorithm to make it happen

number (count)
    probably want total and per-filter
co-added depth (per filter)
means (can save the mean value and the weight, then update weight when adding new visits)
exposure time (if using variable exposure time)
min,max


Things that are harder to track:
median, rms ().  One possibility would be to pre-declare some histogram bins (maybe 10-20), and then
compute median and rms values from the hist.

How to incorporate this? maybe a dict of survey_status objects with keys u,g,r,i,z,y,all?

Decision flow:
which mode to observe in:
1) Check which modes can observe
2) Compare how much each proposal wants to observe (scale between 0-1)
    * maybe just set wide survey to 0.6
    * DD fields peak when they would be observed crossing meridian
    * Maybe each DD as it's own mode? I think that makes the most sense. Then they could dither around differently
3) Use desire and number of observations that have been taken in each mode to pick one. ()

mode then picks a filter. calc a contour to observe in each filter. Can then pick between them based on
1) The ratio of obs in filter_x / total obs in each filter (maybe take a mean?)
2) what the current filter is
3) the percentile conditions in each contour (make sure we don't observe u when moon up even if behind)
(could multi-processor this over 6 chips I suppose, since it will involve smoothing healpix maps)

Pass the winning contour to the tesselate pointings routine

Once pointings are set, can set a starting rotation angle to try and keep that dithered I guess.

For each proposed visit to the queue, can set a max-gap between time and expected execution time,
if gap gets too large, flush queue and re-fill

"""


class HealpixLookup(object):
    def __init__(self, nside=64, radius=1.75):
        hpindex = np.arange(hp.nside2npix(nside))
        ra, dec = _hpid2RaDec(nside, hpindex)

        self._setRad(radius=radius)
        self._buildTree(ra, dec)

    def lookup(self, ra, dec):
        """
        Lookup the healpixels that overlap a pointing
        """
        x, y, z = self._treexyz(ra, dec)
        indices = self.hptree.query_ball_point((x, y, z), self.rad)
        return indices

    def _treexyz(self, ra, dec):
        """Calculate x/y/z values for ra/dec points, ra/dec in radians."""
        # Note ra/dec can be arrays.
        x = np.cos(dec) * np.cos(ra)
        y = np.cos(dec) * np.sin(ra)
        z = np.sin(dec)
        return x, y, z

    def _buildTree(self, ra, dec, leafsize=100):
        """Build KD tree on simDataRA/Dec and set radius (via setRad) for matching.

        simDataRA, simDataDec = RA and Dec values (in radians).
        leafsize = the number of Ra/Dec pointings in each leaf node."""
        if np.any(np.abs(ra) > np.pi*2.0) or np.any(np.abs(dec) > np.pi*2.0):
            raise ValueError('Expecting RA and Dec values to be in radians.')
        x, y, z = self._treexyz(ra, dec)
        data = list(zip(x, y, z))
        if np.size(data) > 0:
            try:
                self.hptree = kdtree(data, leafsize=leafsize, balanced_tree=False, compact_nodes=False)
            except TypeError:
                self.hptree = kdtree(data, leafsize=leafsize)
        else:
            raise ValueError('RA and Dec should have length greater than 0.')

    def _setRad(self, radius=1.75):
        """Set radius (in degrees) for kdtree search.

        kdtree queries will return pointings within rad."""
        x0, y0, z0 = (1, 0, 0)
        x1, y1, z1 = self._treexyz(np.radians(radius), 0)
        self.rad = np.sqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)


def loadSurveyStatus(filename):
    """
    Check which Class it is and restore it.
    """
    pass


class BaseSurveyStatus(object):
    """
    Keep track of various survey statistics continuously
    """


class SurveyStatusSky(BaseSurveyStatus):
    """
    An object that tracks the progress of a survey
    """

    def __init__(self, nside=64, dtype=float, status_name=None):

        self.nside = nside
        self.survey_map = np.zeros(hp.nside2npix(self.nside), dtype=dtype)
        # Keep a count of how many total visits have been added
        self.visit_counter = 0
        # Which attributes must be the same to count two status objects as equal
        self.eq_check = ['__class__', 'nside', 'visit_counter', 'eq_check']

    def __eq__(self, other_survey):
        """
        Check if two surveyStatus objects are equal
        """
        result = False
        check_attrs = []
        for attr in self.eq_check:
            check_attrs.append(getattr(other_survey, attr) == getattr(self, attr))
        check_attrs.append(np.array_equal(self.survey_map, other_survey.survey_map))
        if False not in check_attrs:
            result = True
        return result

    def __ne__(self, other_survey):
        if self == other_survey:
            return False
        else:
            return True

    def add_visit(self, visit, pixels):

        self.visit_counter += 1
        pass

    def save(self, filename=None):
        """
        Save the current status of the survey
        """

    def return_survey_map(self, **kwargs):
        return self.survey_map


class countFilterStatus(SurveyStatusSky):
    def __init__(self, nside=64, dtype=int, status_name=None,
                 filter_name=['r']):
        super(countFilterStatus, self).__init__(nside=nside, dtype=dtype,
                                                status_name=None)
        self.filter_names = filter_name

    def add_visit(self, visit, pixels):
        if visit.filter in self.filter_names:
            self.survey_map[pixels] += 1
        self.visit_counter += 1


class CoaddM5Status(SurveyStatusSky):
    def __init__(self, nside=64, dtype=float, status_name=None,
                 filter_name='r'):

        super(CoaddM5Status, self).__init__(nside=nside, dtype=dtype,
                                            status_name=None)
        self.filter_name = filter_name
        self.flux_map = np.zeros(hp.nside2npix(self.nside), dtype=dtype)

    def add_visit(self, visit, pixels):
        if visit.filter == self.filter_name:
            self.flux_map[pixels] += 10.**(0.8*visit.fiveSigmaDepth)
            self.survey_map[pixels] = 1.25 * np.log10(self.flux_map[pixels])
        self.visit_counter += 1


class HasTemplateStatus(SurveyStatusSky):
    """
    track if there is a good template for a given helpixel
    """
    def __init__(self, nside=64, dtype=bool, status_name=None,
                 filter_name='r', max_seeing=0.8, min_m5=23.):

        super(HasTemplateStatus, self).__init__(nside=nside, dtype=dtype,
                                                status_name=None)
        self.max_seeing = max_seeing
        self.min_m5 = min_m5
        self.filter_name = filter_name

    def add_visit(self, visit, pixels):
        if visit.filter == self.filter_name:
            if (visit.seeing <= self.max_seeing) & (visit.fiveSigmaDepth >= self.min_m5):
                self.survey_map[pixels] = True
        self.visit_counter += 1


class LastObserved(SurveyStatusSky):
    def __init__(self, nside=64, dtype=float, status_name=None,
                 filter_name='any'):
        super(LastObserved, self).__init__(nside=nside, dtype=dtype,
                                           status_name=None)
        self.filter_name = filter_name
        self.survey_map.fill(hp.UNSEEN)

    def add_visit(self, visit, pixels):
        if (visit.filter == self.filter_name) | (self.filter_name == 'any'):
            self.survey_map[pixels] = np.maximum(self.survey_map[pixels], visit.expMJD)
        self.visit_counter += 1


class NightCount(SurveyStatusSky):
    def __init__(self, nside=64, dtype=int, status_name=None,
                 filter_name='any'):
        super(NightCount, self).__init__(nside=nside, dtype=dtype,
                                         status_name=None)
        self.filter_name = filter_name
        self.current_night = -1

    def add_visit(self, visit, pixels):
        if (visit.filter == self.filter_name) | (self.filter_name == 'any'):
            if visit.night != self.current_night:
                self.survey_map *= 0
                self.current_night = visit.night
            self.survey_map[pixels] += 1

    def return_survey_map(self, night=None):
        # If we are no longer on the correct night, return empty map
        if night is not None:
            if night != self.current_night:
                return self.survey_map*0
        return self.survey_map


class WeightedSky(BaseSurveyStatus):

    def __init__(self, nside=64, dtype=float, status_name=None):
        super(WeightedSky, self).__init__(nside=nside, dtype=dtype,
                                          status_name=None)

        self.weight_map = np.zeros(hp.nside2npix(self.nside), dtype=float)

    def add_visit(self, visit, pixels):
        pass


class HistSky(SurveyStatusSky):
    """
    Keep a histogram of values at each point in the sky that can be used
    to compute more compicated statistics.
    """
    def __init__(self, nside=64, bins=None, dtype=float, status_name=None):
        super(HistSky, self).__init__(nside=nside, dtype=dtype,
                                      status_name=None)

        if bins is None:
            self.bins = np.arange(11)
        else:
            self.bins = np.array(bins)
        self.nbins = self.bins.size
        self.survey_map = np.zeros((hp.nside2npix(self.nside), self.nbins-1), dtype=dtype)

    def add_visit(self, visit, pixels):
        pass


class RotationHistSky(HistSky):
    """
    Keep track of what rotation angles have been used
    """

    def __init__(self, nside=64, bins=None, filter_names=['u', 'g', 'r', 'i', 'z', 'y'],
                 dtype=float, status_name=None):
        if bins is None:
            bins = np.linspace(0., 2.*np.pi, 9.)

        super(RotationHistSky, self).__init__(nside=nside, dtype=dtype,
                                              status_name=None, bins=bins)
        self.filter_names = filter_names

    def add_visit(self, visit, pixels):
        if visit.filter in self.filter_names:
            index = np.searchsorted(self.bins, visit.rotSkyPos, side='right') - 1
            self.survey_map[pixels, index] += 1

