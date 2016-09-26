import numpy as np
import healpy as hp
import warnings

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

    def __init__(self, nside=64, dtype=float, status_name=None, footprint='D3.5'):

        self.nside = nside
        self.footprint = footprint
        self.survey_map = np.zeros(hp.nside2npix(self.nside), dtype=dtype)
        # Keep a count of how many total visits have been added
        self.visit_counter = 0
        pass

    def __eq__(self, other_survey):
        """
        Check if two surveyStatus objects are equal
        """
    def __ne__(self, other_survey):
        pass

    def observed_healpixels(self, visit):
        # use self.footprint to decide how to calculate which healpix hit silicon in a visit.
        pass

    def add_visit(self, visit, pixels=None):
        if pixels is None:
            pixels = self.find_pixels(visit)

        self.visit_counter += 1
        pass

    def remove_visit(self, visit, pixels=None):
        # Hmm, this makes sense for things like count, but min/max will be hard since we would need a cache of previous values.
        if pixels is None:
            pixels = self.observed_healpixels(visit)

        self.visit_counter -= 1
        pass

    def save(self, filename=None):
        """
        Save the current status of the survey
        """

    def return_survey_map(self):
        return self.survey_map


class countFilterStatus(SurveyStatusSky):
    def __init__(self, nside=64, dtype=int, status_name=None,
                 filter_name=['r'], footprint='D3.5'):
        super(countFilterStatus, self).__init__(nside=nside, dtype=dtype,
                                                status_name=None)
        self.filter_names = filter_name

    def add_visit(self, visit, pixels):
        if visit.filter_name in self.filter_name:
            self.survey_map[pixels] += 1
        self.visit_counter += 1


class CoaddM5Status(SurveyStatusSky):
    def __init__(self, nside=64, dtype=float, status_name=None,
                 filter_name='r', footprint='D3.5'):

        super(CoaddM5Status, self).__init__(nside=nside, dtype=dtype,
                                            status_name=None)
        self.filter_name = filter_name
        self.flux_map = np.zeros(hp.nside2npix(self.nside), dtype=dtype)

    def add_visit(self, visit, pixels=None):
        if pixels is None:
            pixels = self.find_pixels(visit)
        if visit.filter_name == self.filter_name:
            self.flux_map[pixels] += 10.**(0.8*visit.fivesigmadepth)
            self.survey_map[pixels] = 1.25 * np.log10(self.flux_map[pixels])
        self.visit_counter += 1


class HasTemplateStatus(SurveyStatusSky):
    """
    track if there is a good template for a given helpixel
    """
    def __init__(self, nside=64, dtype=bool, status_name=None,
                 filter_name='r', max_seeing=0.8, min_m5=23.,
                 footprint='D3.5'):

        super(HasTemplateStatus, self).__init__(nside=nside, dtype=dtype,
                                                status_name=None)
        self.max_seeing = max_seeing
        self.min_m5 = min_m5
        self.filter_name = filter_name

    def add_visit(self, visit, pixels=None):
        if pixels is None:
            pixels = self.find_pixels(visit)
        if visit.filter_name == self.filter_name:
            if (visit.seeing <= self.max_seeing) & (visit.fivesigmadepth >= self.min_m5):
                self.survey_map[pixels] = True
        self.visit_counter += 1


class LastObserved(SurveyStatusSky):
    def __init__(self, nside=64, dtype=float, status_name=None,
                 filter_name='any', footprint='D3.5'):
        super(LastObserved, self).__init__(nside=nside, dtype=dtype,
                                           status_name=None)
        self.filter_name = filter_name

    def add_visit(self, visit, pixels=None):
        if pixels is None:
            pixels = self.find_pixels(visit)
        if (visit.filter_name == self.filter_name) | (self.filter_name == 'any'):
            self.survey_map[pixels] = np.maximum(self.survey_map[pixels], visit.mjd)
        self.visit_counter += 1


class WeightedSky(BaseSurveyStatus):

    def __init__(self, nside=64, dtype=float, status_name=None):
        super(WeightedSky, self).__init__(nside=nside, dtype=dtype,
                                          status_name=None)

        self.weight_map = np.zeros(hp.nside2npix(self.nside), dtype=float)

    def add_visit(self, visit, pixels=None):
        pass


class HistSky(BaseSurveyStatus):
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
        self.nbins = bins.size
        self.survey_map = np.zeros((hp.nside2npix(self.nside), self.nbins-1), dtype=dtype)

    def add_visit(self, visit, pixels=None):
        pass

    def return_survey_map(self):
        # Collapse it down to 2-D and return
        pass


class RotationHistSky(HistSky):
    """
    Keep track of what rotation angles have been used
    """

    def __init__(self, nside=64, bins=None, filter_names=['u', 'g', 'r', 'i', 'z', 'y'],
                 dtype=float, status_name=None):
        super(RotationHistSky, self).__init__(nside=nside, dtype=dtype,
                                              status_name=None)
        self.filter_names = filter_names
        if bins is None:
            bins = np.linspace(0., 2.*np.pi, 8.)

    def add_visit(self, visit, pixels=None):
        if pixels is None:
            pixels = self.find_pixels(visit)
        if visit.filter_name in self.filter_names:
            index = np.searchsorted(self.bins, visit.rotskypos, side='right')
            self.survey_map[pixels, index] += 1

