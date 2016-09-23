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
        pass

    def remove_visit(self, visit, pixels=None):
        if pixels is None:
            pixels = self.observed_healpixels(visit)
        pass

    def save(self, filename=None):
        """
        Save the current status of the survey
        """

    def return_survey_map(self):
        return self.survey_map


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

    def add_visit(self, visit, pixels=None):
        pass

    def return_survey_map(self):
        # Collapse it down to 2-D and return
        pass



