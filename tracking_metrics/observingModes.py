import numpy as np
import surveyStatus as ss


def BaseMode(object):
    """
    Class to be the base for the different observing modes
    """

    def __init__(self, **kwargs):
        """
        
        """
        # Make a dict of the status objects that the mode requires to 
        self.requiredStatus = {}

    def visits_to_observe(self, observatory_status):
        pass


def DeepDrillingMode(BaseMode):
    def __init__(self, ra, dec, filter_names=['u', 'g', 'r', 'i', 'z', 'y'], nside=128,
                 **kwargs):
        super(DeepDrillingMode, self).__init__(**kwargs)
        for filter_name in filter_names:
            self.requiredStatus[filter_name+'count'] = ss.countFilterStatus(filter_name=filter_name,
                                                                            nside=nside)
        self.ra = ra
        self.dec = dec

    def visits_to_observe(self, observatory_status):
        # Check altitude of field

        # If the altitude is good (and the hour angle is good), propose a serries of visits
        visits = []
        # Make desireability a function of hour angle.
        desireability = 0.5
        return visits, desireability


def ScanMode(BaseMode):
    """
    Observing mode for raster scanning large contiguous regions of sky. Constructs a map
    of the sky in each filter weighting various considerations.

    Parameters
    ----------
    target_depth_maps : dict of HEALpix arrays
        A dictionary with keys of filter names and values that are
        numpy arrays of a proper healpixel size.
    block_area : float
        The area of sky that should be selected for tesselation (sq degrees)
    smooth_fwhm : float (10.)
        The FWHM of the Gaussian smoothing kernal to use when computing the reward function (degrees).

    """
    def __init__(self, target_depth_maps, block_area=250.,
                 m5_weight=0.5, depth_weight=0.5, smooth_fwhm=10.,
                 pairs=False, nside=128, **kwargs):
        super(ScanMode, self).__init__(**kwargs)

        self.fwhm = np.radians(smooth_fwhm)

        for filter_name in target_depth_maps:
            self.requiredStatus[filter_name+'_m5'] = ss.CoaddM5Status(filter_name = filter_name, nside=nside)

    def visits_to_observe(self, observatory_status, m5_maps):
        """
        Parameters
        ----------
        observatory_status : object
            An object that includes the current filter.
        m5_maps : dict
            A map of the current 5-sigma limiting depth for all filters
        """
        # For each filter, compute the best contour. Find some way to quantify how much I want to observe 
        # each contour--combine how far behind the region is, with the percentile of the observing conditions.

        # with 5-6 contours, pick a contour based on what filter is currently loaded and 
        # how good each contour is and the filter exposure balance.

        # Once map is picked, call code to tesselate the healpixels. Shift pointing around to mimimize
        # power at cosmologically interesting scales


def TwilightMode(BaseMode):
    pass

# Design decisions to make:
# 1) have a mode for each filter, or have the mode pick the filter? I think within the mode makes more sense.
# 2) hard-code visit pairs, or do it with the last-oberved status?

# * need to have a count for observations within a night. 