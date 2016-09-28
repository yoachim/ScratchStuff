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


def m5map2reward(object):
    """
    convert an all-sky five sigma depth map to a reward function
    """
    def __init__(self):
        # load up the cummulative distributions of five-sigma-depths 
        pass

    def return_reward(m5map, filtername):
        """
        
        """
        

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
                 current_m5_weight=0.5, goal_depth_weight=0.5, smooth_fwhm=10.,
                 current_filter_bonus=2., pairs=False, nside=128, **kwargs):
        super(ScanMode, self).__init__(**kwargs)

        self.current_m5_weight = current_m5_weight
        self.goal_depth_weight = goal_depth_weight
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

        # Maybe add a rotation angle check to set the start of the sequence to a rotation angle that
        # helps make things more uniform

        reward_functions = {}
        for filtername in observatory_status.loaded_filters:
            reward_functions[filtername] = m5_maps[filtername]
            # XXX--need to scale it somehow
            # Maybe convert to what percentile depth each healpixel is at in each filter.
            # That gets everything between 0 and 1, and gets things between filters scaled well

            # can make the m5 distribution cumulative, then its an interpolation at every healpixel to
            # get the percentile. Maybe also normalize by fraction of time it is observable at all--nah,
            # look ahead ish things can take care of that I think.
            #

            # Then could convert the coadded depth maps to ~# of observations behind. Normalize by the
            # max behind I guess? Make sure it behaves so everything can go deeper than requested. Maybe put
            # a floor on it too.

        # Add a bonus for the currently loaded filter
        reward_functions[observatory_status.current_filter] *= self.current_filter_bonus




def TwilightMode(BaseMode):
    pass

# Design decisions to make:
# 1) have a mode for each filter, or have the mode pick the filter? I think within the mode makes more sense.
# 2) hard-code visit pairs, or do it with the last-oberved status? Could do ~7.5-minute blocks I suppose.

# * need to have a count for observations within a night. 