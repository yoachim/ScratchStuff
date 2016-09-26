import numpy as np
import healpy as hp
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
    def __init__(self, ra, dec, filter_names=['u', 'g', 'r', 'i', 'z', 'y'], **kwargs):
        super(DeepDrillingMode, self).__init__(**kwargs)
        for filter_name in filter_names:
            self.requiredStatus[filter_name+'count'] = ss.countFilterStatus(filter_name = filter_name)

    def visits_to_observe(self, observatory_status):
        # Check altitude of field

        # If the altitude is good (and the hour angle is good), propose a serries of visits
        visits = []
        # Make desireability a function of airmass
        desireability = 0.5
        return visits, desireability


def ScanMode(BaseMode):
    def __init__(self, ra, dec, filter_names=['u', 'g', 'r', 'i', 'z', 'y'], **kwargs):
        super(ScanMode, self).__init__(**kwargs)
        for filter_name in filter_names:
            self.requiredStatus[filter_name+'count'] = ss.countFilterStatus(filter_name = filter_name)


def TwilightMode(BaseMode):
    pass

