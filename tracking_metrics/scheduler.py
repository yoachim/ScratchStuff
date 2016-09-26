import numpy as np
import healpy as hp


"""
Take the current sky status and survey status and generate a list of pointings
"""


class scheduler(object):

    def __init__(self, modes=None, nside=64):
        """
        Parameters
        ----------
        modes : list
            List of observingModes objects.
        """

        self.nside = nside
        self.modes = modes
        # collect all the required surveyStatus objects
        self.statusList = []
        for mode in self.modes:
            for key in mode.requiredStatus:
                for status in self.statusList:
                    if mode.requiredStatus[key] == status:
                        # If we already have the status object, point the mode to look at that one
                        mode.requiredStatus[key] = status
                    else:
                        self.statusList.append(mode.requiredStatus[key])

    def add_visit(self, visit):
        """
        Update all the status objects when a new visit is made
        """
        for status in self.statusList:
            status.add_visit(visit)

    def propose_visits(self, observatory_status):

        # Get list of observations from each mode. Maybe also an "urgency" from each mode.
        # if more than one mode is observable, need some logic to pick which one to do, probably 
        # just based on fraction of observations

        pass
        