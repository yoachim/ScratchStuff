import numpy as np
import lsst.sims.featureScheduler as fs
from lsst.sims.speedObservatory import Speed_observatory
import matplotlib.pylab as plt
import healpy as hp

# Try to run rolling cadence

# Make a rolling mask that takes out one-third of the WFD at a time.




if __name__ == '__main__':
    nside = fs.set_default_nside(nside=32)
    survey_length = 8 #365.25  #365.25*10  # days

    dec_limits = np.radians([-16.852, -35.545])

    # Define the mask I want to use
    wfd = fs. WFD_healpixels(nside=nside)
    ra, dec = fs.ra_dec_hp_map(nside=nside)
    rolling_mask1 = np.zeros(ra.size, dtype=bool)
    rolling_mask1[np.where((wfd == 1) & (dec > dec_limits[0]))] = True

    rolling_mask2 = np.zeros(ra.size, dtype=bool)
    rolling_mask2[np.where((dec < dec_limits[0]) & (dec > dec_limits[1]))] = True

    rolling_mask3 = np.zeros(ra.size, dtype=bool)
    rolling_mask3[np.where((dec < dec_limits[1]) & (wfd == 1))] = True

    observatory = Speed_observatory(nside=nside)

    # Define what we want the final visit ratio map to look like
    years = np.round(survey_length/365.25)
    target_map = fs.standard_goals(nside=nside)
    filters = ['u', 'g', 'r', 'i', 'z', 'y']
    surveys = []

    for filtername in filters:
        bfs = []
        bfs.append(fs.M5_diff_basis_function(filtername=filtername, nside=nside))
        bfs.append(fs.Target_map_basis_function(filtername=filtername,
                                                target_map=target_map[filtername],
                                                out_of_bounds_val=hp.UNSEEN, nside=nside))

        bfs.append(fs.North_south_patch_basis_function(zenith_min_alt=50., nside=nside))
        bfs.append(fs.Rolling_mask_basis_function(rolling_mask1, year_mod=3, year_offset=0,
                                                  mjd_start=observatory.mjd, nside=nside))
        bfs.append(fs.Rolling_mask_basis_function(rolling_mask2, year_mod=3, year_offset=1,
                                                  mjd_start=observatory.mjd, nside=nside))
        bfs.append(fs.Rolling_mask_basis_function(rolling_mask3, year_mod=3, year_offset=2,
                                                  mjd_start=observatory.mjd, nside=nside))
        #bfs.append(fs.Zenith_mask_basis_function(maxAlt=78., penalty=-100, nside=nside))
        bfs.append(fs.Slewtime_basis_function(filtername=filtername, nside=nside))
        bfs.append(fs.Strict_filter_basis_function(filtername=filtername))

        weights = np.array([3.0, 0.4, 1., 0., 0., 0., 3., 3.])
        surveys.append(fs.Greedy_survey_fields(bfs, weights, block_size=1, filtername=filtername,
                                               dither=True, nside=nside))

    surveys.append(fs.Pairs_survey_scripted([], [], ignore_obs='DD'))

    # Set up the DD
    dd_surveys = fs.generate_dd_surveys()
    surveys.extend(dd_surveys)

    scheduler = fs.Core_scheduler(surveys, nside=nside)
    
    observatory, scheduler, observations = fs.sim_runner(observatory, scheduler,
                                                         survey_length=survey_length,
                                                         filename='feature_rolling_third_%iyrs.db' % years,
                                                         delete_past=True)

