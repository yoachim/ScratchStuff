import numpy as np
import healpy as hp

# Check to see what area is affected by ending early

data = np.load('/Users/yoachim/Astro/Temp/Science_baseline/FullRange/baseline_FullRange_expMJD_HEAL.npz')

mv = data['metricValues'].copy()
mask = data['mask'].copy()

good = np.where(mask == 0)
mv[np.where(mask == 1)] = hp.UNSEEN

area = hp.nside2pixarea(128, degrees=True)

tenyrarea = np.where(mv > 3652-365)[0].size*area
fiveyrarea = np.where(mv[good] < 365.25*5)[0].size*area
