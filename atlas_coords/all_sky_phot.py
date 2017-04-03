import numpy as np
import astropy
from astropy.io import fits
import matplotlib.pylab as plt
from astropy.stats import sigma_clipped_stats
from photutils import DAOStarFinder, aperture_photometry, CircularAperture, CircularAnnulus


# Need to pick up photutils
# conda install -c astropy photutils

# Following example from http://photutils.readthedocs.io/en/stable/photutils/aperture.html

# Let's do an example of generating a catalog of stars from an all-sky image

# Collect relevant kwargs

# FWHM of stars
fwhm = 1.5
# how many sigma above background
threshold = 5.
# Radius for stellar aperture
r_star = 4.
# For sky anulus
r_in = 6.
r_out = 8.

# Load up an image and a few header values
hdulist = fits.open('02k57699o0526w.fits.fz')
mjd = hdulist[1].header['MJD-OBS'] + 0
exptime = hdulist[1].header['EXPTIME'] + 0
image = hdulist[1].data + 0.
try:
    hdulist.close()
except:
    hdulist.close()


# crop down image for easy example
image = image[500:1000, 1000:1500]

# Simple stats of the image
maskval = image[0, 0]
good = np.where(image != maskval)
mean, median, std = sigma_clipped_stats(image[good], sigma=3.0, iters=5)

# include a mask
mask = np.zeros_like(image, dtype=bool)
mask[np.where(image == maskval)] = True

# fwhm in pixels, find sources
daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold*std)
sources = daofind(image - median)


# Do aperature phot
positions = (sources['xcentroid'], sources['ycentroid'])
# aperture for stars
apertures = CircularAperture(positions, r=r_star)
# sky background anulus
annulus_apertures = CircularAnnulus(positions, r_in=r_in, r_out=r_out)
apers = [apertures, annulus_apertures]

phot_table = aperture_photometry(image, apers, mask=mask)
bkg_mean = phot_table['aperture_sum_1'] / annulus_apertures.area()
bkg_sum = bkg_mean * apertures.area()
final_sum = phot_table['aperture_sum_0'] - bkg_sum
phot_table['residual_aperture_sum'] = final_sum



from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
norm = ImageNormalize(stretch=SqrtStretch(), vmin=0, vmax=100)
plt.imshow(image, cmap='Greys', origin='lower', norm=norm)
apertures.plot(color='blue', lw=1.5, alpha=0.5)



