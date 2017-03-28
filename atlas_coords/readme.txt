Looks like I got astropy.net to work if I used k[1000:1500, 1000:1500] region of the image



Center (RA, Dec):   (30.500, 25.912)
Center (RA, hms):   02h 01m 59.929s
Center (Dec, dms):  +25° 54' 44.736"
Size:   36.2 x 36.2 deg
Radius: 25.583 deg
Pixel scale:    260 arcsec/pixel
Orientation:    Up is 22.1 degrees E of N

http://nova.astrometry.net/user_images/1559387#annotated



---

Trying region hdu = fits.PrimaryHDU(ack[1000:1500,300:800])
Doesn't look good--yup, failed

from astropy.io import fits
import numpy as np
import matplotlib.pylab as plt
hdulist = fits.open('02k57699o0526w.fits.fz')
ack = hdulist[1].data.copy()
ack[np.where(ack < 0)] = 0
hdu = fits.PrimaryHDU(ack[1200:1500,500:800])
hdu.writeto('edge_sub_small.fits')

