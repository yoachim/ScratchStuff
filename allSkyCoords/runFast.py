import numpy as np
from astropy.io import fits
import healpy as hp
import ephem

# OK, let's try to read in a fits file, bin it into healpixels, out output the results as fast as possible

# set up the telescope as a global variable
obs = ephem.Observer()
obs.lon = -1.23480997381 #site.longitude_rad
obs.lat = -0.52786436029 #site.latitude_rad
obs.elevation = 2650.0 #site.height
doff = ephem.Date(0)-ephem.Date('1858/11/17')

def raDec2Hpid(nside, ra, dec):
    """
    Assign ra,dec points to the correct healpixel.

    Parameters
    ----------
    nside : int
        Must be a value of 2^N.
    ra : np.array
        RA values to assign to healpixels.
    dec : np.array
        Dec values to assign to healpixels.

    Returns
    -------
    hpids : np.array
        Healpixel IDs for the input positions.
    """
    lat = np.pi/2. - dec
    hpids = hp.ang2pix(nside, lat, ra)
    return hpids

def _hpid2RaDec(nside, hpids):
    """
    Correct for healpy being silly and running dec from 0-180.

    Parameters
    ----------
    nside : int
        Must be a value of 2^N.
    hpids : np.array
        Array (or single value) of healpixel IDs.

    Returns
    -------
    raRet : float (or np.array)
        RA positions of the input healpixel IDs. In radians.
    decRet : float (or np.array)
        Dec positions of the input healpixel IDs. In radians.
    """

    lat, lon = hp.pix2ang(nside, hpids)
    decRet = np.pi/2. - lat
    raRet = lon

    return raRet, decRet

def healbin(ra, dec, values, nside=128, reduceFunc=np.mean, dtype=float):
    """
    Take arrays of ra's, dec's, and value and bin into healpixels. Like numpy.hexbin but for
    bins on a sphere.

    Parameters
    ----------
    ra : np.array
        RA positions of the data points.
    dec : np.array
        Dec positions of the data points.
    values : np.array
        The values at each ra,dec position.
    nside : int
        Healpixel nside resolution. Must be a value of 2^N.
    reduceFunc : function (numpy.mean)
        A function that will return a single value given a subset of `values`.

    Returns
    -------
    mapVals : np.array
        A numpy array that is a valid Healpixel map.
    """

    hpids = raDec2Hpid(nside, ra, dec)

    order = np.argsort(hpids)
    hpids = hpids[order]
    values = values[order]
    pixids = np.unique(hpids)
    pixids = np.arange(hp.nside2npix(nside))

    left = np.searchsorted(hpids, pixids)
    right = np.searchsorted(hpids, pixids, side='right')

    mapVals = np.zeros(hp.nside2npix(nside), dtype=dtype)+hp.UNSEEN

    # Wow, I thought histogram would be faster than the loop, but this has been faster!
    for i, idx in enumerate(pixids):
        mapVals[idx] = reduceFunc(values[left[idx]:right[idx]] )

    # Change any NaNs to healpy mask value
    mapVals[np.isnan(mapVals)] = hp.UNSEEN

    return mapVals


def stupidFast_altAz2RaDec(alt, az, lat, lon, mjd, lmst=None):
    """
    Convert alt, az to RA, Dec without taking into account abberation, precesion, diffraction, ect.
    Parameters
    ----------
    alt : numpy.array
        Altitude, same length as `ra` and `dec`. Radians.
    az : numpy.array
        Azimuth, same length as `ra` and `dec`. Must be same length as `alt`. Radians.
    lat : float
        Latitude of the observatory in radians.
    lon : float
        Longitude of the observatory in radians.
    mjd : float
        Modified Julian Date.
    Returns
    -------
    ra : array_like
        RA, in radians.
    dec : array_like
        Dec, in radians.
    """
    if lmst is None:
        lmst, last = calcLmstLast(mjd, lon)
        lmst = lmst/12.*np.pi  # convert to rad

    dec = np.arcsin(np.sin(lat)*np.sin(alt) + np.cos(lat)*np.cos(alt)*np.cos(az))
    ha = np.arctan2(-np.sin(az)*np.cos(alt), -np.cos(az)*np.sin(lat)*np.cos(alt)+np.sin(alt)*np.cos(lat))
    ra = (lmst-ha)
    raneg = np.where(ra < 0)
    ra[raneg] = ra[raneg] + 2.*np.pi
    return ra, dec



def stupidFast_RaDec2AltAz(ra, dec, lat, lon, mjd, lmst=None):
    """
    Convert Ra,Dec to Altitude and Azimuth.

    Coordinate transformation is killing performance. Just use simple equations to speed it up
    and ignore abberation, precesion, nutation, nutrition, etc.

    Parameters
    ----------
    ra : array_like
        RA, in radians.
    dec : array_like
        Dec, in radians. Must be same length as `ra`.
    lat : float
        Latitude of the observatory in radians.
    lon : float
        Longitude of the observatory in radians.
    mjd : float
        Modified Julian Date.

    Returns
    -------
    alt : numpy.array
        Altitude, same length as `ra` and `dec`. Radians.
    az : numpy.array
        Azimuth, same length as `ra` and `dec`. Radians.
    """
    if lmst is None:
        lmst, last = callablecLmstLast(mjd, lon)
        lmst = lmst/12.*np.pi  # convert to rad
    ha = lmst-ra
    sindec = np.sin(dec)
    sinlat = np.sin(lat)
    coslat = np.cos(lat)
    alt = np.arcsin(sindec*sinlat+np.cos(dec)*coslat*np.cos(ha))
    az = np.arccos((sindec-np.sin(alt)*sinlat)/(np.cos(alt)*coslat))
    signflip = np.where(np.sin(ha) > 0)
    az[signflip] = 2.*np.pi-az[signflip]
    return alt, az



def healbinImage(filename, x,y,alt, az, outfile, nside=32,
                 reduceFunc=np.median, cutoff_level=10.):
    """
    """

    hdulist = fits.open('ut012716.0100.long.M.fits')
    image = hdulist[0].data.copy()
    biasLevel = np.median(image[0:500,0:500])

    image = image - biasLevel

    mjd = hdulist[0].header['MJD-OBS']
    obs.date = mjd - doff
    lst = obs.sidereal_time()

    ra,dec = stupidFast_altAz2RaDec(alt, az, obs.lat, obs.lon, mjd, lmst=lst)

    healMap = healbin(ra,dec, image[y,x],
                      nside=nside, reduceFunc=reduceFunc)


    hpra, hpdec = _hpid2RaDec(nside, np.arange(np.size(healMap)))

    healMap = -2.5*np.log10(healMap)
    cutoff = -2.5*np.log10(cutoff_level)
    healMap[np.isnan(healMap)] = hp.UNSEEN
    healMap[np.where(healMap > cutoff)] = hp.UNSEEN

    good = np.where(healMap != hp.UNSEEN)
    alt, az = stupidFast_RaDec2AltAz(hpra[good], hpdec[good], obs.lat, obs.lon, mjd, lmst=lst)
    airmass = 1./np.cos(np.pi/2-alt)
    hpid = np.arange(np.size(healMap))

    healMap = healMap[good]
    hpid = hpid[good]

    for i,dummy in enumerate(healMap):
         print >>f, '%i, %f, %f, %f' % (hpid[i], healMap[i], airmass[i], mjd)


def checkMap(filename):
    """ Check that I can read the map back in"""
    names = ['hpid', 'rband', 'airmass', 'mjd']
    types = [int,float,float,float]
    ack = np.loadtxt(filename, dtype=zip(names,types), delimiter=',', skiprows=1)
    nside = 32
    healmap = np.zeros(hp.nside2npix(nside), dtype=float)+hp.UNSEEN
    healmap[ack['hpid']] = ack['rband']
    return healmap



if __name__ == '__main__':

    # read in coordinates.  Only do this once since the alt,az coords should be the
    # same in every image
    names = ['x','y','ha','dec']
    dt = [int,int,float,float]
    #XXX Set the single coordinate file to read
    masterCoordFile = 'ut012716.0100.long.M.xxyy'
    masterCoordMJD = 57415.039699
    mjd = masterCoordMJD
    coords = np.loadtxt(masterCoordFile, dtype=zip(names,dt))
    coords['ha'] = np.radians(coords['ha'])
    coords['dec'] = np.radians(coords['dec'])
    # Need to convert this to alt az
    obs.date = mjd - doff
    lst = obs.sidereal_time()
    ra = coords['ha'] - lst
    while ra.min() < 0:
        ra += 2.*np.pi
    ra = ra % (2.*np.pi)

    alt, az = stupidFast_RaDec2AltAz(ra, coords['dec'], obs.lat,
                                     obs.lon, mjd, lmst=lst)
    print 'read in coords'

    # Set up
    f = open('testOut.dat', 'w')
    print >>f, '#healpix id (nside=%i), median R, airmass, mjd '

#    import timeit
#    time = timeit.timeit("healbinImage('ut012716.0100.long.M.fits', coords['x'],coords['y'],alt, az, f)",
#                         number=10, setup="from __main__ import healbinImage, coords,alt,az,f")
#    print 'Time to bin 10 images', time

    #XXX- loop this over all the images. I suppose it might be prudent to make an output file per night
    # to avoid Gb size text files.
    filenames = ['ut012716.0100.long.M.fits']
    for filename in filenames:
        healbinImage(filename, coords['x'],coords['y'], alt, az, f)

    f.close()
