import numpy as np
import sqlalchemy as sqla
import os
import healpy as hp


# grab data from median-binned all-sky database

def medDB(where_clause=None, full_select = None, dtypes=None):
    """
    Read in the median binned all-sky camera data

    Parameters
    ----------
    where_clause: str (None)
        The where-cluase of the sql query
    full_select: str (None)
        A full sql-select string
    dtype: list (None)
        A list of dtypes that match the outputs from `full_select`

    Returns
    -------
    data: numpy.array
        A numpy array with the results from the database
    """

    if full_select is None:
        query = 'select hpindex, M, R, G, B, airmass, mjd from medskybrightness where '+ where_clause+';'
    else:
        query = full_select

    if dtypes is None:
        dtypes = [int, float, float, float, float, float,float]
        names = ['hpindex', 'M', 'R', 'G', 'B', 'airmass', 'mjd']
        dtypes = zip(names,dtypes)

    dbAddress = 'sqlite:///meddata.sqlite'
    engine = sqla.create_engine(dbAddress)
    connection = engine.raw_connection()
    cursor = connection.cursor()
    cursor.execute(query)
    data = cursor.fetchall()
    data = np.asarray(data, dtype=dtypes)

    return data

def single_frame(mjd, filter_frame='R', nside=32):
    result = medDB(where_clause = 'mjd = %f' % mjd)
    hpmap = np.zeros(hp.nside2npix(nside), dtype=float)+hp.UNSEEN
    hpmap[result['hpindex']] = result[filter_frame]

    return hpmap
