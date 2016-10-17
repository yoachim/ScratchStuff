import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic
from lsst.sims.catalogs.db import CatalogDBObject
from lsst.sims.catUtils.baseCatalogModels import *
from lsst.sims.catUtils.exampleCatalogDefinitions import *
from lsst.sims.utils import ObservationMetaData

# Coordinates at the galactic pole
c = SkyCoord(Galactic, l=0.*u.degree, b=-90*u.degree)

# Radius to query, in degrees
boundLength = 3.5/2.

colnames = ['raJ2000', 'decJ2000', 'umag', 'gmag', 'rmag', 'imag', 'zmag', 'ymag']
constraint = 'rmag < 16'

# dbobj = CatalogDBObject.from_objid('allstars')
dbobj = CatalogDBObject.from_objid('brightstars')

obs_metadata = ObservationMetaData(boundType='circle',
                                   pointingRA=c.icrs.ra.deg,
                                   pointingDec=c.icrs.dec.deg,
                                   boundLength=boundLength, mjd=5700)

t = dbobj.getCatalog('ref_catalog_star', obs_metadata=obs_metadata)

stars = t.db_obj.query_columns(colnames=colnames, obs_metadata=obs_metadata,
                               constraint=constraint, limit=1e6, chunk_size=None)
stars = [chunk for chunk in stars][0]

