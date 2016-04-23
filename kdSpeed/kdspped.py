import numpy as np
import scipy
from scipy.spatial import cKDTree as kdtree
import healpy as hp

# Check the speed of the scipy kdtree

def test_tree_time():

	def tree_xyz(lat, lon):
		"""
		convert positions on a sphere to 3d coordinates
		"""
		x = np.cos(lat) * np.cos(lon)
		y = np.cos(lat) * np.sin(lon)
		z = np.sin(lat)
		return x, y, z

	npts = 1e6
	np.random.seed(42)
	# random points on a sphere. 
	lat = np.arccos(2.*np.random.uniform(low=0, high=1, size=npts) - 1.) - np.pi/2.
	lon = np.random.uniform(low=0, high=2*np.pi, size=npts)
	x, y, z = tree_xyz(lat, lon)

	# Postiions to search on
	nside = 128
	position_lat, position_lon = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
	position_lat = np.pi/2. - position_lat

	pos_x, pos_y, pos_z = tree_xyz(position_lat, position_lon)

	leafsize = 100
	tree = kdtree(zip(x,y,z), leafsize=leafsize)

	# Set a search radius
	radius = 2.
	x0, y0, z0 = (1, 0, 0)
	x1, y1, z1 = tree_xyz(np.radians(radius), 0)
	search_rad = np.sqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)
	found_points = 0

	for px, py, pz in zip(pos_x, pos_y, pos_z):
		indices = tree.query_ball_point((px, py, pz), search_rad)
		found_points += len(indices)
	print 'found %i points' % found_points

if __name__ == "__main__":
	import timeit
	print 'scipy version', scipy.__version__
	print timeit.timeit("test_tree_time()", setup="from __main__ import test_tree_time", number=3)