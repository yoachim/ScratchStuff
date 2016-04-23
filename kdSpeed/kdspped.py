import numpy as np
import scipy
from scipy.spatial import cKDTree as kdtree

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

	npts = 2e6
	np.random.seed(42)
	# random points on a sphere. Not uniform.
	lat = np.random.uniform(low=-np.pi/2., high=np.pi/2., size=npts)
	lon = np.random.uniform(low=0, high=2*np.pi, size=npts)
	x, y, z = tree_xyz(lat, lon)

	# Postiions to search on
	position_lon = np.linspace(0., 2.*np.pi, 200.)
	position_lat = np.zeros(position_lon.size)

	pos_x, pos_y, pos_z = tree_xyz(position_lat, position_lon)

	leafsize = 100
	tree = kdtree(zip(x,y,z), leafsize=leafsize)

	# Set a 3 degree search radius
	radius = 3.
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