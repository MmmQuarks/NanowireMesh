from scipy.interpolate import splprep, splev
import pdb
import scipy.stats as stats
import shapely
import numpy as np
from shapely.geometry import LineString, MultiPoint, box
from shapely.ops import substring
from numpy.random import default_rng
import geopandas as gpd
import pandas as pd
import os
import random
from pathlib import Path
from matplotlib import pyplot as plt
import pickle
import skimage
from ImageSegmentation import Algorithms as algos

def velocity_linestring(
	width = 100,
	height = 100,
	n_anchors = 2,
	l_loc = 150,
	l_std = 50,
	v_loc = (0,0),
	v_std = 10,
	t0 = 10
):
	"""
	Generates a LineString by starting at some point r[0] then
	
	For n in range(1, n_anchors):
		r[n] = r[n-1] + v[n] * t
		
	where each v[n] is drawn from the normal velocity distribution with mean=v[n-1] and std=v_loc
	We adjust to make the length equal to the randomly drawn length
	
	Then we spline this with degree = n_anchors - 1 and then convert to linestring
	
	args:
		xlim = 2-tuple of xbounds (xmin, xmax)
		ylim = 2-tuple of ybounds (ymin, ymax)
		n_anchors = number of anchor points in each wire
		l_loc = mean of nanowire length distribution
		l_std = standard deviation of nanowire length distribution
		v_loc = mean of velocity distribution
		v_std = stadard deviation of velocity distribution
		t0 = the initial time step used when making the anchor points
	"""
	
	N = n_anchors
		
	# length of curve to be constructed
	length = stats.truncnorm.rvs(
		a = -l_loc/l_std,
		b = np.inf,
		loc = l_loc,
		scale = l_std
	)
	
	xlim = np.asarray([0,width])
	ylim = np.asarray([0,height])
	# starting position
	rng = default_rng()
	r0 = rng.uniform(
		low = (xlim[0], ylim[0]),
		high = (xlim[1], ylim[1]),
		size = (1,2)
	)
	
	# select N-1 velocities
	v = rng.normal(
		loc = v_loc,
		scale = v_std,
		size = (N-1,2)
	)
	# shifting each velocity so it's the same as being drawn from a gaussian centered on the previous velocity
	for row in range(1,len(v)):
		v[row, :] += v[row-1, :]
		
	# making the linestring of indicated length (we initialize an empty linestring of length 0 first)
	ls = LineString()
	t = t0/1.1
	while ls.length < length:
		t *= 1.1 # t starts at t0 then grows by 10% each iteration
		
		# points
		r = [r0]
		for n in range(1,N):
			r.append(
				r[n-1] + v[n-1] * t
			)
		r = np.vstack(r)
		
		#splining
		tck, u = splprep(
			[r[:,0], r[:,1]],
			k = N - 1
		)
		# evaluating spline at 100 points
		x,y = splev(np.linspace(u.min(), u.max(), 100), tck)
		
		# making linestring
		ls = LineString(zip(x,y))
	
	if ls.length > length:
		ls = substring(ls, 0, length)
		
	return ls


def jagged_hull(xlim, ylim, crop = 0.1, N = 3):
	'''
	Creates a random, irregular cropping of the image with limits xlim, ylim by generating N points in regions on left, right, bottom, top of image and then taking convex hull of all these points
	args
		xlim = the min and max possible x coordinates
		ylim = the min and max possible y coordinates
		crop = the size of the region on the edge of the image (region width = 0.1 * image width) for left and right regions, for example
		N = the number of points to generate in each of left, right, top, bottom
	returns
		convex_hull
	'''
	xlim = np.asarray(xlim)
	ylim = np.asarray(ylim)
	width = xlim.max() - xlim.min()
	height = ylim.max() - ylim.min()
	N = int(N)
	# defining regions
	left = np.array(
		[xlim.min(), xlim.min() + crop * width]
	)
	right = np.array(
		[xlim.max() - crop * width, xlim.max()]
	)
	bottom = np.array(
		[ylim.min(), ylim.min() + crop * height]
	)
	top = np.array(
		[ylim.max() - crop * height, ylim.max()]
	)
	regions = [
		zip(left, ylim),
		zip(right, ylim),
		zip(xlim, bottom),
		zip(xlim, top)
	]
	# making the points that will be along the boundary
	points = []
	rng = default_rng()
	for (low, high) in regions:
		points.append(
			rng.uniform(
				low = low,
				high = high,
				size = (N,2)
			)
		)
	points = np.vstack(points)
	return MultiPoint(points).convex_hull

def param_draws(args, num_images):
	"""
	From the options specified by args, this generates and returns a dataframe with num_images draws of parameters
	args
		args = a dataframe of specifying the distributions to draw parameters from
		num_images = how many draws to make
	returns
		darws = dataframe of parameter draws from distributions
	"""
	# there are two types of options in args
	# range options and choices options
	# columns with "range" in the title specify a range of values that can be taken
	# columns with "choices" specify a set of discrete choices that should be chosen from
	# for our purposes, all ranges can be integers except for those specified in non_integer_range_cols
	
	rng = default_rng()
	non_integer_range_cols = ['d_loc_range', 'd_std_range', 'hull_crop_range']
	
	integer_range_cols = [col for col in args.columns if ('range' in col and col not in non_integer_range_cols)]
	draws = pd.DataFrame()
	for col in integer_range_cols:
		draws[col.replace('_range','')] = rng.integers(
			low = args[col].min(),
			high = args[col].max(),
			size = num_images,
			endpoint = True
		)
	
	for col in non_integer_range_cols:
		draws[col.replace('_range', '')] = rng.uniform(
			low = args[col].min(),
			high = args[col].max(),
			size = num_images
		)
	choices_cols = [col for col in args.columns if 'choices' in col]
	for col in choices_cols:
		draws[col] = [[el for el in args[col] if el == el] for n in range(num_images)]
		
	return draws

def make_sample(folder, params):
	"""Creates a sample datum in directory folder using the options in params
	args
		folder = the folder containing the sample data
		params = a row from the params dataframe
	returns
		None
	"""
	# convenience function to create axes of appropriate size
	def _prep_axes():
		fig, ax = plt.subplots(
			figsize = (params['figsize'],) * 2
		)
		ax.set_xlim([0, params['width']])
		ax.set_ylim([0, params['height']])
		ax.set_axis_off()
		return fig, ax
	
	# convenience function to save current figure
	def _save_fig(path):
		plt.savefig(
			path,
			dpi = params['dpi'],
			pad_inches = 0,
			bbox_inches = 'tight',
			transparent = True
		)
		plt.close()
	
	# making background polygon
	background = box(
		minx = 0,
		miny = 0,
		maxx = params['width'],
		maxy = params['height']
	)
	

	# defining boundaries of image using jagged_hull function
	xlim = [0, params['width']]
	ylim = [0, params['height']]
	hull = jagged_hull(
		xlim = xlim,
		ylim = ylim,
		crop = params['hull_crop'],
		N = params['hull_points']
	)
	
	# storing background and hull
	data = dict(
		geometry = [background, hull],
		type = ['background', 'hull'],
		linestring = [np.nan] * 2,
		diam = [np.nan] * 2
	)       
	
	# generating number of anchors that will in each wire
	n_anchors_list = random.choices(
		params['n_anchors_choices'],
		weights = params['n_anchors_choices_weights'],
		k = params['num_wires']
	)
	
	# generating GeoDataFrame of wires
	for wire_num, n_anchors in enumerate(n_anchors_list):
		ls = velocity_linestring(
			width = params['width'],
			height = params['height'],
			n_anchors = n_anchors,
			l_loc = params['l_loc'],
			l_std = params['l_std'],
			v_loc = 0,
			v_std = params['v_std'],
			t0 = params['t0']
		)
		diam = stats.truncnorm.rvs(
			a = -params['d_loc'] / params['d_std'],
			b = np.inf,
			loc = params['d_loc'],
			scale = params['d_std']
		)
		data['geometry'].append(
			ls.buffer(diam / 2)
		)
		data['type'].append('wire')
		data['linestring'].append(ls)
		data['diam'].append(diam)
		
	gdf = gpd.GeoDataFrame(data)
	
	# removing wires that lie outside the hull entirely
	gdf = gdf[gdf.intersects(hull)]
	# cropping to only include wire portions that are within the hull
	wire_idx = gdf['type'] == 'wire'
	gdf.loc[wire_idx, 'geometry'] = gdf.loc[wire_idx, 'geometry'].intersection(hull)
	
	# making the sample folder for this sample's data
	try:
		os.mkdir(folder)
	except FileExistsError:
		# this folder already made. overwrite contents
		pass
	
	# making the sample image
	fig, ax = _prep_axes()
	gdf['color'] = gdf['type'].apply(lambda x : 1 if x == 'wire' else 0)
	gdf.plot(ax = ax, column = 'color', cmap = 'gray')
	_save_fig(Path(folder, 'sample.png'))
	
	# creating the angle image using radon_angle
	# saving as .npy file (so nans are faithfully preserved)
	sample_img = skimage.io.imread(
		Path(folder, 'sample.png')
	)
	sample_img = skimage.color.rgb2gray(sample_img[:,:,:3])
	angle_img = algos.radon_angle(
		sample_img,
		disable_tqdm = True
	)
	np.save(
		Path(folder, 'angle.npy'),
		angle_img
	)
	
	
	# making the hull image
	fig, ax = _prep_axes()
	gdf['color'] = [1 if is_hull else 0 for is_hull in gdf['type'] == 'hull']
	gdf[gdf['type'].isin(['hull', 'background'])].plot(ax = ax, cmap = 'gray')
	_save_fig(Path(folder, 'hull.png'))
	
	# making the mask for each nanowire
	for idx, row in gdf[gdf['type'] == 'wire'].iterrows():
		gdf['color'] = [1 if iidx == idx else 0 for iidx, rrow in gdf.iterrows()]
		fig, ax = _prep_axes()
		gdf.plot(ax = ax, column = 'color', cmap = 'gray')
		_save_fig(
			Path(folder, 'wire_mask_idx_{}.png'.format(idx))
		)
		
	
	# saving geometries. note that we can only save a gdf with one shapely series 
	# first we save all the data along with geometry
	geometry_cols = [col for col in gdf.columns if col not in ['linestring', 'color']]
	gdf[geometry_cols].to_file(
		Path(folder, 'geometry')
	)
	
	# saving linestrings by themselves
	gdf.set_geometry('linestring')[['linestring']].to_file(
		Path(folder, 'linestring')
	)