import skimage
import pdb
import numpy as np
from scipy import ndimage as ndi
from collections import OrderedDict
import itertools
from copy import deepcopy
import networkx as nx
from scipy.interpolate import splprep, splev, CubicSpline
import pandas as pd
from sklearn.metrics import pairwise_distances
import pdb
from tqdm import tqdm
from p_tqdm import p_map
from sklearn.cluster import DBSCAN
from hdbscan import HDBSCAN
from shapely.geometry import GeometryCollection, LineString, MultiLineString, MultiPoint, Point, Polygon
import matplotlib.pyplot as plt
from sklearn.preprocessing import MinMaxScaler
from skimage.transform import warp

def _idx_to_xy(idx):
	'''
	arg
		idx = 2-tuple of (row, column) indices in an image
	returns
		(x,y) = coordinates of that row and column in xy plane. Note that the positive y axis points down
	'''
	x,y = idx[1], idx[0]
	return x,y

def _xy_to_idx(xy):
	'''
	arg
		xy = 2-tuple, coordinates in xy plane. Note that the positive y axis points down
		(row, col) = 2-tuple of (row, column) indices in an image corresponding to point xy
	'''
	row, col = xy[1], xy[0]
	return row,col




# determine if a given index actually falls within the image
def _idx_in_img(idx, img):
	if 0 <= idx[0] < img.shape[0]:
		if 0 <= idx[1] < img.shape[1]:
			return True
	else:
		return False


def radon_angle(image, theta=None, disable_tqdm = False):
	"""
	Calculates the radon transform of an image given specified
	projection angles.
	Parameters
	----------
	image : array_like
		Input image. The rotation axis will be located in the pixel with
		indices ``(image.shape[0] // 2, image.shape[1] // 2)``.
	theta : array_like, optional
		Projection angles (in degrees). If `None`, the value is set to
		np.arange(180).
	Returns
	-------
	radon_image : ndarray
		Radon transform (sinogram).  The tomography rotation axis will lie
		at the pixel index ``radon_image.shape[0] // 2`` along the 0th
		dimension of ``radon_image``.
	References
	----------
	.. [1] AC Kak, M Slaney, "Principles of Computerized Tomographic
		   Imaging", IEEE Press 1988.
	.. [2] B.R. Ramesh, N. Srinivasa, K. Rajgopal, "An Algorithm for Computing
		   the Discrete Radon Transform With Some Applications", Proceedings of
		   the Fourth IEEE Region 10 International Conference, TENCON '89, 1989
	Notes
	-----
	Based on code of Justin K. Romberg
	(https://www.clear.rice.edu/elec431/projects96/DSP/bpanalysis.html)
	"""
	if image.ndim != 2:
		raise ValueError('The input image must be 2-D')
	if theta is None:
		theta = np.arange(180)

	image = image.astype(float)#convert_to_float(image, preserve_range)
	scaler = MinMaxScaler()
	image = scaler.fit_transform(image) # fits image to interval [0,1]
	image = np.round(image) # now everything is either 0 or 1
	
	# the diagonal length in pixels if the image were square with each
	# side having the length of the largest side
	diagonal = np.sqrt(2) * max(image.shape)
	# how many rows/cols are required to padd the image to make it square 
	# with sidelength of diagonal
	pad = [int(np.ceil(diagonal - s)) for s in image.shape]
	# what the center will be after padding 
	new_center = [(s + p) // 2 for s, p in zip(image.shape, pad)]
	old_center = [s // 2 for s in image.shape]
	pad_before = [nc - oc for oc, nc in zip(old_center, new_center)]
	pad_width = [(pb, p - pb) for pb, p in zip(pad_before, pad)]
	padded_image = np.pad(image, pad_width, mode='constant',
							  constant_values=0)
	def _unpad(img, pad_width):
		return img[
			pad_width[0][0]:-pad_width[0][1],
			pad_width[1][0]:-pad_width[1][1]
		]

	# padded_image is always square
	if padded_image.shape[0] != padded_image.shape[1]:
		raise ValueError('padded_image must be a square')
	center = padded_image.shape[0] // 2
	
	# container for all the angles that have been profiled
	# each entry in this corresponds to a matrix where
	counted_list = [np.ones_like(image) * -1] # the first list of counted is just the background
	angles = [np.nan] # the first angle is background 
	for angle in tqdm(np.deg2rad(theta), disable = disable_tqdm):
		cos_a, sin_a = np.cos(angle), np.sin(angle)
		# note that R is the inverse of the map you want to apply with warp
		R = np.array([[cos_a, sin_a, -center * (cos_a + sin_a - 1)],
					  [-sin_a, cos_a, -center * (cos_a - sin_a - 1)],
					  [0, 0, 1]])
		rotated = np.round(warp(padded_image, R, clip=False))
		
		# labeling different clusters along the profile line
		labeled, num_features = ndi.label(
			rotated,
			structure = np.array([ # this structure means that things are only connected if they in the same column and adjacent
				(0,1,0),
				(0,1,0),
				(0,1,0)
			])
		)
		# now we label the background as -1
		labeled[labeled == 0] = -1
		
		# we only need to count labels within a column because we know that no labels are shared across columns
		for col in range(labeled.shape[1]):
			col_data = labeled[:, col]
			labels, label_counts = np.unique(col_data, return_counts = True)
			# excluding the background label of -1
			label_counts = label_counts[labels > 0]
			labels = labels[labels > 0]
			for n, label in enumerate(labels):
				col_data[col_data == label] = label_counts[n]
				labeled[:, col] = col_data
		# now counted (and labeled) has labeled each cluster by the size of that cluster.
		counted = labeled
		# now we make the matrix to unrotate. This is equivalent to rotating by -angle
		# so cos_a doesn't change but sin_a -> -sin_a
		sin_a *= -1
		R = np.array([[cos_a, sin_a, -center * (cos_a + sin_a - 1)],
					  [-sin_a, cos_a, -center * (cos_a - sin_a - 1)],
					  [0, 0, 1]]
		)       
		# now we unrotate the matrix of cluster sizes
		# each counted array corresponds to an angle.
		# each value in the counted array is the length of the profile line passing through the same pixel in the original image at the given angle
		counted_list.append(
			_unpad(
				warp(
					counted,
					R,
					clip = False,
					cval = -1, # now we use cval of -1 to ensure that background values are not confused with theta = 0
					preserve_range = True # this is important otherwise the values get all fucked
				),
				pad_width = pad_width
			)
		)
		angles.append(angle)
	
	# now we stack the counted_list arrays into a single tensor
	all_counted = np.stack(
		counted_list,
		axis = 2
	)
	max_angle_idx = np.argmax(
		all_counted,
		axis = 2
	)
	# anything where the max_angle_idx is 0 means that is background because the only way this can happen is if 
	# all layers in the stack (along this axis) are -1
	
	angle_img = np.take_along_axis(np.array(angles), max_angle_idx.flatten(), axis = 0).reshape(max_angle_idx.shape)
	angle_img = np.degrees(angle_img)
	angle_img = (angle_img - 90) % 180 # converting from the angle of the normal vector in radon coordinates to the angle of the line segment
	angle_img = 180 - angle_img # changing the direction from cw to ccw
	return angle_img

def line_angle_difference(a1, a2):
	"""
	Calculates the angle difference (in degrees) between two lines at angles a1 and a2 (or pairs of lines if a1,a2 are arrays).
	This accounts for the fact that the max angle difference is 180 because of the twofold degeneracy
	args:
		a1 = angle value(s) of first line(s)
		a2 = angle value(s) of second line(s)
	returns
		delta = the unsigned angle differences between a1 and a2 
	"""
	# ensuring we have arrays
	a1 = np.asarray(a1) % 180
	a2 = np.asarray(a2) % 180
	
	delta1 = (180 - np.abs(a1 - a2)).reshape(-1,1)
	delta2 = np.abs(a1 - a2).reshape(-1,1)
	delta = np.hstack([delta1, delta2]).min(axis = 1)
	if len(delta) == 1:
		return delta.item()
	else:
		return delta
	


def label_nonzero(img):
	'''
	args
		img = np.array
	returns
		labeled = image in which all nonzero values have been reset to be unique and zero values remain zero
	'''
	labeled = np.zeros(img.shape)
	for n, idx in enumerate(np.argwhere(img != 0), start = 1):
		row, col = [int(el) for el in idx]
		labeled[row, col] = n
	return labeled

def diameter_path(g):
	'''
	returns a path between peripheral nodes with path length equal to graph diameter
	note that the diameter path is not necessarily unique
	args
		g = networkx graph. Usually in the form of an skimage RAG
	returns
		path = ordered sequence of nodes in g
	'''
	potential_tips = nx.periphery(g)
	if len(potential_tips) == 2:
		# if only two potential tips, we're done. return these
		t1, t2 = potential_tips
	else:
		diameter = nx.diameter(g)
		#return the first pair we find with distance == diameter
		for t1, t2 in itertools.combinations(potential_tips, 2):
			if nx.shortest_path_length(g, t1, t2) == diameter:
				break
		
	# find shortest path
	return list(nx.shortest_path(g, source = t1, target = t2))



def segment_rag(img, segment):
	'''
	Takes a labeled image where each label corresponds to a segment and produces a rag for the given segment.
	The rag has isolated nodes (pixels) removed and each connected component in the rag is pruned so that only
	the nodes in the diameter path for that component are included in the final rag.

	args
		img = a labeled img
		segment = a single label to be found in img
	returns
		rag = a cleaned RAG corresponding to segment. Note that rag may not be connected. This must be done later.
	'''
	binary = np.where(img == segment, 1, 0)
	binary = np.where(
		skimage.morphology.skeletonize(
			binary,
			method = 'lee'
		),
		1,
		0
	)
	labeled = label_nonzero(binary)
	idx_graph = nx.Graph()
	# add nodes to graph
	# storing the index data temporarily
	for idx in np.argwhere(labeled > 0):
		row, col = [int(el) for el in idx]
		node = labeled[row,col]
		idx_graph.add_node(node)
		idx_graph.nodes[node]['idx'] = tuple(idx)
		idx_graph.nodes[node]['coords'] = _idx_to_xy(idx)
	# creating a rag and passing through the index data
	rag = skimage.future.graph.RAG(
		labeled,
		connectivity = 2,
		data = idx_graph
	)
	# don't forget to remove the 0 background!
	rag.remove_node(0)

	isolates = list(nx.isolates(rag))
	rag.remove_nodes_from(isolates)

	# thinning each component via diameter path to shortest paths
	for component in list(nx.connected_components(rag)):
		component_rag = rag.subgraph(component)
		path = diameter_path(component_rag)

		# find the nodes in component_rag that are not visited in the path
		unvisited_nodes = set(component_rag.nodes) - set(path)
		# remove these nodes from the original rag
		rag.remove_nodes_from( unvisited_nodes )
		

		for node in path:
			assert rag.degree(node) in [1,2], 'Pixel {} in thinned segment {} has degree of not 1 or 2'.format(
				rag.nodes[node]['idx'],
				segment
			)

	return rag


def reconnect_segment_rag(rag):
	'''
	Takes a segment rag and reconnects tips so the rag contains exactly one connected component. 
	Tips are reconnected in order of shortest euclidean distance.	
	args
		rag = a segment rag
	returns
		rag = a segment rag with only one connected component
	'''
	rag = deepcopy(rag)
	tips = [node for node in rag if rag.degree(node) == 1]
	assert len(tips) > 0, 'Number of segment tips ({}) is not greater than 0'.format(len(tips))
	if len(tips) in [1,2]:
		pass # we have nothing to do in this case
	else:
		while not nx.is_connected(rag):
			pair_dists = []
			for n1, n2 in itertools.combinations(tips, 2):
				# only calculate the distance if n1 and n2 are not already connected
				# this excludes nodes that are already on the same segment portion
				# and it precludes comparison of any node with itself
				if not nx.has_path(rag, n1, n2):
					c1, c2 = [np.array(rag.nodes[n]['coords']) for n in [n1,n2]]
					# columns of this will be node1, node2, distance
					pair_dists.append(
						(n1, n2, np.linalg.norm(c1 - c2))
					)
			pair_dists = np.array(pair_dists)
			try:
				min_dist_row_idx = np.argmin(pair_dists[:,2])
			except:
				assert 1==2, 'something horribly wrong here'
			e1, e2 = pair_dists[min_dist_row_idx, 0:2]
			rag.add_edge(e1, e2)
			tips = [el for el in tips if el not in [e1, e2]]
	# now that we're done with the reconnecting, we should have exactly one connected component
	assert (num_comps := nx.number_connected_components(rag)) == 1, 'Graph has {} connected components but should have only 1'.format(num_comps)
	return rag

def cubic_spliner(**kwargs):
	'''
	args
		coords = list of tuples or np.array of ordered coordinates 
		increment = increment when subsetting the coords
	'''


def spliner(coords):
	'''
	args
		coords = list of tuples or np.array. ordered coordinates of line to spline
	returns
		u, tck = tuple to be used in splev
	'''
	# convert to nparray if necessary
	if type(coords[0]) == tuple:
		coords = np.array(coords)
	
	# prepare the spline objects
	if len(coords) <= 10:
		k = 1
		x = [coords[0,0], coords[-1,0]]
		y = [coords[0,1], coords[-1,1]]
		tck, u = splprep(
			[x,y],
			k = k,
			u = list(range(len(x)))
		)	
	else:
		k = 3
		num_anchors = 8 * np.floor(len(coords) / 20 + 1)	# we want 4 anchor points for every 20 pixels
		num_anchors = int(num_anchors)
		anchor_indices = np.linspace(
			start = 0,
			stop = len(coords) - 1,
			num = num_anchors
		)
		anchor_indices = [int(i) for i in anchor_indices]
		anchor_points = [coords[i] for i in anchor_indices]
		x = [c[0] for c in anchor_points]
		y = [c[1] for c in anchor_points]
		tck, u = splprep(
				[x,y],
				s = len(x) - np.sqrt(2 * len(x)), # this smoothing chosen in accordance with manual
				k = k,
				u = list(range(len(x)))
		)
	return (u, tck)

def order_segment_pixels(img_segment_tuple):
	img, segment = img_segment_tuple
	try:
		rag = segment_rag(img, segment)
		rag = reconnect_segment_rag(rag)
		path = diameter_path(rag)
		idx_list = [rag.nodes[node]['idx'] for node in path]
		return segment, idx_list
	except Exception as e:
		msg = 'Error on segment {}: '.format(segment) + e
		print(msg)
		return segment, msg

def labeled_img_to_graph(img, add_electrodes = True, um_per_pix = None, num_cores = 1, increment = 10):
	'''
	args
		img
			labeled and (ideally) skeletonized image of segments
		add_electrodes
			whether to add electrodes on top and bottom
		um_per_pix
			how many microns correspond to a single pixel
		num_cores
			how many cores to use if we want to embarassingly parallelize the most expensive part of this
	returns
		graph of segments with intersections calculated
	'''
	# order segment pixels
	# make splines
	# make graph
	# find intersections

	segment_labels = [n for n in np.unique(img) if n > 0]
	iterable = [(img, segment) for segment in segment_labels]
	def order_segment_pixels(img_segment_tuple):
		img, segment = img_segment_tuple
		try:
			rag = segment_rag(img, segment)
			rag = reconnect_segment_rag(rag)
			path = diameter_path(rag)
			idx_list = [rag.nodes[node]['idx'] for node in path]
			return segment, idx_list
		except Exception as e:
			msg = 'Error on segment {}: '.format(segment) + e
			print(msg)
			return segment, msg
	ordered_segment_pixels_list = p_map(
		order_segment_pixels, 
		iterable, 
		num_cpus = num_cores, 
		desc = 'Ordering Segment Pixels'
	)
	ordered_segment_pixels = {segment : idx_list for (segment, idx_list) in ordered_segment_pixels_list}

	if add_electrodes:
		# add electrodes, making sure that the numbers are how NanowireMesh expects them
		for segment in ordered_segment_pixels.keys():
			assert 0 < segment, 'Segment {} has negative label'.format(segment)
		
		height, width = img.shape
		# make electrodes a little bit inset from the top and bottom
		bottomElectrode = 0
		topElectrode = max(ordered_segment_pixels) + 1
		bottom_row_idx = int(
			round(
				0.99 * (height - 1)
			)
		)
		top_row_idx = int(
			round(
				0.01 * (height - 1)
			)
		)
		ordered_segment_pixels[bottomElectrode] = [(bottom_row_idx, 0), (bottom_row_idx, width - 1)]
		ordered_segment_pixels[topElectrode] = [(top_row_idx, 0), (top_row_idx, width - 1)]

	
	# now make the graph, and insert the splines, linestrings in there
	g = nx.MultiGraph()
	rng = np.random.default_rng()
	for segment, ordered_pixels in ordered_segment_pixels.items():
		try:
			ordered_coords = [_idx_to_xy(idx) for idx in ordered_pixels]
			assert len(ordered_coords) > 1, 'Too few pixels to make segment'
			last_point = ordered_coords[-1]
			ordered_coords = ordered_coords[::increment]
			# making sure last coord of the subset is last coord of the original
			if ordered_coords[-1] != last_point:
				ordered_coords += [last_point]
			ordered_coords = np.asarray(ordered_coords).astype(float)
			# scaling the coordinates to transform from pixel space to 
			# physical space, if necessary
			if um_per_pix is not None:
				ordered_coords *= um_per_pix
				sigma_scale = um_per_pix
			else:
				sigma_scale = 1

			# adding jitter to dramatically decrease the probability that LineStrings
			# cross at anything other than a point
			ordered_coords += rng.normal(
				loc = 0,
				scale = sigma_scale / 100,
				size = ordered_coords.shape
			)
			# t_fit = np.linspace(0,1, ordered_coords.shape[0])
			# interpolator = CubicSpline(
			# 	t_fit, ordered_coords, bc_type = 'natural'
			# )

			# t_eval = np.linspace(0, 1, ordered_coords.shape[0] * 10 * increment)
			# xy_eval = interpolator(t_eval)
			shape = LineString(ordered_coords)
			g.add_node(
				segment,
				idx = ordered_pixels,
				# spline = (u,tck),
				shape = shape, 
			)
		except Exception as e:
			raise Exception(
				'Error on segment {}: {}'.format(segment, e)
			)
	# setting top and bottom electrode attributes (NanowireMesh expects this)
	if add_electrodes:
		g.topElectrode = max(g.nodes)
		g.bottomElectrode = min(g.nodes)
		

	# find intersections
	combos_pbar = tqdm(
		itertools.combinations(g.nodes, 2),
		desc = 'Finding intersections'
	)
	intersection_errors = []
	for n1, n2 in combos_pbar:
		try:
			l1, l2 = [g.nodes[node]['shape'] for node in [n1,n2]]
			X = l1.intersection(l2)
			if X.is_empty:
				pass # don't need to do anything if they don't cross
			else:
				if type(X) == Point:
					g.add_edge(
						n1, n2, 
						shape = X,
						resistanceType = 'cont'
					)
				elif type(X) == MultiPoint:
					for p in X.geoms:
						g.add_edge(
							n1, n2,
							shape = p,
							resistanceType = 'cont'
						)
				else:
					print(type(X))
					raise Exception('Intersection of type {} but must be Point or MultiPoint'.format(type(X)))
		except Exception as e:
			intersection_errors.append(
				'Error calculating intersection between segments {}, {}: {}'.format(n1, n2, e)
			)
	if intersection_errors:
		raise Exception('Intersection errors found\n{}'.format('\n'.join(intersection_errors)))

	for edge in g.edges:
		x,y = g.edges[edge]['shape'].coords[0]
		g.edges[edge]['x'] = x
		g.edges[edge]['y'] = y
		# all resistances at this point are contact resistances
		g.edges[edge]['resistanceType'] = 'cont'
		# but we don't know what resistance values to give them. That will be assigned later
		g.edges[edge]['resistance'] = np.nan
	for node in g.nodes:
		l = g.nodes[node]['shape']
		g.nodes[node]['endpoints'] = [l.coords[0], l.coords[-1]]
		x,y = l.interpolate(0.5, normalized = True).coords[0]
		g.nodes[node]['x'] = x
		g.nodes[node]['y'] = y
		g.nodes[node]['length'] = l.length
		g.nodes[node]['pixel_label'] = node if node not in [g.topElectrode, g.bottomElectrode] else None


	return g

def parallel_edges(g):
	# returns dict of (node1, node2) : [key1, key2, ...]
	return {(edge[0], edge[1]) : list(g[edge[0]][edge[1]].keys()) for edge in g.edges if g.number_of_edges(edge[0], edge[1]) > 1}





def gray_to_rgb(img):
	def _pixel_gray_to_rgb(gray):
		assert gray <= 256**3 - 1, 'Input pixel value must be <= 256^3 - 1'
		r = np.floor(gray / 256**2)
		g = np.floor( (gray - 256**2 * r)/256)
		b = gray - 256**2 * r - 256 * g
		return (r,g,b)

	rescaled = np.round(
		img * (256**3 - 1)/img.max()
	)
	new_img = np.zeros((*img.shape,3))
	for idx in np.argwhere(rescaled > 0):
		row = int(idx[0])
		col = int(idx[1])
		grayval = rescaled[row,col]
		rgb = _pixel_gray_to_rgb(grayval)
		for n in range(3):
			new_img[row,col, n] = rgb[n]

	return new_img



def rgb_to_gray(img):
	def _pixel_rgb_to_gray(rgb):
		r, g, b = rgb
		return 256**2 * r + 256 * g + b 

	r, g, b = [img[:,:,n] for n in range(3)]
	gray_rescaled = _pixel_rgb_to_gray([r,g,b])
	min_nonzero = np.unique(gray_rescaled)[1]
	return np.round(gray_rescaled / min_nonzero)


def nearest_positive_idx(img, x, y):
	row, col = _xy_to_idx(
		np.round([x,y])
	)
	idx = np.array([int(el) for el in [row, col]])
	imgVal = img[idx[0], idx[1]]
	img[idx[0], idx[1]] = 0
	possible_idx = np.argwhere(img > 0)
	# the position of the minimum distance idx in possible_idx
	min_idx_idx = np.linalg.norm(
		possible_idx - idx,
		axis = 1
	).argmin()
	min_idx = possible_idx[min_idx_idx]
	img[idx[0], idx[1]] = imgVal
	return min_idx
