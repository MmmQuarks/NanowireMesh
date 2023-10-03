import numpy as np
import networkx as nx
from shapely.geometry import Point, LineString, MultiPoint
from shapely.ops import substring
from scipy.interpolate import splev, splprep
import pandas as pd
from itertools import combinations
import warnings
warnings.filterwarnings('error')
import shapely
assert '1.8.5' in shapely.__version__, 'This package requires shapely version 1.8.5. Current version is {}'.format(shapely.__version__)



def place_wire(
	xy_distribution,
	L_distribution,
	deviation_number_distribution,
	deviation_size_distribution	
):
	'''
	all distributions are argumentless except deviation_size_distribution which takes size as a kwarg
	'''
	# all vectors here are assumed to be row vectors with shape (1,2)
	p1 = xy_distribution() # the first point of the segment
	# storing data for easy sorting later
	df_rows = [
		dict(
			x = p1[0,0],
			y = p1[0,1],
			normalized_distance = 0
		)
	]
	L = L_distribution() # the target length of the segment
	theta = np.random.uniform(low = 0, high = 2 * np.pi)  # angle to travel to get to endpoint of segment
	p2 = p1 + L * np.array(
		[np.cos(theta), np.sin(theta)]
	).reshape(1,2)
	df_rows.append(
		dict(
			x = p2[0,0],
			y = p2[0,1],
			normalized_distance = 1
		)
	)

	# generating deviations from straightness
	deviation_number = deviation_number_distribution()
	deviation_sizes = deviation_size_distribution(
		size = deviation_number
	)
	for deviation_size in deviation_sizes:
		# generate normalized distance along line
		normalized_distance = np.random.uniform(low = 0, high = 1)
		# find point along line
		s = p1 + normalized_distance * (p2-p1)
		# generate deviation point perpendicular to line from p1 to p2
		phi = theta - np.pi / 2
		D = s + deviation_size * np.array([np.cos(phi), np.sin(phi)]).reshape(1,2)
		df_rows.append(
			dict(
				x = D[0,0],
				y = D[0,1],
				normalized_distance = normalized_distance
			)
		)
	df = pd.DataFrame(df_rows).sort_values(
		by = 'normalized_distance', 
		ascending = True
	)
	tck, u = splprep(
		[df.x.values, df.y.values],
		ub = 0,
		ue = 1,
		k = 2 if deviation_number == 1 else 3,
		s = len(df.x) + np.sqrt(2 * len(df.x))
	)
	x_points, y_points = splev(
		np.linspace(0,1,100),
		tck = tck
	)[:2]
	ls = LineString(
		zip(x_points,y_points)
	)
	# assert ls.length >= L
	# # calculating distance to cut from each end
	# # assuming ls.length - 2 * cut = L
	# cut = (ls.length - L)/2
	# ls = substring(
	# 	ls,
	# 	start_dist = cut,
	# 	end_dist = ls.length - cut,
	# 	normalized = False
	# )
	spline = (u,tck)
	return ls, spline
	
def make_curvy_graph(
	num_wires,
	sample_size,
	xy_distribution,
	L_distribution,
	deviation_number_distribution,
	deviation_size_distribution
):
	g = nx.MultiGraph()
	while len(g) < num_wires:
		try:
			ls, spline = place_wire(
				xy_distribution = xy_distribution,
				L_distribution = L_distribution,
				deviation_number_distribution = deviation_number_distribution,
				deviation_size_distribution = deviation_size_distribution
			)
		except RuntimeWarning as e:
			# if we get a splining problem, we just discard this wire and generate a new one
			print('Discarding wire because {}'.format(e))
			continue
		node = len(g) + 1
		g.add_node(
			node,
			spline = spline,
			shape = ls
		)

	# adding electrodes
	electrodes = []
	for idx in range(2):
		x_coords = [0, sample_size[0]] # x coords are the same for both
		y_coords = [idx * sample_size[1]] * 2  # y coords are either 0 or sample height
		tck, u = splprep(
			[x_coords, y_coords],
			ub = 0,
			ue = 1,
			k = 1,
			s = 0
		)
		ls = LineString(
			[(x_coords[n], y_coords[n]) for n in range(len(x_coords))]
		)
		node = min(g.nodes) - 1 if idx == 0 else max(g.nodes) + 1
		electrodes.append(
			(node, dict(shape = ls, spline = (u,tck)))
		)
	g.add_nodes_from(electrodes)
	g.bottomElectrode = electrodes[0][0]
	g.topElectrode = electrodes[1][0]
				

	# adding intersections
	for n1, n2 in combinations(g.nodes, 2):
		l1, l2 = [g.nodes[node]['shape'] for node in [n1,n2]]
		try:
			intersection = l1.intersection(l2)
		except RuntimeWarning as e:
			print('Error with intersections between nodes {},{}: {}'.format(n1,n2,e))
		if intersection.is_empty:
			pass # don't need to do anything if they don't cross
		else:
			if type(intersection) == Point:
				g.add_edge(
					n1, n2,
					shape = intersection,
					resistanceType = 'cont'
				)
			elif type(intersection) == MultiPoint:
				for pt in intersection.geoms:
					g.add_edge(
						n1, n2,
						shape = pt,
						resistanceType = 'cont'
					)
			else:
				raise Exception(
					'Intersection between nodes {},{} is of type {} but only allowed to be Point or MultiPoint'.format(
						n1,
						n2,
						type(intersection)
					)
				)

	# making sure every intersection has what it needs to be processed by to_nanowire_mesh
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
	
	return g
