import networkx as nx
from sortedcontainers import SortedDict
from copy import deepcopy
from shapely.ops import substring
from shapely.geometry import Point, LineString
from pandas import DataFrame 
import pickle
from tqdm import tqdm
import pdb
import numpy as np
import sys
sys.path.append('/Users/adamtrebach/Documents/Research/TPV/Nanowires')
import NanowireMesh as nwm

def get_edges(g, node):
	edges = []
	for neighbor in g[node]:
		for key in g[node][neighbor]:
			edges.append(
				(node, neighbor, key) # this order must be preserved. we rely on node to come before neighbor
			)
	return edges

def resultant_segments(g, node):
	if node in [max(g.nodes), min(g.nodes)]:
		return 1
	else:
		edges = get_edges(g,node)
		return max(len(edges), 1)

def scale_lengths(g, pixels, um):
	um_per_pixel = um / pixels
	# nodes
	for node in g:
		shape = g.nodes[node]['shape']
		coords = np.array(shape.coords[:])
		# scaling up coordinates
		coords = um_per_pixel * coords 
		g.nodes[node]['shape'] = shape = LineString(coords)
		# modifying properties
		g.nodes[node]['endpoints'] = [shape.coords[0], shape.coords[-1]]
		g.nodes[node]['x'], g.nodes[node]['y'] = shape.interpolate(distance = 0.5, normalized = True).coords[0]
		g.nodes[node]['length'] = shape.length

	# edges
	for edge in g.edges:
		shape = g.edges[edge]['shape']
		coords = np.array( shape.coords[0] )
		coords = um_per_pixel * coords
		g.edges[edge]['shape'] = shape = Point(coords)
		g.edges[edge]['x'] = shape.x
		g.edges[edge]['y'] = shape.y
		if g.edges[edge]['resistanceType'] == 'int':
			g.edges[edge]['length'] = shape.length



def to_nanowire_mesh(g, pixels, um, nwDiam, initialTemp):
	'''
	Function to convert a multigraph into a NanowireMesh object
	args
		g = multigraph or graph
		pixels = the size (in pixels) of the scale bar
		um = the size (in um) of the scale bar
		nwDiam = nanowire diameter
		initialTemp = system temperature

	the object g must have the following properties
		topElectrode = max(g.nodes)
		bottomElectrode = min(g.nodes)
		all nodes and edges must have 'shape' attribute that is a shapely geometry
	'''
	g = deepcopy(g)
	final_node_count = sum(
		[resultant_segments(g,node) for node in g]
	)
	# checking labels for compatibility with algorithm
	for node in g:
		assert node >= 0, 'Node {} has negative label'.format(node)

	
	# iterating through each node
	topElectrode = max(g.nodes)
	bottomElectrode = min(g.nodes)
	assert bottomElectrode == 0, 'Bottom electrode must be set to 0'
	availNodes = set(g.nodes) - {topElectrode, bottomElectrode}
	# iterating through all nodes except the electrodes
	tqdmAvailNodes = tqdm(availNodes)
	tqdmAvailNodes.set_description('Adding Internal Resistance')
	for node in tqdmAvailNodes:
		original_segment = g.nodes[node]['shape']
		neighbors_df = []
		for edge in get_edges(g,node):
			contact_point = g.edges[edge]['shape']
			dist = original_segment.project(contact_point)
			neighbors_df.append(
				{'dist' : dist, 'edge' : edge}
			)
		# if the node doesn't have to be subdivided, simply relabel it
		if len(neighbors_df) in [0,1]:
			new_label = min(g.nodes) - 1
			g = nx.relabel_nodes(g, {node : new_label}, copy = False)
		else: # if we have subdividing to do (most cases)
			neighbors_df = DataFrame(neighbors_df) # creating df
			neighbors_df = neighbors_df.sort_values( # sorting df by distances
				by = 'dist',
				ascending = True
			).reset_index( # important to reset the index after sorting
				drop = True
			)

			# right now neighbors_df has columns (rows sorted by increasing dist)
			# dist | edge 

	
			# creating the new node labels and adding to the graph
			# remember that each contact neighbor corresponds to one new segment nodes
			new_nodes = []
			for idx, row in neighbors_df.iterrows():
				new_node = min(g.nodes) - 1
				g.add_node(new_node, **deepcopy(g.nodes[node]))
				new_nodes.append(new_node)
			neighbors_df['new node'] = new_nodes

			# right now neighbors_df has columns (rows sorted by increasing dist)
			# dist | edge | new node 
			# new node is the new node label that will soon be connected to the contact neighbor


			# creating a list of distances that will list all the endpoints of the new segments in order
			# segment one will have endpoints at distance 0 and 1,
			# segment two will have endpoints at distances 1, 2 
			# etc
			# we also determine the what the new edge labels will be here
			endpoint_distances = [0]
			internals_df = []
			for idx, row in neighbors_df.iterrows():
				# stop before the last iteration
				if idx == len(neighbors_df) - 1:
					break
				this_edge = row['edge']
				next_edge = neighbors_df.loc[idx+1, 'edge']
				d1 = original_segment.project(g.edges[this_edge]['shape'])
				d2 = original_segment.project(g.edges[next_edge]['shape'])
				endpoint_distances.append(
					(d1 + d2)/2
				)
				this_node = row['new node']
				next_node = neighbors_df.loc[idx + 1, 'new node']
				# calculcating internal resistor positions
				int_x, int_y = original_segment.interpolate((d1 + d2) / 2).coords[0]
				internals_df.append(
					{
						'edge' : (this_node, next_node),
						'length' : d2 - d1,
						'x' : int_x,
						'y' : int_y
					}
				)
			# now have a new df called internals_df that's for the internal resistors. has columns
			# edge | length | x | y
			internals_df = DataFrame(internals_df)

			endpoint_distances.append(original_segment.length)
			# now make the endpoint distances in the nodesdf
			neighbors_df['ep1_dist'] = endpoint_distances[:-1]
			neighbors_df['ep2_dist'] = endpoint_distances[1:]

			# neighbors_df has columns
			# dist | edge | new node | ep1_dist | ep2_dist
			# ep_dist means the distance **along the segment** to this endpoint of the sub segment

			# we now have everything we need to make the new corrected data


			#resetting the new node properties now
			# and adding the contact edges
			for idx, row in neighbors_df.iterrows():
				# deleting inaccurate (and hard to correct) properties
				new_data = g.nodes[row['new node']]
				del new_data['idx']
				del new_data['spline']
				# correcting other properties
				shape = substring(
					new_data['shape'], 
					start_dist = row['ep1_dist'], 
					end_dist = row['ep2_dist']
				)	
				new_data['shape'] = shape
				new_data['endpoints'] = [shape.coords[0], shape.coords[-1]]
				new_data['x'], new_data['y'] = shape.interpolate(distance = 0.5, normalized = True).coords[0]
				new_data['length'] = shape.length

				# adding contact edge
				old_edge = neighbors_df.loc[idx, 'edge']
				old_edge_data = g.edges[old_edge]
				contact_neighbor = old_edge[1]
				g.add_edge(row['new node'], contact_neighbor, **old_edge_data)

			# adding internal resistors
			for idx, row in internals_df.iterrows():
				node1, node2 = row['edge']
				segment_1 = g.nodes[node1]['shape']
				segment_2 = g.nodes[node2]['shape']
				intersection = segment_1.intersection(segment_2).coords[0]
				# sanity check below
				assert intersection == (row['x'], row['y'])
				g.add_edge(
					node1,
					node2,
					length = row['length'],
					shape = Point((row['x'], row['y'])),
					x = row['x'],
					y = row['y'],
					resistanceType = 'int'
				)

			# deleting the original node
			g.remove_node(node)

	# relabeling the top electrode
	nx.relabel_nodes(g, {topElectrode: min(g.nodes) - 1}, copy = False)

	# now relabeling all the nodes to positive integers
	nx.relabel_nodes(
		g,
		{node : -node for node in g},
		copy = False
	)
	assert len(g) == final_node_count, 'Final graph should have {} nodes but actually has {}'.format(final_node_count, len(g))

	# scaling the lengths of everything
	scale_lengths(g, pixels = pixels, um = um)

	# calculating the actual resistance values for internal resistors
	for edge in g.edges:
		if g.edges[edge]['resistanceType'] != 'int':
			continue
		L = g.edges[edge]['length']
		R = nwDiam / 2
		A = np.pi * R**2
		rho = nwm.calc_internal_resistivity(diam = nwDiam, temp = initialTemp)
		g.edges[edge]['resistance'] = rho * L / A

	for edge in g.edges:
		g.edges[edge]['name'] = '_'.join(['r' + g.edges[edge]['resistanceType'], str(edge[0]), str(edge[1])])


	h = nx.Graph(g)
	assert len(h) == len(g), 'Conversion from MultiGraph to Graph failed. MultiGraph has {} nodes while Graph has {}'.format(
		len(g),
		len(h)
	)
	final = nwm.NanowireMesh(
		makeEmpty = False,
		inGraph = h,
		removeNonPercolating = False,
		rcMean = 'bellew'	
	)
	return final