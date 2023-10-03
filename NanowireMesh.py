import numpy as np
import networkx as nx
import re
import time
import subprocess
import Timer
import DrawFromDistribution as dfd
import math
import statistics as stat
import pickle
import csv
import os
import numexpr as ne
import sys
from scipy.interpolate import griddata
import scipy.integrate as integrate
import scipy
from scipy import optimize, interpolate
import matplotlib.pyplot as plt
import matplotlib.collections as mc
import pylab as pl
from tqdm import tqdm
from copy import deepcopy
from sortedcontainers import SortedDict
import pdb
from scipy.special import jv, jvp, hankel1, h1vp # needed for transparency calculations
import random
import itertools
import types
from shapely.geometry import LineString, Point
import pandas as pd
import geopandas as gpd
from shapely.ops import substring


def calc_internal_resistivity(diam, temp):
	# calculate effective Debye Temperature
	# T1 and T2 are the effective debye temperatures to interpolate between
	# d1 and d2 are the diameters (x coords) corresponding to these temps
	# note that all lengths are in um
	#print(' interpolating to find T_Debye')
	if diam < .030:
		print('Diameter too small to confidently calculate resistivity')
		TD = np.nan
	elif diam >= .030 and diam < .200:
		diams = [.030, .100, .200]
		debyeTemps = [170, 174, 187]
		TD = np.interp(diam, diams, debyeTemps)
	else:
		print('Diameter too large to confidently calculate resistivity (for right now)')
		TD = np.nan
	#print('defining alpha')
	alpha = 4.6 * 10**(-2) # this is just the alpha for electron phonon contribution. units of ohm * um
	# in this whole calculation we are assuming that the residual resistance is negligible compared to electron
	# phonon contribution
	n = 5
	#print('TD = ' + str(TD))
	#print('alpha = ' + str(alpha))
	# below we calculate an integral that returns a tuple of (answer , upper bound on error)
	integral = integrate.fixed_quad(lambda x : x**n / ((np.exp(x) - 1) * (1 - np.exp(-x))), 0, TD/temp)
	rho = alpha * (temp / TD)**n * integral[0]
	# rho returned in units of ohm * um
	return rho




class NanowireMesh(nx.Graph):

	# Initialization Functions
#	def __init__(self, ):
#		super().__init__(self)

	def __init__(self,
		width=150,
		height=150,
		nwLength=10,
		nwLengthSD = 0,
		nwDiam = .15,
		nwDiamSD = 0,
		percMultiple=1.5,
		buffer = 2,
		initialTemp = 298.15,
		inFolder = None,
		inPickle = None,
		inGraph = None,
		makeEmpty = True,
		removeNonPercolating = False,
		addInternalResistance = True,
		useFastJunctionFinder = True,
		rcMean = 50,
		rcSD = 5,
		pdf = lambda x, y : 1,
		angleFunc = lambda : np.random.uniform(low = -np.pi/2, high = np.pi/2), 
		wireDict = None,
		disableTQDM = False):
		# width = sample width
		# height = sample height
		# nwLength = nanowire width
		# percMultiple = normalized nanowire density in terms of the 50% percolation threshold
		# buffer = extra area for which we also generate wires to prevent weird edge effects. I.e. 
		# if buffer = 2 then we generate nanowires over all area PLUS 2 extra nanowire lengths in each direction
		# note that all lengths are in microns
		# angleList is the list of angles to be used for non electrode wires (so angles can be non-random)
		
		#adding all function inputs to the properties dictionary of self
	
		super().__init__(self)

		t0 = time.time()
		

		#list of custom properties that are initialized from / written to csv files in csv function
		self.propsToWrite = ['width', 'height', 'nwLength', 'nwLengthSD', 'percMultiple', 'buffer', 'topElectrode', 'bottomElectrode', 'rcMean', 'rcSD', 'removeNonPercolating']
			
		self.silverSpecificHeat = 235 # J / (kg K)

		# timer object I created to track performance of different parts of the initialization.
		self.timer = Timer.Timer()
		
		# If there is no inFile, we generate a network.
		if makeEmpty:
			pass
		elif inFolder != None and inPickle != None:
			raise Exception('inFolder and inPickle cannot both be used to initialize network. One must be set to None')
		elif inFolder != None and inPickle == None:
			self.init_from_csv(inFolder, disableTQDM = disableTQDM)
		elif inFolder == None and inPickle != None:
			self.init_from_pickle(inPickle, disableTQDM = disableTQDM)
		else:
			# adding all initialization function inputs as properties of the network
			for key,val in locals().items():
				if key not in ['self', '__class__']:
					# we sometimes get errors when pickling if 
					# we store functions as attributes 
					# this fxi is a kludge. We really shouldn't be storing any attributes this way at all
					if type(val) != types.FunctionType:
						setattr(self, key, val)
			self.silverDensity = 1.049 * 10**(-14) # in kg / um^3
			self.silverDensityUnits = 'kg/um^3'

			if inGraph is not None:
				# if inGraph is not None, we have to take special care if wires cross more than once
				# inGraph.edges(node) lists all edges including node. if there are multiedges, 
				# you will get a listing like (0,1), (0,1), (0,2), ...
				# so testing whether a set with no repeated elements has same length as a list with
				# possibly repeating elements tells you whether you have multiple crosses or not
				
				# here we handle multi-edge problems
				# this will require internaledges to be added in the add_internal_resistance function
				# for any segment that is subdivided here
				# use the graph.nodes[node]['pre-segmented'] flag to identify whether a wire was chopped up this way
				# and use graph.nodes[node]['pixel_label'] to determine which segments go together
				inGraph = deepcopy(inGraph)
				edge_counts = dict.fromkeys(
					inGraph.edges(keys = False),
					0
				)
				for edge in inGraph.edges:
					edge_counts[tuple(edge[:2])] += 1
				while max(edge_counts.values()) > 1:
					df = pd.DataFrame()
					df['edges'] = edge_counts.keys()
					df['count'] = edge_counts.values()
					df = df[df['count'] > 1].reset_index(drop = True)
					
					# determine which node you're going to split up
					edge = df.loc[0, 'edges']
					thisnode = edge[0] if edge[0] not in [max(inGraph.nodes), min(inGraph.nodes)] else edge[1]
					ls = inGraph.nodes[thisnode]['shape']
					gdf = pd.DataFrame() # we actually don't need a GeoDataFrame here
					gdf['ewk'] =inGraph.edges(thisnode, keys = True) # ewk means edges with keys
					gdf['edges'] = gdf['ewk'].apply(lambda x : (x[0], x[1]))
					gdf['keys'] = gdf['ewk'].apply(lambda x : x[2])
					gdf['neighbor'] = gdf['edges'].apply(lambda x : x[1])
					gdf['geometry'] = gdf['ewk'].apply(lambda x : inGraph.edges[x]['shape'])
					gdf = gdf.set_geometry('geometry')
					gdf['dist'] = gdf['geometry'].apply(lambda x : ls.project(x))
					gdf = gdf.sort_values('dist').reset_index(drop = True)
					

					# find the midpoints where thisnode will be split
					midpoint_dists = list(0.5 * (gdf.dist.values[:-1] + gdf.dist.values[1:]))
					subsegment_bounds = [0] + midpoint_dists + [ls.length]
					subsegment_starts = subsegment_bounds[:-1]
					subsegment_stops = subsegment_bounds[1:]
					segments = []
					for start, stop in zip(subsegment_starts, subsegment_stops):
						segments.append(
							substring(ls, start, stop)
						)
						
					gdf['segments'] = segments
					gdf['new_node'] = max(inGraph.nodes) + 1 + gdf.index
					

					# adding the new nodes and contact edges to the graph, complete with relevant data
					for idx, row in gdf.iterrows():
						new_data = deepcopy(inGraph.nodes[thisnode])
						# removing data that is no longer accurate
						for key in ['idx']:
							del new_data[key]
						shape = row['segments']
						new_data['parent_shape'] = inGraph.nodes[thisnode]['shape']
						new_data['shape'] = shape
						new_data['endpoints'] = [shape.coords[0], shape.coords[-1]]
						midpoint = shape.interpolate(0.5, normalized = True)
						new_data['x'] = midpoint.x
						new_data['y'] = midpoint.y
						new_data['length'] = shape.length
						new_data['pre-segmented'] = True
						inGraph.add_node(row['new_node'], **new_data)
						
						# now the edge data
						old_edge = row['ewk']
						old_edge_data = inGraph.edges[old_edge]
						inGraph.add_edge(
							row['new_node'],
							row['neighbor'],
							**old_edge_data
						)
						
					# adding the internal resistors between the new nodes
					idf = pd.DataFrame()
					idf['n1'] = gdf.new_node.values[:-1]
					idf['n2'] = gdf.new_node.values[1:]
					idf['length'] = gdf.dist.values[1:] - gdf.dist.values[:-1]
					idf['geometry'] = idf.apply(
						lambda x : inGraph.nodes[x.n1]['shape'].intersection(
							inGraph.nodes[x.n2]['shape']
						),
						axis =1
					)
					# making internal resistor data to add to inGraph
					try:
						intResistorData = idf.apply(
							lambda row : (
								row.n1,
								row.n2,
								dict(
									shape = row.geometry,
									x = row.geometry.x,
									y = row.geometry.y,
									resistanceType = 'int',
									mass = 0,
									temp = initialTemp,
									power = 0,
									length = row['length']
								)
							),
							axis = 1
						)
						inGraph.add_edges_from(list(intResistorData))
					except AttributeError:
						print(row)
					
					
					# delete the old node
					inGraph.remove_node(thisnode)
						
					# now we test to see if we are done
					edge_counts = dict.fromkeys(
						inGraph.edges(keys = False),
						0
					)
					for edge in inGraph.edges:
						edge_counts[tuple(edge[:2])] += 1       
						
				# relabeling all the nodes to be 0, -1, -2, ... -N
				N = len(inGraph) - 1
				non_electrodes = set(inGraph.nodes) - set([inGraph.topElectrode, inGraph.bottomElectrode])
				mapping = {inGraph.topElectrode : -N, inGraph.bottomElectrode : 0}
				for n, node in enumerate(non_electrodes):
					mapping[node] = -n-1
				inGraph = nx.relabel_nodes(inGraph, mapping, copy = False)
				mapping = {node : -node for node in inGraph}
				inGraph = nx.relabel_nodes(inGraph, mapping, copy = False)
				inGraph.topElectrode = max(inGraph.nodes)
				inGraph.bottomElectrode = min(inGraph.nodes)

					
				# add wires from inGraph
				self.add_nodes_from(inGraph.nodes(data = True))
				# add intersections from inGraph
				self.add_edges_from(inGraph.edges(data = True))

				# apply some important properties that are not calculated in the creation of a graph from an image
				for node in self.nodes:
					self.nodes[node]['temp'] = initialTemp
					self.nodes[node]['diam'] = nwDiam
					ls = self.nodes[node]['shape']
					self.nodes[node]['endpoints'] = [ls.coords[0], ls.coords[-1]]
					x,y = ls.interpolate(0.5, normalized = True).coords[0]
					self.nodes[node]['x'] = x
					self.nodes[node]['y'] = y
					self.nodes[node]['length'] = ls.length

				for edge in self.edges:
					self.edges[edge]['x'] = self.edges[edge]['shape'].x
					self.edges[edge]['y'] = self.edges[edge]['shape'].y
					# setting contact resistances
					if self.edges[edge]['resistanceType'] == 'cont':
						if self.rcMean == 'bellew':
							resistance = gen_bellew_resistances(k = 1)
						else:
							resistance = np.random.normal(self.rcMean, self.rcSD)
						self.edges[edge]['resistance'] = resistance


				# setting top electrode and bottom electrode
				self.topElectrode = max(inGraph.nodes)
				self.bottomElectrode = min(inGraph.nodes)

				# marking whether this bad boy is percolating
				self.isPercolating = nx.has_path(inGraph, self.topElectrode, self.bottomElectrode)
				assert self.isPercolating, 'The input graph is not percolating'
			else:
				self.isPercolating = False
				failCount = 0

				while not self.isPercolating:
					self.add_wires(pdf, angleFunc = angleFunc, wireDict = wireDict, disableTQDM = disableTQDM)
					if useFastJunctionFinder:
						self.fast_add_junctions(disableTQDM = disableTQDM)
					else:
						self.robust_add_junctions(disableTQDM = disableTQDM)
					self.isPercolating = nx.has_path(self, self.topElectrode, self.bottomElectrode)
					if not self.isPercolating:
						failCount += 1
						print('Failed to create a percolating network ' + str(failCount) + ' times')
						print('Clearing Network and trying again.')
						self.remove_nodes_from(list(self.nodes))

			# marking the initial number of wires and their initial total length
			wires = set(self.nodes) - {self.topElectrode, self.bottomElectrode}
			self.N0 = len(wires)
			self.lT0 = sum([self.nodes[node]['length'] for node in wires])
			

			# adding internal resistance (only done once we have a percolating network)
			if addInternalResistance:
				self.add_internal_resistance(disableTQDM = disableTQDM, inGraph = inGraph)

			# find the actual nodes in the percolating cluster for future use
			self.find_percolating_cluster()

			if removeNonPercolating:
				print('Removing non-percolating wires.')
				wires = np.array(self.nodes)
				nonPercolatingWires = wires[np.invert(np.isin(wires, self.percolatingCluster))]
				self.remove_nodes_from(nonPercolatingWires)

#		# If there is an inFile, we read from that file to initialize the graph
#		elif type(inFolder) == str and :
#			self.init_from_csv(inFolder, disableTQDM = disableTQDM)
#		else:
#			self.init_from_pickle(inPickle, disableTQDM = disableTQDM)

		self.check_names()
		t1 = time.time()
		print("Runtime: " + str(round(t1-t0,2)))


		# Utility function: initializes the graph by reading from a CSV
	def init_from_csv(self, folderName, disableTQDM = False):
		print(''.join(['Initializing network with data stored in ',folderName]))
		junctionFileName = '/'.join([folderName, 'junctionData.csv'])
		wireFileName = '/'.join([folderName, 'wireData.csv'])
		networkFileName = '/'.join([folderName, 'networkData.csv'])

		# read in junction data
		print('Initializing junction data.')
		file = open(junctionFileName)
		table = list(csv.reader(file, delimiter = ','))
		file.close()
		for counter, row in enumerate(tqdm(table, disable = disableTQDM)):
			if counter == 0:
				attributeKeys = row[2:]
			else:
				w1 = float(row[0]) if isNumber(row[0]) else row[0]
				w2 = float(row[1]) if isNumber(row[1]) else row[1]
				# The above gets wire numbers to floats but allows for non-numeric wire identifiers.

				self.add_edge(w1, w2)
				col = 2
				attributes = {}
				# this adds key value pairs to the attributes dictionary
				for key in attributeKeys:
					attributes.update([(key, float(row[col]) if isNumber(row[col]) else row[col] ) ] )
					col += 1
				self[w1][w2].update(attributes)
		print('Completed.')

		# read in node data
		print('Initializing node data.')
		file = open(wireFileName)
		table = list(csv.reader(file, delimiter = ','))
		file.close()
		counter = 0
		#for counter, row in enumerate(table):
		for counter, row in enumerate(tqdm(table), disable = disableTQDM):
			if counter == 0:
				attributeKeys = row[1:] # all but first col because first col is node number
			else:
				wire = float(row[0]) if isNumber(row[0]) else row[0]
				# add wire to graph
				self.add_node(wire)
				# add attributes to wire
				for col, key in enumerate(attributeKeys, 1):
					value = float(row[col]) if isNumber(row[col]) else row[col]
					self.nodes[wire].update({key : value})
		print('Completed.')

		# reprocess the endpoints data because it was just read in as a string
		print('Reformatting endpoints data from string to list of tuples.')
		ep = nx.get_node_attributes(self, 'endpoints')
		ep = { key : val.replace('[','').replace(']','').replace('(','').replace(')','').split() for key, val in ep.items()}
		ep = {key : {'endpoints' : [(float(val[0].replace(',','')), float(val[1].replace(',',''))), (float(val[2].replace(',','')), float(val[3].replace(',','')))] } for key,val in ep.items()}
		nx.set_node_attributes(self, ep)

		# read in total network data
		print('Fetching simulations parameters.')
		file = open(networkFileName)
		table = list(csv.reader(file, delimiter = ','))
		file.close()
		# table only has two rows
		for rc, row in enumerate(table):
			if rc == 0:
				headers = list(row)
			elif rc == 1:
				values = list(row)
		for i in range(len(headers)):
			value = values[i]
			value = float(value) if isNumber(value) else value
			name = headers[i]
			setattr(self, name, value)
		# we could store the percolating cluster somewhere, but it's SO fast to find it. Doesn't seem necessary.
		print('Completed.')
		self.find_percolating_cluster()
		print('Initialization from files complete.')


	def init_from_pickle(self, inFileName = 'netpickle.p', disableTQDM = False):
		# commenting out the below name check because it certainly seems to cause more harm than good
#		if inFileName[-2::] != '.p':
#			inFileName = inFileName + '.p'
		data = pickle.load( open(inFileName, 'rb'))
		# write global network properties
		for key, val in data['networkProperties'].items():
			setattr(self, key, val)
		# write node properties
		self.add_nodes_from(data['nodes'].keys())
		nx.set_node_attributes(self, data['nodes'])
		#write edge properties
		self.add_edges_from(data['edges'].keys())
		nx.set_edge_attributes(self, data['edges'])
		# finding percolating cluster
		#self.find_percolating_cluster()

		# setting percolating cluster attributes
		try:
			self.percolatingCluster = data['percolatingCluster']
		except KeyError:
			print('No percolating cluster stored in pickle.')
			self.find_percolating_cluster()
		self.isPercolating = self.bottomElectrode in self.percolatingCluster
		

	def mark_to_show(self, nodesToShow):
		for node in self.nodes:
			self.nodes[node]['to_show'] = 1 if node in nodesToShow else 0

	def _get_nodes(self, nodes):
		# gets nodes from network into a list such that they can be added
		# back in with self.add_nodes_from
		# should return a list of nodes in format
		# [ (node, {key : val}), (node2, ...]

		# converting nodes to a list of 1 node
		# if it is not a list already
		try:
			iter(nodes)
		except TypeError:
			nodes = [nodes]

		return [(node, self.nodes[node]) for node in nodes]


	def _get_edges(self, edges):
		# gets edges from network into a list such that
		# they can be added back in with self.add_edges_from
		# should return a list of edges in the format
		# [(node1, node2, {key : val}) ... ]

		# converting edges to proper format
		if type(edges) == tuple:
			edges = [edges]
		else:
			try:
				iter(edges)
			except TypeError:
				raise Exception('Edges must be iterable container of edges or a single eddge')
		return [(edge[0], edge[1], self.edges[edge]) for edge in edges]

	def pop_nodes(self, nodes):
		nodesPopped = self._get_nodes(nodes)
		nodesToRemove = [elem[0] for elem in nodesPopped]
		self.remove_nodes_from(nodesToRemove)
		return nodesPopped

	def pop_edges(self, edges):
		edgesPopped = self._get_edges(edges)
		edgesToRemove = [(elem[0], elem[1]) for elem in edgesPopped]
		self.remove_edges_from(edgesToRemove)
		return edgesPopped

	def pop_nodes_with_edges(self, nodes):
		nodesPopped = self._get_nodes(nodes)
		edges = []
		for nodeData in nodesPopped:
			node = nodeData[0]
			edges = edges + list(self.edges(node))
		edgesPopped = self._get_edges(edges)
		nodesToRemove = [elem[0] for elem in nodesPopped]
		edgesToRemove = [(elem[0], elem[1]) for elem in edgesPopped]
		self.remove_nodes_from(nodesToRemove)
		self.remove_edges_from(edgesToRemove)
		return nodesPopped, edgesPopped

#	def subgraph(self, nodes):
#		print('making subgraph')
#		sg = self.__class__(makeEmpty = True)
#		nodeList = list(self.nbunch_iter(nodes))
#		sg.add_nodes_from(
#			(n, self.nodes[n]) for n in nodeList # using this comprehension w/o enclosing parentheses makes it a generator
#			)
#		if len(nodeList) > 1:
#			sg.add_edges_from(
#				(edge[0], edge[1], self.edges[edge]) for edge in itertools.combinations(nodeList, 2) if edge in self.edges
#				)
#		return sg


		

	# i'm keeping this in here purely for potential future use. 
	# if we ever want to do heat transfer, this may be a useful way of setting up our networks.
	def add_internal_resistance_by_uniform_segmentation(self, disableTQDM = False):
		nonElectrodeNodes = set(self.nodes) - {self.topElectrode, self.bottomElectrode}
		newNodeIndex = -1
		numSubdivisions = 10

		for node in tqdm(nonElectrodeNodes, disable = disableTQDM):
			# subdividing the wire
			nodeStart = np.array(self.nodes[node]['endpoints'][0])
			nodeStop = np.array(self.nodes[node]['endpoints'][1])
			newNodeBounds = np.linspace(nodeStart, nodeStop, numSubdivisions + 1)
			# calculating distances of each of these node bounds from the start point
			distanceBins = [np.linalg.norm(bound - nodeStart) for bound in newNodeBounds]
			nodesToAdd = []
			for n in range(numSubdivisions):
				# deepcopy necessary to stop from assigning same values over again
				newNodeData = deepcopy(self.nodes[node])
				newNodeStart = newNodeBounds[n]
				newNodeStop = newNodeBounds[n+1]
				newNodeData['endpoints'] = [tuple(newNodeStart), tuple(newNodeStop)]
				newNodeData['length'] = np.linalg.norm(newNodeStop - newNodeStart)
				newNodeData['mass'] = self.calc_wire_mass(newNodeData['diam'], newNodeData['length'])
				newNodeData['x'], newNodeData['y'] = (newNodeStart + newNodeStop)/2
				nodesToAdd.append((newNodeIndex, newNodeData))
				newNodeIndex -= 1
			# at this point, all new nodes have been added.

			# now we add internal resistors among the newly created nodes
			edgesToAdd = []
			# total node resistance
			rho = self.calc_internal_resistivity(self.nodes[node]['diam'], self.initialTemp)
			A = np.pi * (self.nodes[node]['diam']/2)**2
			internalResistance = rho * self.nodes[node]['length'] / A
			numInternalResistors = numSubdivisions - 1
			for n in range(numInternalResistors):
				# the internal resistors should be at the mutual points of adjacent wire segments
				newEdge = (nodesToAdd[n][0],  # get node labels to make edge
						nodesToAdd[n+1][0],
						{'x' : newNodeBounds[n+1][0], # setting edge attributes
						'y' : newNodeBounds[n+1][1],
						'resistance' : internalResistance / numInternalResistors,
						'resistanceType' : 'int',
						'mass' : 0,
						'temp' : self.initialTemp,
						'power': 0,
						'length' : self.nodes[node]['length'] / numInternalResistors,
						'diam' : self.nodes[node]['diam'],
						'name' : '_'.join(['rint', str(nodesToAdd[n][0]), str(nodesToAdd[n+1][0])])
						}
					)
				edgesToAdd.append(newEdge)
				
			# now we iterate through the contacts and connect them the appropriate new nodes
			# first we calculate the 
			for neighbor in self.neighbors(node):
				# confirm that we are only dealing with contacts and not internals
				# deepcopy necessary to avoid reassigning the same values over again
				edgeData = deepcopy(self.edges[node, neighbor])
				assert edgeData['resistanceType'] == 'cont', 'All edges should be contact resistors at this point'
				contactPoint = np.array((edgeData['x'] , edgeData['y']))
				# calculate the contact distance from the start point
				contactDist = np.linalg.norm(contactPoint - nodeStart)
				# finding the node that contains this position
				# the -1 at the end is because the bins are labeled with 1 being the first bin
				# and 0 being less than the first bin
				newNodeBin = np.digitize(contactDist, distanceBins) - 1
				try:
					newContactNode = nodesToAdd[newNodeBin][0]
				except IndexError:
					pdb.set_trace()
				# adding the new edge to our list
				edgesToAdd.append((newContactNode,
							neighbor,
							edgeData))

			# remove old node and all associated edges
			self.remove_node(node)
			# add in new sub nodes and edges
			self.add_nodes_from(nodesToAdd)
			self.add_edges_from(edgesToAdd)
			
		# once all nodes have been subdivided, relabel the nodes
		mapping = {node : -node for node in self if node not in [self.topElectrode, self.bottomElectrode]}
		mapping[self.bottomElectrode] = self.bottomElectrode # bottom elec = 0
		mapping[self.topElectrode] = len(mapping) # top elec should be max number 
		self.topElectrode = mapping[self.topElectrode]
		self.bottomElectrode = mapping[self.bottomElectrode]
		nx.relabel_nodes(self, mapping, copy = False)

	def add_internal_resistance(self, disableTQDM = False, inGraph = None):
		print('Adding internal resistance')
		def _clone_node(self, node, newIndex):
			# making the new node data
			newNodeList = [(newIndex, self.nodes[node])]
			self.add_nodes_from(newNodeList)

			# making the new edge data
			neighbors = list(self.neighbors(node))
			newEdgeList = [(newIndex, neighb, self.edges[node, neighb]) for neighb in neighbors]
			self.add_edges_from(newEdgeList)

		# return the coordinates of an edge as nparray
		def _edge_np(edge):
			return np.array( (self.edges[edge]['x'], self.edges[edge]['y']))
	
		# a sorted dictionary with prev and next functions
		class _PathDict(SortedDict):
			def __init__(self,*args, **kwargs):
				super().__init__(*args, **kwargs)

			def next(self, key):
				keysList = list(self.keys())
				index = self.index(key)
				if index + 1 < len(self):
					# returning the value of the next key
					return self[keysList[index + 1]]
				else:
					raise Exception('Error in add_internal_resistance', 'The last element in a dictionary has no next element')

			def prev(self, key):
				keysList = list(self.keys())
				index = self.index(key)
				if index - 1 >= 0:
					# returning the value of the next key
					return self[keysList[index - 1]]
				else:
					raise Exception('Error in add_internal_resistance', 'The first element in a dictionary has no previous element')

		# before we do anything else, record the initial labeling of each node
		for node in self:
			self.nodes[node]['originalLabel'] = deepcopy(node)
			

		availNodes = set(self.nodes) - {self.topElectrode, self.bottomElectrode}
		# the total number of nodes at the end of this function will 2 (for the electrodes) 
		# plus 1 for each node with no neighbors 
		# plus the number of contacts for each node with neighbors
		# this count must exclude nodes that are neighbors via internal resistors
		finalNumberOfNodes = 2
		for node in availNodes:
			contact_edges = [edge for edge in self.edges(node) if self.edges[edge]['resistanceType'] == 'cont']
			finalNumberOfNodes += max(len(contact_edges), 1)
			
			
		# calculating internal resistances for the already extant internal resistors
		# these resistors will only appear when we are converting from image to graph
		# and have to coerce the image to be a graph rather than a multigraph (because curvy
		# wires can cross more than once)
		internal_edges = [edge for edge in self.edges if self.edges[edge]['resistanceType'] == 'int']
		for edge in internal_edges:
			n1,n2 = edge
			diam1 = self.nodes[n1]['diam']
			diam2 = self.nodes[n2]['diam']
			assert diam1 == diam2, "Error on edge {}: diameters of wires participating in internal resistance are unequal".format(edge)
			rho = self.calc_internal_resistivity(diam1, self.edges[edge]['temp'])
			L = self.edges[edge]['length']
			A = np.pi * (diam1 / 2)**2
			self.edges[edge]['resistance'] = rho * L / A

		# make an array of new node indices that we can pop from. We only want non Electrode nodes so
		# we exclude the the largest negative index (top electrode) and 0 (the bottom electrode)
		newNodeIndices = list(range(-finalNumberOfNodes + 1, 0))

		#iterate through the availNodes and divide them, inserting internal resistances
		for node in tqdm(availNodes, disable = disableTQDM):
			# any node here will only have contact edges, so we can assume that all edges are contacts
			# we must handle two cases separately: 1 neighbor or fewer vs more than 1 neighbor

			contact_edges = [edge for edge in self.edges(node) if self.edges[edge]['resistanceType'] == 'cont'] 
			# note that neighbors is really only contact neighbors
			neighbors = [e[1] for e in contact_edges]
			if len(neighbors) <= 1:
				# handle the case of dangling ends and isolates
				# since these don't need any internal resistors, we can just relabel them
				# copy = False sets this to relabel nodes in place
				nx.relabel_nodes(self, {node : newNodeIndices.pop()}, copy = False)
			else:
				# handle the case requiring internal resistors
				# create dict of neighbors sorted by contact distance from node endpoint
				origin = np.array(self.nodes[node]['endpoints'][0])
				if inGraph == None: # i.e. if this isn't graph of curvy wires
					sortedNeighbors = _PathDict([(np.linalg.norm(origin - _edge_np((node, neighb))) , neighb) for neighb in neighbors])
				else: # if this is a curvy graph, we have to be more careful with the distance calculations
					nodeObject = self.nodes[node]['shape']
					distanceTuples = []
					for neighbor in neighbors:
						intersection = self.edges[node,neighbor]['shape']
						distance = nodeObject.project(intersection)
						distanceTuples.append(
							(distance, neighbor)
						)
					sortedNeighbors = _PathDict(distanceTuples)


				sortedNewIndices = _PathDict({dist : newNodeIndices.pop() for dist in sortedNeighbors.keys()})

				# using the structures above, which share the same key, we create new nodes, with identical properties
				# to the original node, except that each of these is only connect to one of the original 
				# node's neighbors
				edgesToAdd = []
				nodesToAdd = []
				for n, key in enumerate(sortedNeighbors):
					neighbor = sortedNeighbors[key]
					newNode = sortedNewIndices[key]
					edgeData = deepcopy(self.edges[node, neighbor]) 
					edgesToAdd.append((newNode, neighbor, edgeData))
					
					newNodeData = deepcopy(self.nodes[node])
					# correcting the node data (endpoints then length and  mass)

					# correcting the first endpoint of the node
					if n != 0:  # this is not done if we are on the first iteration
						# the beginning of this segment should be the end of the previous segment
						# these are stored in the newNodes list which has the form
						# [ (nodeIndex, dictOfAttributes), (nodeIndex2, dictOfAttributes2), ...]
						newNodeData['endpoints'][0] = nodesToAdd[n-1][1]['endpoints'][1]
					
					# correcting the second endpoint of the node
					if n < len(sortedNeighbors) - 1: # this is not done if we are on the last node
						# the new node endpoint should be the midpoint between two adjacent contact points 
						contactPos1 = _edge_np((node, neighbor))
						try:
							nextNeighbor = sortedNeighbors.next(key)
						except ValueError:
							a = sortedNeighbors
							pdb.set_trace()
						contactPos2 = _edge_np((node, nextNeighbor))
						if inGraph == None:
							newNodeData['endpoints'][1] = tuple( 1/2 * (contactPos1 + contactPos2))
							internalResistorLength = np.linalg.norm(contactPos1 - contactPos2)
						else:
							# if we are initializing from an image we must calculate the midpoint 
							# between the contact points more carefully.
							# in this case, this point lies halfway between the contact points 
							# when travelling along the nodeObject
							nodeObject = self.nodes[node]['shape']
							distanceTo1 = nodeObject.project(Point(contactPos1))
							distanceTo2 = nodeObject.project(Point(contactPos2))
							internalResistorLength = np.abs(distanceTo1 - distanceTo2)
							distanceToMidpoint = 1/2*(distanceTo1 + distanceTo2)
							midpoint = nodeObject.interpolate(distanceToMidpoint).coords[0]
							newNodeData['endpoints'][1] = midpoint


						# making the internal resistor between this new node and next new node
						
						intResistorData = {'x' : newNodeData['endpoints'][1][0],
									'y' : newNodeData['endpoints'][1][1],
									'resistanceType' : 'int',
									'mass' : 0,
									'temp' : newNodeData['temp'],
									'power' : 0,
									'length' : internalResistorLength, # this is calculated correctly for both straight and curvy segments
									'diam' : newNodeData['diam']}
						rho = self.calc_internal_resistivity(intResistorData['diam'], intResistorData['temp'])
						L = intResistorData['length']
						A = np.pi * (intResistorData['diam'] / 2)**2
						intResistorData['resistance'] = rho * L / A
						
						# adding to the list of edges to add later
						edgesToAdd.append( (sortedNewIndices[key], sortedNewIndices.next(key), intResistorData) )

					# correcting length
					if inGraph == None:
						ep1 = np.array(newNodeData['endpoints'][0])
						ep2 = np.array(newNodeData['endpoints'][1])
						newNodeData['length'] = np.linalg.norm( ep2 - ep1) 
						newNodeData['x'], newNodeData['y'] = 1/2 * (ep1 + ep2) 
					else:
						# if we are taking curvy segments, we must recalculate the node object
						# and and then use that to calculate lengths and positions (rather than assuming everything is a straight line)
						ep1,ep2 = newNodeData['endpoints']
						# find the distances along the original node
						origNodeObject = self.nodes[node]['shape']
						dist1 = origNodeObject.project(Point(ep1))
						dist2 = origNodeObject.project(Point(ep2))
						# because of the way the preceding code is structured, we are guaranteed that
						# ep1 will lie before ep2 when traveling along origNodeObject.
						# all the same, let's be sure
						assert dist1 < dist2, 'Something appears to have gone wrong with the ordering of junctions along a node object in the add_internal_resistance function'
						# making the segment nodes
						newNodeData['shape'] = substring(origNodeObject, start_dist = dist1, end_dist = dist2)
						newNodeData['length'] = newNodeData['shape'].length
						midpoint = newNodeData['shape'].interpolate(0.5, normalized = True).coords[0]
						newNodeData['x'], newNodeData['y'] = midpoint
					radius = newNodeData['diam'] / 2
					newNodeData['mass'] = self.silverDensity * np.pi * radius**2 * newNodeData['length']

					nodesToAdd.append((newNode, newNodeData))

				self.add_edges_from(edgesToAdd)
				self.add_nodes_from(nodesToAdd)
				

				# asserting that the new wires have lengths that sum to the original length
				oldLength = self.nodes[node]['length']
				newLength = sum([newNode[1]['length'] for newNode in nodesToAdd])
				assert abs(oldLength - newLength)/oldLength <= 0.01, 'Nanowire length is not conserved. Something is wrong'
				
				# removing the node (and all its incident edges) now that we have created wire segments
				self.remove_node(node)
				
		
		# inserting geometries for internal resistances if we have an inGraph
		if inGraph is not None:
			for edge in self.edges:
				if self.edges[edge]['resistanceType'] == 'int':
					s1 = self.nodes[edge[0]]['shape']
					s2 = self.nodes[edge[1]]['shape']
					isect = s1.intersection(s2)
					self.edges[edge]['shape'] = isect
					assert type(isect) == Point, 'Internal resistor {} has type {}'.format(edge, type(isect))
					
	
		# re-numbering all the nodes so they have positive indices
		nonElectrodes = set(self.nodes) - {self.topElectrode, self.bottomElectrode}
		mapping = {node :int( -node) for node in nonElectrodes}
		mapping.update({self.bottomElectrode : int(0),self.topElectrode : len(self.nodes) - 1})
		self.topElectrode = mapping[self.topElectrode]
		nx.relabel_nodes(self, mapping, copy = False)
		# we want to coerce all node labels to be integers, but because 
		for edge in self.edges:
			self.edges[edge]['name'] = '_'.join(['r' + self.edges[edge]['resistanceType'], str(edge[0]), str(edge[1])])

	# making sure that all edges are named correctly and that all node labels are integers
	# note that networkx treats 0.0 and 0 as the same node label, so we have to do an intermediate step
	# to convert from float to int. The easiest way to is map 0.0 -> '0.0' -> 0
	def check_names(self):
		print('Checking that all node names are integers.')
		#mapping1 goes from float -> str
		mapping1 = {}
		for node in self:
			if type(node) != int:
				# if the node can be recast as an integer, we can work with it
				if node == int(node):
					mapping1[node] = str(node)
				else:
					raise Exception(' '.join(['Node', str(node), 'cannot be recast as integer']))
		nx.relabel_nodes(self, mapping1, copy = False)
		#mapping2 goes from str -> int
		mapping2 = {node : int(float(node)) for node in mapping1.values()}
		nx.relabel_nodes(self, mapping2, copy = False)

		# confirm all node labels are integers now
		for node in self:
			assert type(node) == int, ' '.join(['Node', str(node), 'is still not integer even after check_names tried to fix.'])

		# now just for extra special safekeeping, run through all the edge names again
		for edge in self.edges:
			self.edges[edge]['name'] = '_'.join(['r' + self.edges[edge]['resistanceType'], str(edge[0]), str(edge[1])])

# Network Creation Functions
	
	def fast_add_junctions(self, disableTQDM = False):
		# -------- Helper functions -------

		# function to determine if a line segment contains a point
		# that does lie somewhere on the infinite line
		def _contains(segment, point, tolerance = 10**(-10)):
			x1, y1 = segment[0]
			x2, y2 = segment[1]

			x = point[0]
			y = point[1]

			if min(x1,x2) - tolerance <= x <= max(x1,x2) + tolerance:
				if min(y1,y2) - tolerance <= y <= max(y1,y2) + tolerance:
					return True
			return False
		
		# defining function to calculate intersections between line segments
		# linge segments are given as np arrays (or lists of tuples) where each row is one endpoint
		def _find_intersection(seg0, seg1, stepwise = False):
			if stepwise:
				pdb.set_trace()
			tolerance = 10**(-10)
			x1, y1 = seg0[0]
			x2, y2 = seg0[1]

			x3, y3 = seg1[0]
			x4, y4 = seg1[1]

			# make the coefficients matrix A
			# and the output vector b for the equation
			# Ax = b
			A = np.array([ (-(y2 - y1), x2 - x1),
					(-(y4 - y3), x4 - x3)])
			b = np.array([ -(y2 - y1) * x1 + (x2 - x1) * y1,
					-(y4 - y3) * x3 + (x4 - x3) * y3])

			try:
				x, y = np.linalg.solve(A, b)
			except np.linalg.LinAlgError:
				# means lines are parallel
				return None
			else:
				# check to make sure the result is on the finite extent of the lines
				if _contains(seg0, (x,y)) and _contains(seg1, (x,y)):
					return (x,y)
				return []

		# -------- end helper functions
		print('Finding intersections with fast algorithm.')

		# assemble numpy array of node properties
		nodeData = [[node, self.nodes[node]['x'], self.nodes[node]['y'], self.nodes[node]['length']] for node in self.nodes]
		nodeData = np.array(nodeData) 
		# indices of these objects are
		nodeCol, xCol, yCol, lenCol = range(4)

		# tqdm shit
		pbar = tqdm(total = len(nodeData) - 1)
		while len(nodeData) > 1:
			thisNode = nodeData[0]
			nodeData = nodeData[1:]

			# calculating the distances between all the centers
			# of segments
			centerDisplacements = nodeData[:, [xCol, yCol]] - thisNode[[xCol, yCol]]
			centerDisplacements = nodeData[:,[xCol, yCol]] - thisNode[[xCol, yCol]]
			centerDistances = np.sqrt( 
						np.sum( 
							np.square(centerDisplacements),
							axis = 1
							)
						)
			
			# calculating the max distance between centers which is (L1 + L2)/2
			maxDistances = 1/2 * (thisNode[lenCol] + nodeData[:,lenCol])

			# getting the wires that this node may potentially cross
			candidateCrosses = nodeData[:, nodeCol]
			candidateCrosses = candidateCrosses[ centerDistances <= maxDistances]


			# calculating the actual intersections if they exist
			node = thisNode[nodeCol]
			for candidate in candidateCrosses:
				nodeAngle = self.nodes[node]['angle']
				stepwise = False
				if nodeAngle in {0, np.pi/2, -np.pi/2}:
					candAngle = self.nodes[candidate]['angle']
					if candAngle in {0, np.pi/2, -np.pi/2} - {nodeAngle}:
						stepwise = False
				cross = _find_intersection(self.nodes[node]['endpoints'],
							self.nodes[candidate]['endpoints'],
							stepwise = stepwise)

				if cross:
					eps = 10**(-5)
					if 0 - eps <= cross[0] <= self.width + eps and 0 - eps <= cross[1] <= self.height + eps:
						resistance = gen_bellew_resistances(k = 1) if self.rcMean == 'bellew' else np.random.normal(self.rcMean, self.rcSD)
						self.add_edge(node, candidate,
								x = cross[0],
								y = cross[1],
								resistance = resistance,
								resistanceType = 'cont',
								mass = 0,
								temp = self.initialTemp,
								power = 0,
								length = 'N/A',
								diam = 'N/A',
								name = '_'.join(['rcont', str(node), str(candidate)])
								)
			pbar.update(1)
		pbar.close()

					

	def robust_add_junctions(self, disableTQDM = False):
		nodeList = list(self.nodes)
		print('Finding Junctions')
		for n in tqdm(range(len(nodeList)), disable = disableTQDM):
			for m in range(n,len(nodeList)):
				w1 = self.nodes[nodeList[n]]
				w2 = self.nodes[nodeList[m]]
				if w1['angle'] == w2['angle']:
					pass
					# this ensures that we don't get any linear algebra errors if the wires are parallel (which the electrodes are)
				else:
					A = np.array( [ [- np.tan(w1['angle']), 1], [-np.tan(w2['angle']), 1]])
					b = np.array( [-np.tan(w1['angle']) * w1['x'] + w1['y'], -np.tan(w2['angle']) * w2['x'] + w2['y']])
					cross = np.linalg.solve(A,b)
					# test if the cross point is on the sample with some small tolerance for finite precision
					# arithmetic
					if 0 - 0.01 <= cross[0] <= self.width + 0.01:
						if 0 - 0.01 <= cross[1] <= self.height + 0.01:
							ep1 = w1['endpoints']
							ep2 = w2['endpoints']
							xBounds1 = [ep1[0][0], ep1[1][0]]
							xBounds2 = [ep2[0][0], ep2[1][0]] 
							# making sure that the cross lies on the finite extent of the lines
							if min(xBounds1) <= cross[0] <= max(xBounds1):
								if min(xBounds2) <= cross[0] <= max(xBounds2):
									attr = {'x' : cross[0],
										'y' : cross[1]}
									self.add_edge(nodeList[n], nodeList[m], **attr)
		attributes = {edge : {'resistance' : np.random.normal(self.rcMean, self.rcSD), 'resistanceType' : 'cont', 'mass' : 0, 'temp' : self.initialTemp, 'power' : 0, 'length' : 'N/A', 'diam' : 'N/A', 'name' : '_'.join(['rcont', str(edge[0]), str(edge[1])])} for edge in self.edges}
		nx.set_edge_attributes(self, attributes)
						

	def get_number_of_non_electrode_wires(self):
		area  = (self.width + 2* self.buffer * self.nwLength) * (self.height + 2* self.buffer * self.nwLength)
		numWires = int(self.percMultiple * 5.63726 * area / (self.nwLength**2))
		return numWires
		

	# Utility function: Generates the wire positions and angles that will form the network and stores
	# it to the self.wires property.
	def add_wires(self, pdf, angleFunc, disableTQDM = False, wireDict = None):
		print('Making wires')
		# here we make enough wires to cover our sample and the buffer region around it
		if wireDict == None:
			numWires = self.get_number_of_non_electrode_wires()
			self.topElectrode = numWires + 1
			self.bottomElectrode = 0 # below the +1 is to make sure that our non-electrode wires don't have the number 0 (which is reserved for bottom electrode)
			wireNumArray = np.hstack( (np.arange(numWires)+ 1, self.bottomElectrode, self.topElectrode) )
	
			# wire coordinates are generated such that the sample goes from x = 0 to x = height and y = 0 to y = width
	#		xArray = np.hstack( (- self.buffer * self.nwLength + np.random.ranf(numWires) * (self.width + 2* self.buffer * self.nwLength), self.width/2 , self.width/2 ) )
			# note that self.buffer is the number of wirelengths beyond the sample we want to also cover with wires to prevent edge effects
			xBounds = [-self.buffer * self.nwLength, self.width + self.buffer * self.nwLength]
			yBounds = [-self.buffer * self.nwLength, self.height + self.buffer * self.nwLength]
			xArray, yArray = dfd.draw(pdf, xBounds = xBounds, yBounds = yBounds, numDraws = numWires)
			# adding positions of electrodes as last elements of lists
			electrodesX = [self.width/2., self.width/2.]
			electrodesY = [0, self.height]
			xArray = xArray + electrodesX
			yArray = yArray + electrodesY
			# note that we have offset one of the electrodes slightly so that its endpoints will not overlap with the other electrode. This is important for proper functioning of the junction finder
			# ^^^ this doesn't appear to be true anymore
	
	#		yArray = np.hstack( (- self.buffer * self.nwLength + np.random.ranf(numWires) * (self.height + 2* self.buffer * self.nwLength), self.height, 0) )
	
			# making the list of angles using the function angleFunc
			angleArray = np.array([angleFunc() for n in range(numWires)])
	#		if angleList == None:
	#			# if not given a list of angles, draw angles from uniform distribution
	#			angleArray = -np.pi/2 + np.random.ranf(numWires) * np.pi
	##			angleArray = np.hstack( (-np.pi/2 + np.random.ranf(numWires) * np.pi, 0, 0) )
	#			# last two components are rotation angles of electrodes (which are horizontal)
	#		else:
	#			# if given a list of angles, these should be angles of all non electrode wires
	#			if len(angleList) != numWires:
	#				raise Exception('List of angles has the incorrect length. Make sure length of list is equal to output of get_number_of_non_electrode_wires')
	#			angleArray = np.array(angleList)
	
			angleArray = np.hstack((angleArray, 0, 0))
	
			nwLengthArray = self.positive_random_normal(self.nwLength, self.nwLengthSD, len(wireNumArray)-2)
			nwLengthArray = np.hstack((nwLengthArray, self.width, self.width))
			# last two components are widths of electrodes
	
			#nwDiamGlobal = 0.079 this is for Thomas's data
			nwDiamArray = self.positive_random_normal(self.nwDiam, self.nwDiamSD, len(wireNumArray))
			# we're just pretending that the electrodes have the same 'diameter' as the rest of the wires.
			# this property is never used so it won't matter	
			# note that we add a 'tested' column to be used when finding junctions. this column marks whether a wire has already been analyzed to find intersections
			self.wires = [(wireNumArray[n], xArray[n], yArray[n], angleArray[n], nwLengthArray[n], nwDiamArray[n], 0) for n in range(len(wireNumArray)) ]
			self.wires = np.array(self.wires, dtype = [('wireNum', 'f8'), ('x', 'f8'), ('y', 'f8'), ('angle', 'f8'),  ('length', 'f8'), ('diam', 'f8'), ('tested','i4')])
			# now we add these wire properties to the graph
			for wire in tqdm(self.wires, disable = disableTQDM): #i in range(len(self.wires)): # in self.wires:
				x1 = wire['x'] - wire['length'] / 2 * np.cos(wire['angle'])
				x2 = wire['x'] + wire['length'] / 2 * np.cos(wire['angle'])
				y1 = wire['y'] - wire['length'] / 2 * np.sin(wire['angle'])
				y2 = wire['y'] + wire['length'] / 2 * np.sin(wire['angle'])
				
				# approximately eliminating any wires that lie completely out of the sample
				xVals = np.linspace(x1, x2, 10, endpoint = True)
				yVals = np.linspace(y1, y2, 10, endpoint = True)
				xInSample = any(np.logical_and(0 < xVals, xVals < self.width))
				yInSample = any(np.logical_and(0 < yVals, yVals < self.height))
				if (xInSample and yInSample) or (wire['wireNum'] in [self.topElectrode, self.bottomElectrode]):
					# making first pair bottom endpoint and second pair top endpoint
					# the add_internal_resistance function used to depend on this but it 
					# doesn't any longer
					endpoints = [(x1, y1), (x2, y2)] if y1 < y2 else [(x2, y2), (x1, y1)]
					# adding nodes with attributes
					self.add_node(wire['wireNum'], 
						length = wire['length'],
						diam = wire['diam'],
						angle = wire['angle'],
						endpoints = endpoints,
						x = wire['x'],
						y = wire['y'],
						voltage = 0,
						mass = self.calc_wire_mass(wire['diam'], wire['length']),
						temp = self.initialTemp
						)
		else:
			# private class that allows tuples to be sorted
			# first from lowest y to highest y
			# then from lowest x to highest x
			# this allows the code to handle horizontal line segments
			# which it could not do before because it was just sorting by ypositions
			# in a dict but for horizontal line segments all junctions will have same
			# y position
			class _OrderedTuple(tuple):
				def __init__(self, inTuple):
					super().__init__()
	
				def __lt__(self, other):
					if self[1] < other[1]:
						return True
					elif self[1] > other[1]:
						return False
					elif self[1] == other[1]:
						if self[0] < other[0]:
							return True
						elif self[0] > other[0]:
							return False
						elif self[0] == other[0]:
							return False

			print('Using user defined wires rather than generating random')
			print('Verifying that the user defined wires have all the necessary attributes to create network.')
			requiredAttrs = ['length', 'x', 'y', 'angle']
			for key in wireDict.keys():
				for attr in requiredAttrs:
					try:
						temp = wireDict[key][attr]
					except KeyError:
						msg = ' '.join(['Wire with key', str(key), 'lacks required attribute', attr])
						raise Exception(msg)

			print('Finding top and bottom electrodes')
			for key in wireDict.keys():
				try:
					isTop = wireDict[key]['isTopElectrode']
					if isTop:
						self.topElectrode = key
						# we delete this entry so future exports and imports don't get weird
						del wireDict[key]['isTopElectrode']
				except KeyError:
					pass

			for key in wireDict.keys():
				try:
					isBottom = wireDict[key]['isBottomElectrode']
					if isBottom:
						self.bottomElectrode = key
						# we delete this entry so future exports and imports don't get weird
						del wireDict[key]['isBottomElectrode']
				except KeyError:
					pass

			print('Adding automatic attributes for wires')
			newAttrs = {}
			# first we calculate the new properties
			for key in wireDict.keys():
				attr = wireDict[key]
				x1 = attr['x'] - attr['length'] / 2 * np.cos(attr['angle'])
				x2 = attr['x'] + attr['length'] / 2 * np.cos(attr['angle'])
				y1 = attr['y'] - attr['length'] / 2 * np.sin(attr['angle'])
				y2 = attr['y'] + attr['length'] / 2 * np.sin(attr['angle'])
				# making first pair bottom endpoint and second pair top endpoint
				# the add_internal_resistance function depends on this (thought it didn't but it still does) 
				newAttrs[key] = {'endpoints' : [(x1, y1), (x2, y2)] if _OrderedTuple((x1, y1)) < _OrderedTuple((x2,y2))  else [(x2, y2), (x1, y1)],
							'diam' : self.nwDiam,
							'voltage' : 0,
							'mass' : self.calc_wire_mass(self.nwDiam, attr['length']),
							'temp' : attr['temp']}

			# then we add them to the existing wireDict but only if they are not overwriting anything
			for nodeKey in newAttrs.keys():
				for attrKey in newAttrs[nodeKey].keys():
					if attrKey in wireDict[nodeKey].keys():
						print('Calculated attribute was already specified for node', nodeKey, '. Using specified value')
					else:
						wireDict[nodeKey][attrKey] = newAttrs[nodeKey][attrKey]

			# now we can add the wires with all their attributes to the network
			wiresList = [(key, wireDict[key]) for key in wireDict.keys()]
			self.add_nodes_from(wiresList)


			# actually this isn't necessary! The wires will be re-numbered in the add_internal_resistance section
			# making the numbering consistent with expectations
#			print('Re-numbering nodes')
#			nonElectrodeNodes = set(self.nodes) - {self.topElectrode, self.bottomElectrode}
#			mapping = {node : n + 1 for n, node in enumerate(nonElectrodeNodes)}
#			mapping[self.bottomElectrode] = 0
#			self.bottomElectrode = 0
#			mapping[self.topElectrode] = len(self.nodes) - 1
#			self.topElectrode = len(self.nodes) - 1
#			self.relabel(mapping)
		# marking the number of wires at the beginning
		print('Wires completed.')


# Network Modification Functions
	def modify_contact_resistances(self, rcMean, rcSD):
		rand_norm = np.random.normal
		newAttributes = {edge : {'resistance' : rand_norm(rcMean, rcSD)} for edge in self.edges if self.edges[edge]['resistanceType'] == 'cont' }
		nx.set_edge_attributes(self, newAttributes)

	def remove_failing_components(self):
		if 'failedJunctions' not in dir(self):
			self.failedJunctions = {}
		for edge in list(self.edges):
			powerMassRatio = self.edges[edge]['power'] / self.edges[edge]['mass']
			if powerMassRatio >=  0.3 * 10**12:
				self.failedJunctions.update({edge : self.edges[edge]})
				self.remove_edge(*edge)
		self.find_percolating_cluster()
		print('Number of failed junctions =', len(self.failedJunctions))

	def add_shunt_resistors(self, rShunt = 10**(12)):
		# making a list of the dangling ends
		# but excluding dangling ends that have the bottom electrode as their 
		# sole neighbor
		danglingEnds = [node for node in self if self.degree(node) == 1 and self.bottomElectrode not in self.neighbors(node)]
		for node in danglingEnds:
			self.add_edge(node, self.bottomElectrode,
					name = '_'.join(['rshunt', str(int(node)), str(int(self.bottomElectrode))]),
					resistance = rShunt,
					resistanceType = 'shunt')

	def remove_shunt_resistors(self):
		shuntEdges = [edge for edge in self.edges if self.edges[edge]['resistanceType'] == 'shunt']
		self.remove_edges_from(shuntEdges)

	def solve_circuit_using_xyce(self,
					xycePath = 'xyce', 
					netlistName = 'netlist',
					voltage = 1,
					verbosity = 'WARN'):
		assert verbosity in ['DEBUG', 'INFO', 'WARN', 'ERROR'], 'unrecognized option {} for verbosity. Must be one of DEBUG, INFO, WARN, ERROR'.format(verbosity)
		print('Solving circuit using xyce')

		# the steps here are 
		# 1. write the netlist
		# 2. solve with xyce
		# 3. update network with xyce output
		
		# 1. Writing the netlist
		# first make a list where each element corresponds
		# to a line of the netlist

		# some safety procedures to make sure that
			# the network has a path from top to bottom
			# all the node labels are integers
			# the network does not include non percolating wires when written to netlist because
			# this will cause a xyce error.
			# all node names are integers
		# making sure top and bottom electrodes have a path between them
		if not nx.has_path(self, self.topElectrode, self.bottomElectrode):
			raise Exception('Top and bottom electrodes are not connected. Aborting solve.')
		# making sure that all nodes have integer labels
		self.check_names()
		self.find_percolating_cluster()
		nonPercolatingNodes = set(self.nodes) - set(self.percolatingCluster) # note we don't have to remove dangling ends
		# maybe we do have to remove dangling ends?
		poppedNodes, poppedEdges = self.pop_nodes_with_edges(nonPercolatingNodes)
#		danglingEnds = [node for node in self if self.degree(node) == 1]
#		while danglingEnds:
#			pNodes, pEdges = self.pop_nodes_with_edges(danglingEnds)
#			poppedNodes += pNodes
#			poppedEdges += pEdges
#			danglingEnds = [node for node in self if self.degree(node) == 1]
#
		# making header and voltage source
		netlist = ['Header Row']

		# adding shunt resistors where necessary
#		self.add_shunt_resistors()

		# voltage source line should have structure
		# VIN topElectrode bottomElectrode DC voltage
		voltageIn = ' '.join(['VIN', str(int(self.topElectrode)), str(int(self.bottomElectrode)), str(voltage)])
		netlist.append(voltageIn)

		# making the resistors which should have structure
		# name node1 node2 resistance
		# note that name must begin with r or R
		resistors = [' '.join([self.edges[edge]['name'], 
				str(int(edge[0])),
				str(int(edge[1])),
				str(self.edges[edge]['resistance'])]) for edge in self.edges] 
		netlist = netlist + resistors

		# adding operating point analysis
		netlist.append('.OP')

		# adding print line
		printLine = ['.PRINT DC FORMAT=CSV']
		printPowers = [''.join(['P(', self.edges[edge]['name'],')']) for edge in self.edges]
		printVoltages = [''.join(['V(', str(int(node)),')']) for node in self.nodes]
		printCurrent = ['I(VIN)']
		printLine = printLine + printPowers + printVoltages + printCurrent
		printLine = ' '.join(printLine)
		netlist.append(printLine)

		# adding end line
		netlist.append('.END')

		with open(netlistName, 'w') as fileHandle:
			fileHandle.writelines("%s\n" % line for line in netlist)

		# 2. solving circuit with xyce
		command = [xycePath, netlistName, '-l', verbosity]
		#cp is a CompletedProcess object
		cp = subprocess.run(command, 
					capture_output = True)
		# logging output 
		outputName = netlistName + '_out'
		with open(outputName, 'w') as fileHandle:
				fileHandle.write(cp.stdout.decode('utf-8'))
		print('subprocess stdout data')
		print(cp.stdout.decode('utf-8'))
		print('subprocess stderr data')
		print(cp.stderr.decode('utf-8'))

		#logging error only if error is returned
		if cp.returncode != 0:
			errName = netlistName + '_err'
			with open(errName, 'w') as fileHandle:
				fileHandle.write(cp.stdout.decode('utf-8'))
			exceptionMsg = 'Xyce Error. See ' + outputName + ' or ' + errName + ' for more information.'
			raise Exception(exceptionMsg)

		# 3. updating network with xyce output
		# first open the file and convert to a list
		file = open(netlistName + '.csv')
		table = list(csv.reader(file, delimiter = ','))
		file.close()

		# table has two rows.
		# first row is the name of the value, second row is the value
		# all names will be either
		# P(rtype_node_othernode) OR V(node) OR I(VIN)

		# making a mapping so it's easy to assign values back to the right edges
		# when reading from the xyce output csv
		nameToEdgeMapping = {self.edges[edge]['name'].upper(): edge for edge in self.edges}
		# note that we have converted all edge names to upper case because that's how they'll
		# be printed in Xyce's output


		printNames = table[0]
		vals = table[1]
		for col, printName in enumerate(printNames):
			# getting string inside parentheses
			insideParens = re.search(r'\((.*)\)', printName)
			insideParens = insideParens.group(0)
			# removing parentheses
			elementName = insideParens[1:-1]
			if elementName[0] == 'R':
				# then this is a power through resistor
				# marked as P(RTYPE_NODE_OTHERNODE)
				edge = nameToEdgeMapping[elementName]
				val = float(vals[col])
				self.edges[edge]['power'] = val
			elif elementName[0] == 'V':
				# then this is the current through the whole system
				# marked as I(VIN)
				current = float(vals[col])
			else:
				# if name is neither of the above, it must be a voltage
				# marked as V(NODE)
				# we check this and raise an exception if
				# the name is not a voltage
				if printName[0] != 'V':
					exceptionMsg = ' '.join(['Unexpected value found in column',
								str(col),
								'of xyce output file',
								netlistName + '.csv'])
					raise Exception(exceptionMsg)
				node = int(elementName)
				val = float(vals[col])
				self.nodes[node]['voltage'] = val

		# setting the sheet resistance of the system
		self.sheetResistance = voltage / abs(current)

		# removing the shunt resistors we added for xyce
#		self.remove_shunt_resistors()
		
		# calculating the currents for edges
		for edge in self.edges:
			power = self.edges[edge]['power']
			resistance = self.edges[edge]['resistance']
			current = np.sqrt(power / resistance)
			self.edges[edge]['current'] = current

		# adding back in the popped nodes and edges
		self.add_nodes_from(poppedNodes)
		self.add_edges_from(poppedEdges)

		print('Xyce successfully solved circuit')

	# updates self with properties read from output of xyce simulations
	def update_with_xyce_output(self, inFile, disableTQDM = False):
		data = read_xyce_output(inFile, disableTQDM = disableTQDM, topElectrode = self.topElectrode)
		self.sheetResistance = data['sheetResistance']

		print('Writing node data')
		for key in tqdm(data['nodes'], disable = disableTQDM):
			nx.set_node_attributes(self, {key : data['nodes'][key]})

		print('Writing edge data')
		for key in tqdm(data['edges'], disable = disableTQDM):
			nx.set_edge_attributes(self, {key : data['edges'][key]})
		print('Network update complete.')

	def propagate_time(self, tstep, option):
		sqrt = np.sqrt
		cbrt = np.cbrt
		def T_root(i, p1, p2):
			a = 4 * p1 * cbrt(2/3)
			# note: b/c ignoring conduction, the argument of cbrt is guaranteed to be real. We take the real part below 
			# just so numpy doesn't freak out
			b = cbrt( np.real(sqrt(3) * sqrt(27 * p2**4 - 256 * p1**3 + 0j) + 9 * p2**2  ))
			c = a / b + b / cbrt(18)
			s = 1 if i in [1,2] else -1
			return s * 1/2 * sqrt(c + 0j) + (-1)**i * 1/2 * sqrt(-c - s * 2 * p2 / sqrt(c + 0j) + 0j)
		def find_new_temp(node, p1 , p2 , p3p, option, tstep = tstep):
			# note p1, p2 are unprimed but p3p is primed.
			def ugly_prod(T):
				return np.prod([T - T_root(mu, p1, p2) ** (1/ np.prod([T_root(mu, p1, p2) - T_root(gamma, p1, p2) for gamma in [1,2,3,4] if gamma != mu])) for mu in [1,2,3,4] ])
			T_old = self.nodes[node]['temp']
			m = self.nodes[node]['mass']
			C = self.silverSpecificHeat
			F = lambda T : np.exp(p3p * tstep / (m * C)) * ugly_prod(T_old) - ugly_prod(T)
			if option ==1:
				return optimize.broyden1(F, T_old)
			elif option == 2:
				dTdt = p3p/(m * C) * (p1 + p2 * T_old + T_old**4)
				if dTdt <0:
					print('cools down at node:', node)
				return T_old + dTdt * tstep
		# set constants
		sigma = 5.67 * 10**(-20) # stefan-boltzmann in W / (um^2 K^4)
		eps_NW = 0.025 # silver nanowire emissivity
		T_air = 298.15
		h_air = 50 * 10**(-12) # heat transfer coefficient in W / (um^2 K)

		# iterating over all nodes with temperatures (not electrodes)
		print('Solving for Temperatures')
		newNodeTemps = {}
		for node in tqdm(self.percolatingCluster, disable = disableTQDM):
			if node not in [self.topElectrode, self.bottomElectrode]:
				A_surf =  np.pi * self.nodes[node]['diam'] * self.nodes[node]['length']
				
				p3p = -sigma * eps_NW * A_surf
				# print('p3p = ',p3p)
				#calculating p1
				jouleHeating = 1/2 * sum([self.edges[node, neighbor]['power'] for neighbor in self.neighbors(node)])
				# print(jouleHeating)
				convectionFromAir = 1/2 * h_air * A_surf * T_air
				blackbodyFromAir = sigma * eps_NW * A_surf * T_air **4
				p1 = (jouleHeating + convectionFromAir + blackbodyFromAir) / p3p

				#calculating p2
				p2 = h_air/(2 * sigma * eps_NW)

				# solve ugly fucker
				# print('eval ugly fucker')
				# Tlist = find_new_temp(node, p1, p2, p3p, option = option)
				Tlist = [T_root(i, p1, p2) for i in [1,2,3,4]]
				Tlist = [temp for temp in Tlist if (np.imag(temp) == 0 and np.real(temp) >= 0)]
				newTemp = np.real(Tlist[0])
				newNodeTemps.update({node : {'temp' : newTemp}})
		# recording old node temperatures for use when updating contact resistances
		oldNodeTemps = {node : self.nodes[node]['temp'] for node in self.nodes}
		# writing new node temperatures
		nx.set_node_attributes(self, newNodeTemps)
		# remove failing nodes
		print('Removing Failing Nodes')
		failedNodes = [node for node in self.nodes if self.nodes[node]['temp'] >= 473.15]
		if len(failedNodes) > 0:
			self.remove_nodes_from(failedNodes)

		print('Updating Electrical Properties')
		newProps = {}
		for edge in tqdm(self.edges(data = True), disable = disableTQDM):
			edgeTemp = np.average([self.nodes[edge[0]]['temp'], self.nodes[edge[1]]['temp'] ])			
			if edge[2]['resistanceType'] == 'int':
				rho = self.calc_internal_resistivity(edge[2]['diam'], edgeTemp)
				newResistance = rho * edge[2]['length'] / (np.pi * (edge[2]['diam']/2)**2 )
			elif edge[2]['resistanceType'] == 'cont':
				oldEdgeTemp = np.average([oldNodeTemps[edge[0]], oldNodeTemps[edge[1]] ])
				edgeDiam = np.average([self.nodes[edge[0]]['diam'], self.nodes[edge[1]]['diam'] ])
				rho = self.calc_internal_resistivity(edgeDiam, edgeTemp)
				oldRho = self.calc_internal_resistivity(edgeDiam, oldEdgeTemp)
				newResistance = rho / oldRho * self.edges[edge[0],edge[1]]['resistance']
			newProps.update({(edge[0], edge[1]) : {'resistance' : newResistance}})
			if newResistance < 0:
				print('fucked up resistance', edge[0], edge[1])
		nx.set_edge_attributes(self, newProps)
		self.find_percolating_cluster()


# visualizations
	def plot_gnp(self, fileName = 'plot',
		show = True,
		resolution = 8192):


		dataFileName = r'_'.join([fileName, 'data'])
		scriptFileName = r'_'.join([fileName, 'script'])

		# making data file for gnuplot to read 
		with open(dataFileName, 'w') as dataFile:
			dataFile.write('#x y\n')
			for node in self.nodes:
				ep = self.nodes[node]['endpoints']
				randomNumber = np.random.random()
				dataFile.write(' '.join([str(ep[0][0]), str(ep[0][1]), str(randomNumber) + '\n']))
				dataFile.write(' '.join([str(ep[1][0]), str(ep[1][1]), str(randomNumber) + '\n']))
				dataFile.write('\n')
			dataFile.close()

		# making gnuplot script
		with open(scriptFileName, 'w') as scriptFile:
			wr = scriptFile.write
			wr('set size square\n')
			wr('set term png size ' + str(resolution) + ',' + str(resolution) + '\n')

			# setting output name
			wr('set output \"' + fileName + '.png\"\n')

			# set plot title
			wr('set title \"' + fileName + '\" font \",250\"\n')

			# set plot top and bottom margins
			wr('set tmargin 80\n')
			wr('set bmargin 80\n')

			# setting plot bounds 
			xBounds = [-0.05 * self.width, 1.05 * self.width]
			wr('set xrange [' + str(xBounds[0]) + ':' + str(xBounds[1]) + ']\n')
			yBounds = [-0.05 * self.height, 1.05 * self.height]
			wr('set yrange [' + str(yBounds[0]) + ':' + str(yBounds[1]) + ']\n')

			# setting x and y labels
			wr('set xlabel \"Position in um\" font \",160\" offset 0,-30\n')
			wr('set ylabel \"Position in um\" font \",160\" offset -30,0\n')
			
			# setting tic labels
			wr('set xtics font \",100\" offset 0,-10\n')
			wr('set ytics font \",100\"\n')
			
			# other settings
			wr('set key off\n')
			wr('set pm3d map\n')
			wr('set palette color\n')
			wr('set cbrange [0:1] # [0:1] noreverse writeback\n')
			wr('set cbtics font \",100\"\n')

			#plot
			wr('plot \"' + dataFileName + '\" using 1:2:3 with lines palette')

			scriptFile.close()

		# running gnuplot script
		os.system('gnuplot ' + scriptFileName)

		if show:
			os.system('open ' + fileName + '.png')


	def plot(self, showJunctions = False):
		lines = []
		colors = []
		endpoints = [self.nodes[node]['endpoints'] for node in self.nodes]
		fig, ax = pl.subplots()
		ax.set_aspect('equal')
		# options for style

		# hilighting percolating wires
		colors = ['blue' if node in self.percolatingCluster else 'black' for node in self.nodes]

		# making the electrodes thicker 
		lineWidths = [4 if node in [self.topElectrode, self.bottomElectrode] else .5 for node in self.nodes]
	
		# plotting the line segments
		lc = mc.LineCollection(endpoints, colors = colors, linewidths = lineWidths)
		ax.add_collection(lc)

		# plotting the junctions
		if showJunctions:
			junctionsX = nx.get_edge_attributes(self, 'x').values()
			junctionsY = nx.get_edge_attributes(self, 'y').values()
			plt.scatter(junctionsX, junctionsY, c = 'green', s = 5)
		ax.set_xlim(0, self.width)
		ax.set_ylim(0, self.height)
		plt.xlabel('um')
		plt.ylabel('um')
		plt.show()	



	def dumb_plot(self, highlightList, centralNodes):
		lines = []
		colors = []
		endpoints = [self.nodes[node]['endpoints'] for node in self.nodes]
		fig, ax = pl.subplots()
		ax.set_aspect('equal')
		# options for style

		# hilighting percolating wires
#		colors = ['red' if node in highlightList else 'black' for node in self.nodes]
		colors = {node : 'red' if node in highlightList else 'black' for node in self.nodes}
		for node in centralNodes:
			colors[node] = 'blue'
		colors = list(colors.values())

		# making the electrodes thicker 
		lineWidths = [2 if node in highlightList + centralNodes  else .5 for node in self.nodes]
	
		# plotting the line segments
		lc = mc.LineCollection(endpoints, colors = colors, linewidths = lineWidths)
		ax.add_collection(lc)

		# plotting the junctions
		#junctionsX = nx.get_edge_attributes(self, 'x').values()
		#junctionsY = nx.get_edge_attributes(self, 'y').values()
		#plt.scatter(junctionsX, junctionsY, c = 'green', s = 1)
		ax.set_xlim(0, self.width)
		ax.set_ylim(0, self.height)
		plt.xlabel('um')
		plt.ylabel('um')
		plt.show()	


# File writing functions
	def to_pickle(self, outFileName = 'netpickle.p'):
		networkProperties = {key : getattr(self, key) for key in dir(self) if (key not in dir(nx.Graph()) and not callable(getattr(self,key)) )}
		nodes = {node : self.nodes[node] for node in self.nodes}
		edges = {edge : self.edges[edge] for edge in self.edges}
		outputDict = {'networkProperties' : networkProperties, 
			'nodes' : nodes,
			'edges' : edges,
			'percolatingCluster' : self.percolatingCluster}
		pickle.dump(outputDict, open(outFileName, 'wb'))

	# saves current state of network to csv
	def to_csv(self, folderName = 'output'):
		try:
			os.mkdir(folderName)
		except OSError:
			print('Destination folder ' + folderName + ' already exists. Moving to folder and overwriting contents.')
		else:
			print('successfully made folder: ' + folderName)

		networkFileName = '/'.join([folderName, 'networkData.csv'])

		dataTypes = ['junction', 'wire']
		# this loops is so we don't have to write the same code twice if I cnahge how the attributes are stored
		# for nodes/edges
		for dType in dataTypes:
			fileName = folderName + '/' + dType + 'Data.csv'
			with open(fileName, mode = 'w') as file:
				writer = csv.writer(file, delimiter = ',')
				header = ['wire1','wire2'] if dType == 'junction' else ['wire']
				if dType == 'junction':
					data = self.edges(data = True)
				else:
					data = self.nodes(data = True)

				attributesIndex = 2 if dType == 'junction' else 1
				# kludge to get dict of attribute keys for each graph object. assumes all graph objects
				# have same number of attributes
				for item in data:
					keys = item[attributesIndex].keys()
					break

				# make row of headers
				for key in keys:
					header.append(key)
				writer.writerow(header)

				# write rows of data
				attributeKeys = []
				for item in data:
					row = list(item[0:attributesIndex])
					theseKeys = []
					for key in item[attributesIndex].keys():
						row.append(item[attributesIndex][key])
						theseKeys.append(key)
					writer.writerow(row)
					# we are keeping track of the order in which the keys are written in every row. This will allow us to correct data
					# that may be generated with inconsistent columns without totally regenerating the data.
					attributeKeys.append(theseKeys)

			attributeKeysMatch = all([ x == attributeKeys[0] for x in attributeKeys])
			# if they don't all match, we store what they are and throw an error
			if not attributeKeysMatch:
				print('Error: ' + dType + ' column conventions are inconsistent.')
				attributeKeysFileName = folderName + '/ERROR_' + dType + 'Keys.csv'
				with open(attributeKeysFileName, mode = 'w') as file:
					writer = csv.writer(file, delimiter = ',')
					for row in attributeKeys:
						writer.writerow(row)

		# writing the total network data to file. Note that we may not know the sheet resistance yet so we have to account for that.
		with open(networkFileName, mode = 'w') as file:
			writer = csv.writer(file, delimiter = ',')
			header = ['width' , 'height', 'nwLength', 'percMultiple' , 'buffer', 'topElectrode', 'bottomElectrode']
			data = [self.width, self.height, self.nwLength, self.percMultiple, self.buffer, self.topElectrode, self.bottomElectrode]
			try:
				self.sheetResistance
				header.append('sheetResistance')
				data.append(self.sheetResistance)
			except AttributeError:
				print('Sheet Resistance currently unknown, so no sheet resistance data saved')
			writer.writerow(header)
			writer.writerow(data)

	# finds node whose center is closest to the given point
	def find_closest_node(self, x,y):
		data = [(node, self.nodes[node]['x'], self.nodes[node]['y']) for node in self.nodes if node not in [self.topElectrode, self.bottomElectrode]]
		data = np.array(data)
		pos = data[:,1:]
		point = np.array([(x,y)])
		displacements = pos - point
		distances = np.linalg.norm(displacements, axis = 1)
		minDistanceIndex = np.argmin(distances)
		return data[minDistanceIndex, 0]


	def to_netlist(self,voltage = 1, netlistName = 'netlist'):
		if type(voltage) == str:
			raise TypeError("Voltage must be numeric. The first argument, if not using kwargs, is the voltage. second argument is netlistName")
		xw = XyceWriter(self)
		xw.make_header()
		xw.make_source(self.topElectrode, self.bottomElectrode, 'VIN', voltage)
		xw.make_topology()
		xw.make_analysis_statement()
		# choosing what to print
		currentsList = ['VIN']
		percolatingSet = {node for node in self.percolatingCluster}
		voltagesList = [str(int(node)) for node in percolatingSet]
		powersList = [edge[2]['name'] for edge in self.edges(data = True) if edge[0] in percolatingSet and edge[2]['resistance'] != 0]
		xw.make_print_statement(powersList=powersList, voltagesList=voltagesList, currentsList = currentsList)
		xw.write_netlist(outFileName = netlistName)	
	
	def to_gnp_data(self, 
		outFile = 'plot_data',
		z = None,
		showElectrodes = False,
		showJunctions = False):
		# z is the third column of data that *can* be used to color code lines in the network
		# if z is none, no third column is used
		# if z is a dict, it is treated as though the keys are nodes and their
		# values are the z values
		# if z is a string, it's treated as though that string is the attribute key of node data stored in graph

		# showJunctions can be three things: True, False, or a List of which junctions to show if don't want to
		# show all of them
		# removing electrodes
		nodeList = set(self.nodes)
		if not showElectrodes:
			nodeList = nodeList - {self.topElectrode, self.bottomElectrode}
		nodeList = list(nodeList)

		# getting the z data
		dataDict = {}
		if z is not None:
			for node in nodeList:
				if type(z) == dict:
					# making sure that if not all nodes are in our z dictionary this does not pass silently
					if node in z:
						dataDict[node] = z[node]
					else:
						print('Node', node, 'not found in dictionary of z values. Skipping this node.')
				elif type(z) == str:
					# if z is a string, it is treated as a key for the node attribute with that name
					# if a node does not have that attribute, a warning message is shown
					try:
						dataDict[node] = self.nodes[node][z]
					except KeyError:
						print('Node', node, 'does not have attribute named \"', z + '\". Skipping this node')
			# making the nparray of data with nodes sorted by their z attribute
			data = np.empty((len(dataDict),2))
			nodeKeys = np.array(list(dataDict.keys()))
			zValues = np.array(list(dataDict.values()))
			data[:,0] = nodeKeys
			data[:,1] = zValues
			sortedIndices = data[:,1].argsort()
			data = data[sortedIndices]
			nodeList = list(data[:,0])

		# making data file for gnuplot to read
		with open(outFile, 'w') as f:
		#	nodeSet = set(self.nodes)
		#	# removing the electrodes from the list of wires
		#	if not includeElectrodes:
		#		nodeSet = nodeSet - {self.topElectrode, self.bottomElectrode}
			
			for node in nodeList:
				# the data must be writte in the form
				# x1 y1 z1
				# x2 y2 z2
				#
				# x1' y1' z1'
				# .....
				# getting endpoints and converting to string
				skipNode = False
				ep1, ep2 = self.nodes[node]['endpoints']
				ep1 = (str(ep1[0]), str(ep1[1]))
				ep2 = (str(ep2[0]), str(ep2[1]))
				line1 = ' '.join(ep1)
				line2 = ' '.join(ep2)

				if type(z) == dict:
					# making sure that if not all nodes are in our z dictionary this does
					# not pass silently
#					if node in z:
					line1 = ' '.join([line1, str(z[node]) ] )
					line2 = ' '.join([line2, str(z[node]) ] )
#					else:
#						print('Node', node, 'not found in dictionary of z values. Skipping this node')
#						skipNode = True
				elif type(z) == str:
					# if z is a string, it is treated as a key for the node attribute dictionaries
					# if a node doesn't have this attribute a warning message is shown
				#	if z in self.nodes[node]:
					line1 = ' '.join([line1, str(self.nodes[node][z]) ])
					line2 = ' '.join([line2, str(self.nodes[node][z]) ])
				#	else:
						# if the node doesn't have z attribute
				#		print('Node', node, 'has no attribute\"', z + '\"', '- Skipping this node')
				#		skipNode = True
				#if not skipNode:
				f.write(line1 + '\n')
				f.write(line2 + '\n')
				f.write('\n')
		
		# making a separate data file with junctions
		if showJunctions or type(showJunctions) == list:
			with open('_'.join([outFile, 'junctions']), 'w') as f:
				f.write('#junctionX junctionY\n')
				edgesToUse = self.edges if showJunctions == True else showJunctions 
				for edge in edgesToUse:
					x = self.edges[edge]['x']
					y = self.edges[edge]['y']
					line = ' '.join([str(x), str(y)])
					f.write(line + '\n')
			

	# show percolating wires
	def show_percolating(self, output = 'nw_img.png', openImage = False):
		print('Gnuplotlib was deemed an expensive dependency and was removed. Need to rewrite this function using matploblib collections')
		# self.find_percolating_cluster()
		# electrodes = {self.topElectrode, self.bottomElectrode}
		# wires = set(self.nodes) - electrodes

		# # container for all curves
		# curves = []

		# # making the wire segments
		# for node in wires:
		# 	x = np.array([self.nodes[node]['endpoints'][n][0] for n in range(2)])
		# 	y = np.array([self.nodes[node]['endpoints'][n][1] for n in range(2)])
		# 	if node in self.percolatingCluster:
		# 		curveOptions = {'with' : 'lines lc \"blue\" lw 6'}
		# 	else:
		# 		curveOptions = {'with' : 'lines lc \"black\" lw 6'}
		# 	curves.append( (x, y, curveOptions) )

		# # making electrodes
		# electrodeThickness = self.nwLength
		# for node in electrodes:
		# 	x = np.array([self.nodes[node]['endpoints'][n][0] for n in range(2)])
		# 	y = np.array([self.nodes[node]['endpoints'][n][1] for n in range(2)])
		# 	# note the alpha value (first two digits after #) is the transparency
		# 	# higher numbers mean more transparent
		# 	curveOptions = {'with' : 'filledcurves lc rgb \"#60454545\"',
		# 				'tuplesize' : 3}
		# 	yLowerBound = y if node == self.topElectrode else y - electrodeThickness
		# 	yUpperBound = y + electrodeThickness if node == self.topElectrode else y
		# 	curves.append( (x, yLowerBound, yUpperBound, curveOptions) )
			
		# #plot options
		# plotOptions = {'xrange' : ''.join([str(-1.1 * self.nwLength),
		# 					':',
		# 					str(self.width + 1.1 * self.nwLength)]),
		# 		'yrange' : ''.join([str(-1.1 * self.nwLength),
		# 					':',
		# 					str(self.height + 1.1 * self.nwLength)]),
		# 		'terminal' : 'pngcairo size 4096,4096',
		# 		'hardcopy' : output if output[-4:] == '.png' else output + '.png',
		# 		'set' : ['lmargin 0',
		# 				'tmargin 0',
		# 				'rmargin 0',
		# 				'bmargin 0'],
		# 		'unset' : ['grid',
		# 				'xtics',
		# 				'ytics']
		# 		}
							

		# gp.plot(*curves, **plotOptions)
		# if openImage:
		# 	os.system('open ' + output)
		

	# creating the image of the network from the gnuplot data
	def to_img(self, 
		outFile = 'plot',
		z = None,
		title = 'Plot', 
		xLabel = 'X Axis',
		yLabel = 'Y Axis',
		zLabel = 'Z Value',
		lineWidth = 7,
		pointSize = 3,
		showJunctions = False,
		showElectrodes = True,
		extraCode = [],
		cleanup = False,
		openImage = True):

		# make intermediate file names
		dataName = outFile + '_data'
		scriptName = outFile + '_script'


		# make the source data
		self.to_gnp_data(outFile = dataName, z = z, showJunctions = showJunctions, showElectrodes = showElectrodes)
		
		# make the script text
		scriptText = ['set size square',
				'set term png size 8192, 8192 linewidth ' +str(lineWidth),
				'set output \"' + outFile + '.png\"',
				'set title \"' + title + '\" font \",250\"',
				'set tmargin 80',
				'set bmargin 80',
				'set xrange [' + str(-self.buffer * self.nwLength) + ':' + str(self.width + self.buffer*self.nwLength) + ']',
				'set yrange [' + str(-self.buffer * self.nwLength) + ':' + str(self.height + self.buffer*self.nwLength) + ']',
				'set xlabel \"' + xLabel + '\" font \",160\" offset 0,-30',
				'set ylabel \"' + yLabel + '\" font \",160\" offset -30,0',
				'set xtics font \",100\" offset 0,-10',
				'set ytics font \",100\"',
				'set key off']
		# only adding these lines to script text if there exists a third column of data to use as color
		# plotCommand sets the appropriate plot command for how many columns of data there are
		if z != None:
			scriptText = scriptText + ['stats \"' + dataName + '\" using 3',
				'p0 = STATS_min',
				'p1 = STATS_lo_quartile',
				'p2 = STATS_median',
				'p3 = STATS_up_quartile',
				'p4 = STATS_max',
				'set palette defined (p0 \"gray20\", p1 \"midnight-blue\", p2 \"blue\", p3 \"magenta\", p4 \"red\")',
				'set cbrange [STATS_min:STATS_max]',
				'set cblabel \"' + zLabel + '\" font \",160\" offset 70,0', 
				'set cbtics font \",100\"']
			plotCommand = 'plot \"' + dataName + '\" using 1:2:3 with lines palette'
		else:
			plotCommand = 'plot \"' + dataName + '\" using 1:2 with lines lc \"black\"'
			
		# adding command to show junctions
		if showJunctions:
			plotCommand +=  ', \"' + '_'.join([dataName, 'junctions']) + '\" with points pointtype 7 pointsize ' + str(pointSize) + ' lc \"black\"'


		# adding any extra code
		extraCode = extraCode if type(extraCode) == list else [extraCode]
		scriptText = scriptText + extraCode

		# adding the plot command into the script text
		scriptText = scriptText + [plotCommand]

		# adding new line characters
		scriptText = '\n'.join(scriptText)

		# write script to file
		with open(scriptName, 'w') as f:
			f.writelines(scriptText)

		# run script to make plot
		try:
			command = ['gnuplot', scriptName]
			subprocess.check_call(command, stdout = subprocess.DEVNULL)
			if openImage:
				openCommand = ['open', outFile + '.png']
				subprocess.check_call(openCommand, stdout = subprocess.DEVNULL)
		except subprocess.CalledProcessError as err:
			print('Error in plotting:')
			print('Error message:', err.output)

		if cleanup:
			os.remove(scriptName)
			os.remove(dataName)
			
		
	def solve_circuit(self, voltage = 1, rShunt = 10**(12), use = None):
		# we attempt to solve the circuit using the matrix equation
		# A^T K A x = A^T K v
		# A is a matrix denoting the drop in voltage across a resistor
			# in terms of node potentials
		# K is a matrix with 1/R_i on the main diagonal
		# v is a column vector which has value V if 
		# the row corresponds to a resistor touching the top electrode
		# and is 0 otherwise
		def _matrix_operations(A, K, v, use):
			if use == 'matmul':
				ATK = np.matmul( np.transpose(A), K)
				ATKA = np.matmul( ATK, A)
				ATKv = np.matmul( ATK, v)
			elif use == 'dot':
				ATK = A.T.dot(K)
				ATKA = ATK.dot(A)
				ATKv = ATK.dot(v)
			elif use == 'blas':
				ATK = scipy.linalg.blas.sgemm(1, A.T, K)
				ATKA = scipy.linalg.blas.sgemm(1, ATK, A)
				ATKv = scipy.linalg.blas.sgemv(1, ATK, v)
			elif use == 'numexpr':
				ATK = ne.evaluate("np.matmul( np.transpose(A), K)")
				ATKA = ne.evaluate("np.matmul( ATK, A)")
				ATKv = nne.evaluate("np.matmul( ATK, v)")
			
			return ATKA, ATKv

			
		self.timer.start('solve_circuit')


		# first we add shunt resistors to any nodes of degree 1
		# these connect dangling ends to ground
		# we call these dangling end
		danglingEnds = [node for node in self if self.degree(node) == 1]
		addedEdges = []
		for node in danglingEnds:
			self.add_edge(node, self.bottomElectrode,
					resistance = rShunt,
					resistanceType = 'shunt')
			addedEdges.append((node, self.bottomElectrode))

		# making  matrix A 
		# rows = number of resistors
		# columns = number of nodes - 2 (don't count top and bottom electrodes
		# b/c their voltages are already known)
		edgesList = list(self.edges)
		nodesWithoutElectrodes = set(self.nodes) - {self.topElectrode, self.bottomElectrode}
		nodesWithoutElectrodes = list(nodesWithoutElectrodes)
		nodeIndices = {node : n for n,node in enumerate(nodesWithoutElectrodes)}
		A = np.zeros( (len(edgesList), len(nodesWithoutElectrodes)))
		for row, edge in enumerate(edgesList):
			if self.topElectrode in edge and self.bottomElectrode in edge:
				raise Exception('Top and Bottom electrodes are directly connected. Problem.')
			elif self.topElectrode in edge and self.bottomElectrode not in edge:
				# put a +1 in the correct place for non-electrode node
				# get index of non electrode node
				nonElectrodeNode = edge[ edge.index(self.topElectrode) - 1]
#				col = nodesWithoutElectrodes.index(nonElectrodeNode)
				col = nodeIndices[nonElectrodeNode]
				A[row, col] = 1
			elif self.topElectrode not in edge and self.bottomElectrode in edge:
				# put a -1 in the correct place for non-electrode node
				# get index of non electrode node
				nonElectrodeNode = edge[ edge.index(self.bottomElectrode) - 1]
#				col = nodesWithoutElectrodes.index(nonElectrodeNode)
				col = nodeIndices[nonElectrodeNode]
				A[row, col] = -1
			else:
				# this block is entered when neither electrode is 
				# part of the edge
				for n,node in enumerate(edge):
					sign = (-1)**n
#					col = nodesWithoutElectrodes.index(node)
					col = nodeIndices[node]
					A[row, col] = sign
		A = A.astype(dtype = np.float32)

		# Making vector v
		v = np.zeros( (len(edgesList), 1))
		for row, edge in enumerate(edgesList):
			if self.topElectrode in edge:
				v[row, 0] = voltage
		v = v.astype(dtype = np.float32)

		# Making matrix K
		# has 1/R_i on main diagonal.
		# all other elements zero
		self.timer.start('matrix operations')
		diagElements = [1/self.edges[edge]['resistance'] for edge in edgesList]
		K = np.diag(diagElements)
		K = K.astype(dtype = np.float32)
		
		# setting up the matrix equation
		# A^T K A x = A^T K v
#		ATK = np.matmul( np.transpose(A), K)
#		ATKA = np.matmul( ATK, A)
#
#		ATKv = np.matmul( ATK, v)
		ATKA, ATKv = _matrix_operations(A, K, v, use = use)
		self.timer.stop('matrix operations')
		self.timer.stop('solve_circuit')


		return ATKA, ATKv

# Misc Functions			
	# Utilify function: This uses NetworkX to find the cluster connected to the top electrode. It also
	# determines whether the network is percolating (by checking if the aforementioned cluster contains 
	# both top and bottom electrodes) and writes this to self.isPercolating.
	def find_percolating_cluster(self, verbose = True):
		if verbose:
			print('Finding percolating cluster.')
		T = nx.dfs_tree(self, self.topElectrode)
		# self.percolatingCluster = np.array(list(T.nodes)).astype(int)
		self.percolatingCluster = list(T.nodes)
		# the percolating cluster starts from the top electrode, so the top electrode is in it by definition. Therefore only necessary to test if bottom electrode is in it.
		self.isPercolating = self.bottomElectrode in self.percolatingCluster
		if verbose:
			print('Percolating cluster found.')

	# get actual density as a fraction of the percolation threshold
	def get_perc_multiple(self):
		# NL = total length of wires on sample excluding electrodes
		# A = sample area
		# L = length of single wire
		nodesExcludingElectrodes = set(self.nodes) - {self.topElectrode, self.bottomElectrode}
		nodeLengths = [self.nodes[node]['length'] for node in nodesExcludingElectrodes]
		NL = sum(nodeLengths)
		A = self.width * self.height
		L = self.nwLength
		return NL * L / (5.63726 * A)


	# Utility function: gets the edges and power of the edge with maximum power dissipation.
	def get_max_power(self):
		power = nx.get_edge_attributes(self, 'power')

		# I don't understand how this line works but it works
		maxKey = max(power, key = power.get)
		return {maxKey : self.get_edge_data( maxKey[0], maxKey[1])['power']}

		nodeLengths = [self.nodes[node]['length'] for node in nodesExcludingElectrodes]
		NL = sum(nodeLengths)
		A = self.width * self.height
		L = self.nwLength
		return NL * L / (5.63726 * A)

	# Utility function: gets the edges and power of the edge with maximum power dissipation.
	def get_max_power(self):
		power = nx.get_edge_attributes(self, 'power')

		# I don't understand how this line works but it works
		maxKey = max(power, key = power.get)
		return {maxKey : self.get_edge_data( maxKey[0], maxKey[1])['power']}

	# Calculate the tortuosity of a network, optional plotting ability as well. 
	def calculate_tortuosity(self, meshSize = [], plotOption = False):
		# meshSize =  float(round(np.sqrt(self.width*self.height)/(np.sqrt((5.63726*self.percMultiple)/(self.nwLength**2))*10)))
		if not meshSize: 
			meshSize =  float(round(self.percMultiple*np.sqrt(self.width*self.height)/(self.nwLength*2)))
		
		#Load the node,voltage, and position data into wireVolt
		wireVolt = []
		attributesToGet = ['voltage' , 'x' , 'y']
		for node in self.nodes():
			row = [node]
			attributeDictionary = self.nodes[node]
			for name in attributesToGet:
				row.append(attributeDictionary[name])
			wireVolt.append(row)
		wireVolt = np.asarray(wireVolt)
		#remove the extra wires outside the network boundaries 
		wireVolt = wireVolt[wireVolt[:,2] <= self.width, :]
		wireVolt = wireVolt[wireVolt[:,2] >= 0, :] 
		wireVolt = wireVolt[wireVolt[:,3] <= self.height, :] 
		wireVolt = wireVolt[wireVolt[:,3] >= 0, :]

		#only keep percolating wires
		wireVolt = wireVolt[np.isin(wireVolt[:,0], self.percolatingCluster), :]

		#Create a mesh grid and interpolation of all the voltages based on the meshSize defined earlier
		# and use it to create a contour object
		X,Y = np.meshgrid(np.linspace(0.,self.width,meshSize),np.linspace(0.,self.height,meshSize))	#create mesh grid of the x,y points
		Z = griddata((wireVolt[:,2],wireVolt[:,3]), wireVolt[:,1], (X,Y), method = 'nearest')	#interpolate to find the closest 'Z' point i.e. the voltage
		# make sure that this contour plot isn't shown until the end and only if plotOption == True
		plt.ioff()
		conPlot = plt.contour(X,Y,Z, int(0.5*meshSize))		#create the contour object
		#plt.scatter(wireVolt[:,2],wireVolt[:,3],c=wireVolt[:,1],s=2 )

		#Finds the lines that comprise the contour plot and manipulates the individual points to find a tortuosity
		contours = []
		segmentTorts = []
    	# for each contour line
		for cc in conPlot.collections:
		    paths = []
		    # for each separate section of the contour line
		    for pp in cc.get_paths():
		        xy = []
		        # for each segment of that section
		        for vv in pp.iter_segments():
		            xy.append(vv[0])
		        xy = np.asarray(xy)
		        numberPts = len(xy[:,0])
		        if numberPts > 2 and round(displacement(xy), 5) != 0: #avoiding floating point errors in the displacement
		        	# print(displacement(xy))
		        	segmentTorts.append([arc_length(xy[:,0],xy[:,1])/displacement(xy), numberPts])	#length over displacement with a number of points index as well
		# print(segmentTorts)
		segmentTorts = np.asarray(segmentTorts)		#convert to numpy array
		averageTortuosity = (np.sum(np.multiply(segmentTorts[:,0],segmentTorts[:,1])))/(np.sum(segmentTorts[:,1])) #weighted average tortuosity
		if plotOption: #Display contour plot if you so desire
			plt.show()
		return averageTortuosity

	# Relabels nodes where mapping = {old_node_label : new_node_label}
	# can relabel arbitrary numbers of nodes
	# only necessary because relabeling the graph otherwise erases all the custom properties
	def relabel(self, mapping):
		tempG = self
		# it is necessary to set copy = False in order to do a partial remapping
		nx.relabel_nodes(self, mapping, copy = False)
		for key in vars(tempG).keys():
			if key not in vars(self).keys() :
				setattr(self, key, getattr(tempG, key))

	def calc_internal_resistivity(self, diam, temp):
		# alias of module function to preserve backwards compatibility
		return calc_internal_resistivity(diam, temp)

	def calc_wire_mass(self, diam, length):
		return float(np.pi * (diam/2)**2 * length * self.silverDensity)

	# temporary fix. May not need this function at all.
	def calc_junction_mass(self, diam1, diam2, resistanceType = 'cont'):
		# approximating all wires as having square cross-sections
		# for these types of square exxtrusions, the volume enclosed by the juntion
		# when both have same 'diameter' d is
		# V_enclosed = sqrt(2) d^3
		# if cylinders have different diameters, we just use the above equation but take that average of the two diameters
		# for the value of d
		if resistanceType == 'cont':
			d = np.average((diam1, diam2))
			correctionFactor = np.pi/4 # this is the ratio of areas of a circle of diameter d to a square with sidelength d
			# we use it to correct the fact that our masses are not quite right for cylinders
			mass = correctionFactor * np.sqrt(2) * d**3 * self.silverDensity
		elif resistanceType == 'int':
			# here we claim the internal resistors are 1 diameter long
			mass  = self.calc_wire_mass(diam1, diam1)
			# note we only use one of the diameters because for an internal resistor the diameters will be equal
		return mass

	def positive_random_normal(self, mean, SD, sampleSize):
		# this generates a random normal distribution BUT whenever a negative number is drawn it redraws those numbers
		sample = np.random.normal(mean, SD, sampleSize)
		while min(sample) <= 0:
			sample[sample <= 0] = np.random.normal(mean, SD, len(sample[sample<=0]))
		return sample
	
	def gdf(self):
		df_list = []
		assumption_warning = False
		for node in self.nodes:
			row = deepcopy(self.nodes[node])
			row['key'] = node
			row['el_type'] = 'wire' if node not in [self.topElectrode, self.bottomElectrode] else 'electrode'
			row['nx_type'] = 'node'
			if 'shape' not in row:
				row['shape'] = LineString(row['endpoints'])
				assumption_warning = True
			df_list.append(row)
		for edge in self.edges:
			row = deepcopy(self.edges[edge])
			row['key'] = edge
			row['nx_type'] = 'edge'
			row['el_type'] = row['resistanceType']
			if 'shape' not in row:
				row['shape'] = Point(row['x'], row['y'])
			df_list.append(row)

		if assumption_warning:
			print('Assuming all nanowires are perfectly straight to generate shapely geometries')
		
		df = pd.DataFrame(df_list).rename(columns = {'shape' : 'geometry'})
		return gpd.GeoDataFrame(df)

		# gdf = gpd.GeoDataFrame(pd.DataFrame(df_list))
		# if 'geometry' not in gdf.columns:
		# 	line_idx = gdf.el_type.isin(['wire','electrode'])
		# 	point_idx = ~line_idx
		# 	gdf.loc[line_idx, 'geometry'] = gdf.loc[line_idx, 'endpoints'].apply(
		# 		lambda x : LineString(x)
		# 		)
		# 	gdf.loc[point_idx,'geometry'] = gdf.loc[point_idx, ['x','y']].apply(
		# 		lambda row : Point(row.x, row.y),
		# 		axis = 1
		# 	)
		# 	return gdf.set_geometry('geometry')

class XyceWriter():
	def __init__(self, network):
		self.network = network
		self.orderedSections = ['header', 'sources', 'topology', 'analysis', 'print']
		self.sections = {name : [] for name in self.orderedSections}


	def make_header(self, header = None):
		if header == None:
			networkProperties = {key : getattr(self, key) for key in dir(self) if (key not in dir(nx.Graph()) and type(getattr(self, key) != 'method'))}
			headerList = [key + ' = ' + str(val) for key,val in networkProperties.items()]
			header = ' | '.join(headerList)
		self.sections['header'].append(header)

	def make_source(self, node1, node2, name, value):
		newSource = [name, str(int(node1)), str(int(node2)), 'DC', str(value)]
		newSource = ' '.join(newSource)
		self.sections['sources'].append(newSource)

	def make_topology(self):
		# making real resistors
		percolatingSet = {node for node in self.network.percolatingCluster}
		realResistorList = [ [edge[2]['name'], str(int(edge[0])), str(int(edge[1])), str(edge[2]['resistance']) ] for edge in self.network.edges(data = True) if edge[0] in percolatingSet ]
		topology = [' '.join(resistor) for resistor in realResistorList]

		#making shunt resistors
		danglingEndsList = [node for node in percolatingSet if self.network.degree(node) == 1]
		shuntResistorList = [ ['_'.join(['rshunt', str(int(node)), str(int(self.network.bottomElectrode)) ]), str(int(node)), str(int(self.network.bottomElectrode)), '10e9'] for node in danglingEndsList]
		topology += [' '.join(resistor) for resistor in shuntResistorList]
		self.sections['topology'] = topology

	def make_analysis_statement(self, options = None):
		self.sections['analysis'].append('.OP')
		if options != None:
			self.sections['analysis'].append('.OPTIONS ' + options)
	
	def make_print_statement(self, powersList = None, voltagesList = None, currentsList = None):
		# adding sources to print statement
		printStatement = ['.PRINT DC FORMAT=CSV']
		if powersList != None:
			printStatement += ['P(' + element + ')' for element in powersList]
		if voltagesList != None:
			printStatement += ['V(' + element + ')' for element in voltagesList]
		if currentsList != None:
			printStatement += ['I(' + element + ')' for element in currentsList]
		self.sections['print'] = [' '.join(printStatement)] # these brackets necessary because the write_netlist function expects each section to be a list

	def write_netlist(self, outFileName = 'netlist'):
		lines = []
		for sectionName in self.orderedSections:
			lines += self.sections[sectionName]
		lines += ['.END']
		with open(outFileName, "w") as file:
			file.writelines("%s\n" % line for line in lines)
	
	

	# def to_netlist(self, netlistName = 'netlist', voltage = 1):
	# 	#Note that the output file from xyce/sice will be netlistName.prn and it will be a space separated file
	# 	# print(' '.join(['Opening file with name:', netlistName]))
	# 	# edges = self.edges(data = True)
	# 	# percolatingCluster = self.percolatingCluster

	# 	#title line (not read by spice/xyce)
	# 	# header = 'width = ' + str(self.width) + ' | height = ' + str(self.height) + ' | nwLength = ' + str(self.nwLength)
	# 	# header = header + ' | percMultiple = ' + str(self.percMultiple) + ' | buffer = ' + str(self.buffer)

	# 	# initialize text list object where each item in the list will eventually become a line in the netlist
	# 	# textList  = [header]

	# 	#voltage source
	# 	# textList.append('VIN ' + str( int(self.topElectrode) ) + ' ' + str( int(self.bottomElectrode) ) + ' DC ' + str(voltage))	

	# 	# real resistors	
	# 	print('Writing real resistors into netlist.')
	# 	# this is the solution below. we use numpy to figure out which resistors are percolating.
	# 	# isPercolatingArray = np.isin(self.edges, percolatingCluster)
	# 	# realResistorList = [' '.join([edge[2]['name'], str(int(edge[0])), str(int(edge[1])), str(edge[2]['resistance'])]) for counter, edge in enumerate(edges) if isPercolatingArray[counter][0]]
	# 	# textList = textList + realResistorList

	# 	# Adding shunt resistors which connect dangling nodes to ground
	# 	# print('Adding shunt resistances to dangling ends.')
	# 	# shuntResistorList = [' '.join([''.join(['rshunt', str(int(node))]), str(int(node)), '0 10e9']) for node in percolatingCluster if self.degree(node) == 1 ]
	# 	# textList = textList + shuntResistorList


	# 	# Tell spice what kind of analysis to perform (operating point)
	# 	textList.append('.OP')

	# 	# Tell spice what data to print and in what format.
	# 	print('Writing print statement.')
	# 	# here we store the print statement as a list of strings with no spaces between them. The spaces will be added during writing.
	# 	printStatement = ['.PRINT DC FORMAT=CSV']

	# 	#tell spice to print current through entire circuit
	# 	printStatement.append('I(VIN)')

	# 	#tell spice to print voltages at all nodes
	# 	print('Writing voltages in print statement.')
	# 	voltagesList = [''.join(['V(', str(int(node)), ')']) for node in percolatingCluster]
	# 	printStatement = printStatement + voltagesList

	# 	# Tell spice to print the power dissipated in all (non shunt) resistors. 
	# 	print('Writing power dissipation in print statement.')
	# 	powersList = [''.join(['P(', edge[2]['name'], ')']) for counter, edge in enumerate(edges) if isPercolatingArray[counter][0] ]
	# 	printStatement = printStatement + powersList
	# 	printStatement = ' '.join(printStatement)
	# 	textList.append(printStatement)
	# 	textList.append('.END')

	# 	with open(netlistName, "w") as file:
	# 		file.writelines("%s\n" % line for line in textList)
	# 	print('Netlist complete.')		


		

##########
# All functions below this point are *not* methods of the NanowireMesh class, rather they are 
# functions defined in the NanowireMesh module file. Inside this file, the below functions can be 
# called by name. Outside this file, these functions must bed called by nwm.function_name() 
# (assuming you have imported NanowireMesh as nwm).

# In general, routines that may be useful even without a Nanowire_Mesh object should be included as functions
# here. An example of this is a function that reads output from spice simulations. There may be occasions
# where we want to read output from spice without instantiating a Nanowire_Mesh object.

def copy(g, as_view = False):
	"""Essentially the same as networkx graph.copy() method but also copies graph attributes"""
	h = g.copy(as_view)
	for attr in g.propsToWrite:
		obj = getattr(g, attr)
		if not callable(obj):
			setattr(h, attr, deepcopy(obj))
	return h


# Determines whether a string contains purely numerical information
def isNumber(string):
	try:
		float(string)
		return True
	except ValueError:
		return False

def displacement(data): 		#returns scalar displacement from the contour lines
	a, b = [data[0,:], data[-1,:]] 
	a_min_b = b - a 
	return np.sqrt(np.einsum('i,i',a_min_b,a_min_b)) 

def arc_length(x, y): #arc length used for contour lines
	npts = len(x)
	arc = 0
	for k in range(0, npts-1):
		arc = arc + np.sqrt((x[k+1] - x[k])**2 + (y[k+1] - y[k])**2)
	return arc

# the below is a function to be used with a solver. This version of the function has no syntactic sugar. 
# There is another version that I'm trying to use that is much more readable. 
# def makeResistanceFunctionToSolve(D, Rw, Rs):
# 	na = 0.20276 * np.pi * D / 2
# 	return lambda Rc : Rs - Rw / D * (1/2 * rm - np.sqrt(Rc * rm / (2 * Rw * na)) * np.tanh(np.sqrt(Rw * na * rm / (2 * Rc))) )**(-1)

def forro_sheet_resistance(percMultiple, Rw, Rc):
	D = 5.63726 * percMultiple
	na = 0.20276 * np.pi * D / 2
	rm = (na - 1 + Rc * (Rc + Rw / (na + 1) )**(-1) )/(na + 1)
	Rs = Rw / D * (1/2 * rm - np.sqrt(Rc * rm / (2 * Rw * na)) * np.tanh(np.sqrt(Rw * na * rm / (2 * Rc) ) ) )**(-1)
	return Rs

def makeResistanceFunctionToSolve(percMultiple, Rw, Rs):
	return lambda Rc : Rs - rSheet(percMultiple, Rw, Rc)

def solveForContactResistance(percMultiple, Rw, Rs):
	F = lambda Rc : Rs - rSheet(percMultiple, Rw, Rc)
	return optimize.broyden1(F, 40)

# reads the csv file that xyce spits out
# note: top electrode must be specified for sheet resistance to be calculated
def read_xyce_output(inFile, disableTQDM = False, topElectrode = None):
	print('Reading Xyce output.')
	with open(inFile) as file:
		table = csv.reader(file, delimiter = ',')
		
		# making sure that the memory limit accommodates the table
		# whil also preventing overflow errors
		maxSize = sys.maxsize
		while True:
			# try to set the csv field size to max int.
			# if it fails, keep decreasing field size by order of magnitude
			# until it doesn't fail
			try:
				# if this works, exit the while loop
				csv.field_size_limit(maxSize)
				break
			except OverflowError:
				#otherwise reduce maxSize by order of magnitude
				maxSize = int(maxSize / 10)
		data = list(table)

	# converting numeric strings to floats
	data[1] = [float(x) for x in data[1]]
	# convert data to a dictionary with two elements: nodes and edges. 
	# The edges entry will be a list with the following structure

	# establishing the name patterns to correspond to each type of data stored in spice output
	voltagePattern = re.compile(r'V\(\d+\)')
	powerPattern = re.compile(r'P\(R.+_\d+_\d+\)')
	currentPattern = re.compile(r'I\(VIN\)') # note this is current through entire circuit
	#note that in the below the parentheses are not escaped because they are included in the string before the regex is compiled
	topElectrodeVoltageName = 'V(' + str(int(topElectrode)) + ')' if topElectrode is not None else 'No_Top_Electrode_Specified'
	if topElectrode is None:
		print('Top electrode is not specificed so sheet resistance cannot be calculated')
	else:
		print('Top electrode voltage is named', topElectrodeVoltageName)

	numberPattern = re.compile(r'\d+') # generic pattern for matching integers.


	# The structure of output: output is a dictionary with three entries. 
	#The key 'nodes' points to a dict where each entry has the node number as a key and a dictionary of a property name  and value in column 2.
	# The key 'edges' points to a list where each row has a list of the two vertices in column 1 and a dictionary of a property name and value 
	# in column 2.
	# The key 'sheetResistance' points to the resistance.
	output = {'nodes' : {},
		'edges' : {},
		'sheetResistance': 'blank'}

	# the structure of the xyce out put is two big ass rows.
	# first row is the name of the data and second row is the value for that data.
	for col in tqdm(range(len(data[0])), disable = disableTQDM):
		name = data[0][col]
		value = data[1][col]

		#note: the logic below assumes we have a graph and not a multigraph. if nodes can touch more than once, the update
		# function will not work correctly for a dictionary.
		if voltagePattern.match(name) and name != topElectrodeVoltageName:
			nodeNumber = numberPattern.search(name).group(0)
			nodeNumber = int(nodeNumber)
			output['nodes'].update({nodeNumber : {'voltage' : value}})
		elif powerPattern.match(name):
			nodeNumbers = numberPattern.findall(name)
			nodeNumbers = [int(x) for x in nodeNumbers]
			nodeNumbers = tuple(nodeNumbers)
			output['edges'].update({nodeNumbers : {'power' : value}})
		elif currentPattern.match(name):
			print('Found source current')
			currentThroughSource = value
		elif name == topElectrodeVoltageName:
			# not sure why this is different from the first block but it's here...
			print('Found source voltage')
			voltageAcrossSource = value
			nodeNumber = numberPattern.search(name).group(0)
			nodeNumber = int(nodeNumber)
			output['nodes'].update({nodeNumber : {'voltage' : value}})
			
	if topElectrode is not None:
		if currentThroughSource == 0:
			sheetResistance = np.inf
		else:
			sheetResistance = abs(voltageAcrossSource / currentThroughSource)
		output['sheetResistance'] = sheetResistance
	return output

# updates g with properties read from output of xyce simulations
def update_with_xyce_output(g, inFile, disableTQDM = False):
	data = read_xyce_output(inFile, disableTQDM = disableTQDM)
	g.sheetResistance = data['sheetResistance']
	print('Writing node data')
	nx.set_node_attributes(g, data['nodes'])
	print('Writing edge data')
	nx.set_edge_attributes(g, data['edges'])
	print('Network update complete.')

# make network into directed graph where edges point down potential
def to_directed(g, potential = 'voltage', weight = 'resistance'):
	# make directed graph from G using the potential attribute
	gDirected = nx.DiGraph()
	# make nodes with potentials
	try:
		nodesWithPotentials = [(node , {potential : g.nodes[node][potential]}) for node in g.nodes]
		gDirected.add_nodes_from(nodesWithPotentials)
	except KeyError:
		# throw error if not all nodes have voltages
		print('Error when converting the graph to a directed graph. Not all nodes have voltages assigned')
		nodesWithoutPotentials = [node for node in g.nodes if potential not in g.nodes[node].keys()]
		print('First 10 Nodes lacking a potential named', potential, 'are:', nodesWithoutPotentials)
		raise Exception('Encountered Fatal Error. Aborting run.')
	#make list of directed edges. First node in edge is source and second node is target. Listshould include higherp otential node first and the lower potential node second
	for edge in g.edges:
		node1 = edge[0]
		node2 = edge[1]
		potential1 = g.nodes[node1][potential]
		potential2 = g.nodes[node2][potential]
		if potential1 > potential2:
			source = node1
			target = node2
		elif potential1 < potential2:
			source = node2
			target = node1
		else:
			source = target = None
			#Do not add an edge if the source and target are at the same potential
		if source != None and target != None:
			# assigning the edge weights to the new directed graph
			gDirected.add_edge(source, target)
			gDirected.edges[source, target][weight] = g.edges[source, target][weight]
	# adding potentials to the new directed graph
	for node in gDirected.nodes:
		gDirected.nodes[node][potential] = g.nodes[node][potential]
	return gDirected
#	except KeyError:
#		print('Node found without attribute:', potential)
#		print('Aborting opretation')

def electrode_centrality(g, potential = 'voltage', weight = 'resistance', normalized = True):
	# making the directed graph version
	gDirected = to_directed(g, potential = potential, weight = weight)
	if gDirected is not None:
		sources = list(gDirected.neighbors(g.topElectrode))
		targets = list(gDirected.predecessors(g.bottomElectrode))
		electrodeCentrality = nx.betweenness_centrality_subset(gDirected, 
									sources = sources, 
									targets =  targets, 
									normalized = normalized, 
									weight = weight)
	else:
		electrodeCentrality = None
	return electrodeCentrality

def edge_electrode_centrality(g, potential = 'voltage', weight = 'resistance', normalized = True):
	# making the directed graph version
	gDirected = to_directed(g, potential = potential, weight = weight)
	if gDirected is not None:
		sources = gDirected.neighbors(g.topElectrode)
		targets = gDirected.predecessors(g.bottomElectrode)
		edgeElectrodeCentrality = nx.edge_betweenness_centrality_subset(gDirected, 
									sources = sources, 
									targets =  targets, 
									normalized = normalized, 
									weight = weight)
	else:
		edgeElectrodeCentrality = None
	return edgeElectrodeCentrality


def percolation_centrality(g, attribute = 'voltage', weight = 'resistance',  downPotentialOnly = False):
	# making a directed version of our graph if necessary
	gTemp = to_directed(g, potential = attribute, weight = weight) if downPotentialOnly  else g

	# testing to make sure that all nodes have potentials assigned to them
	try:
		potentials = [g.nodes[node][attribute] for node in g.nodes]
	except KeyError:
		# throw error if not all nodes have voltages
		print('Error: Not all nodes have voltages assigned')
		nodesWithoutPotentials = [node for node in g.nodes if attribute not in g.nodes[node].keys()]
		print('First 10 Nodes lacking a potential named', attribute, 'are:', nodesWithoutPotentials)
		raise Exception('Encountered Fatal Error. Aborting run.')

	return nx.percolation_centrality(gTemp, attribute = attribute, weight = weight)


def _series_approximation(func, nStart, tolerance):
	"""
	Approximates the value of a (potentially partial or infinite) series.
	The terms in the series are defined by func[n]
	If we define S_n to be the partial sum of func[nStart] + ... func[n], then we continue summing
	until (func[n+1] - S_n) / S_n < tolerance

	args
		func = function to be evaluated at different integers and summed
		nStart = the integer at which to start the summation
		tolerance = desired accuracy
	returns
		sum of func[nStart], func[nStart + 1] + ... to within desired tolerance
	"""

	# initializing the list of terms with the first term
	n = nStart
	terms = [func(n)]
	relativeChange = 1
	while relativeChange > tolerance or n < 50:
		n += 1
		newTerm = func(n)
		relativeChange = ( np.abs(newTerm) - np.abs(sum(terms)) ) / np.abs(sum(terms))
		terms.append(newTerm)
	return sum(terms)

	
		
def extinction_coefficient(radius, mr = 0.055 + 3.32j, wavelength = 550E-9):
	"""
	Calculates wavelength-depdnecnt extinction cross_section PER METER NANOWIRE LENGTH 
	args
		radius = radius of nanowire in meters
		mr = ratio of refractive index of nanowire material to refractive index of surroundings at wavelength.
			we can use the refractive index of the nanowire material with high precision if the surroundings are air
		wavelength = wavelength
	returns
		C_ext = extinction cross section per meter of nanowire length
	"""
	# calculating derived quantities
	k = 2 * np.pi / wavelength
	xr = k * radius
	# a little syntactic sugar
	def Jn(n, z): # Bessel function of first kind
		return jv(n, z)
	
	def Jnp(n, z): # First deriv of Bessel function of first kind
		return jvp(n, z)
	
	def Hn(n, z): # Hankel function of first kind
		return hankel1(n, z)
	
	def Hnp(n, z): # First deriv of Hankel function of first kind
		return h1vp(n, z)
	
	def extinction_coefficient(wavelength, nwLength, mr, xr):
		# a terms are for C_s
		# b terms are for C_p
		real = np.real
	
		def _a_n(n, mr = mr, xr = xr): # this yields Real{a_n + b_n}
			aNumerator =  Jn(n, xr) * Jnp(n, mr * xr) - mr * Jnp(n, xr) * Jn(n, mr * xr) 
			aDenominator = Hn(n, xr) * Jnp(n, mr * xr) - mr * Hnp(n, xr) * Jn(n, mr * xr)
			a_n = aNumerator / aDenominator
			return real(a_n)
	
		def _b_n(n, mr = mr, xr = xr): # this yields Real{b_n}
			bNumerator =  mr * Jn(n, xr) * Jnp(n, mr * xr) - Jnp(n, xr) * Jn(n, mr * xr) 
			bDenominator = mr * Hn(n, xr) * Jnp(n, mr * xr) -  Hnp(n, xr) * Jn(n, mr * xr)
			b_n = bNumerator / bDenominator
			return real(b_n)
	
		def _nthTerm(n, mr = mr, xr = xr): # this yields Real{a_n + b_n}
	#		aNumerator =  Jn(n, xr) * Jnp(n, mr * xr) - mr * Jnp(n, xr) * Jn(n, mr * xr) 
	#		aDenominator = Hn(n, xr) * Jnp(n, mr * xr) - mr * Hnp(n, xr) * Jn(n, mr * xr)
	#		a_n = aNumerator / aDenominator
	#
	#		bNumerator =  mr * Jn(n, xr) * Jnp(n, mr * xr) - Jnp(n, xr) * Jn(n, mr * xr) 
	#		bDenominator = mr * Hn(n, xr) * Jnp(n, mr * xr) -  Hnp(n, xr) * Jn(n, mr * xr)
	#		b_n = bNumerator / bDenominator
	
			return _a_n(n, mr = mr, xr = xr) + _b_n(n, mr = mr, xr = xr)
	
	
		# the sum of Real{a_0 + b_0}
		#a0b0 = _nthTerm(n = 0)
		a0 = _a_n(n = 0)
		b0 = _b_n(n = 0)
	
		# now calculating sum for Real(a_n + b_n) for 1 < n < infty
		#anbn = _series_approximation(func = _nthTerm, nStart = 1, tolerance = 0.000001)
		an = _series_approximation(func = _a_n, nStart = 1, tolerance = 0.000001)
		bn = _series_approximation(func = _b_n, nStart = 1, tolerance = 0.000001)
	
		C_s = 2 * wavelength * nwLength / np.pi * (a0 + 2 * an)
		C_p = 2 *  wavelength * nwLength / np.pi * (b0 + 2 * bn)
		#C_ext = wavelength * nwLength /np.pi * (a0b0 + 2 * anbn)
	
		C_ext = 1/2 * (C_s + C_p)
		return C_ext 

	k = 2 * np.pi / wavelength
	xr = k * radius
	C_ext = extinction_coefficient(wavelength = wavelength, nwLength = nwLength, mr = mr, xr = xr)
	return np.exp(- n_s * C_ext)
	
	
# values taken from the Bellew paper 
def gen_bellew_resistances(k = 1):
	contactResistances = 		[2.5, 10, 20, 30, 40, 50, 60, 220, 280]
	contactResistanceWeights =	[6,   15,  3,  1,  2,  1,  1,   1,   1]
	if k == 1:
		return random.choices(contactResistances, weights = contactResistanceWeights, k = k)[0]
	else:
		return random.choices(contactResistances, weights = contactResistanceWeights, k = k)

def current_flow_betweenness_centrality_samples_required(
	errorTolerance = None,
	failureProbability = None,
	numberOfNodes = None
	) :

	assert type(errorTolerance) == float, 'errorTolerance must be a float between 0 and 1'
	assert 0 < errorTolerance < 1, 'errorTolerance must be a float between 0 and 1'
	assert type(failureProbability) == float, 'failureProbablity must be a float between 0 and 1'
	assert 0 < failureProbability < 1, 'failureProbability must be a float between 0 and 1'
	assert type(numberOfNodes) in [int, float], 'numberOfNodes must be a number greater than 2'
	assert numberOfNodes > 2, 'numberOfNodes must be a number greater than 2'

	eps = errorTolerance
	t = failureProbability 
	n = numberOfNodes

	return np.ceil(
		-1/2 * n**2 / (eps**2 * (n-2)**2) * np.log(t / 2)
		)

def linear_system(g,
				  current = 1
				 ):
	'''
	args
		g = input graph/circuit
		current = current injected into top electrode and withdrawn from bottom electrode
		
	returns
		adm = admittance matrix
		I = current injection vector
		nodelist = list of nodes with positions corresponding to their positions in adm and I
	'''
	sparse = scipy.sparse
	# removing non_percolating wires
	g.find_percolating_cluster(verbose = False)
	nodelist = g.percolatingCluster
	# note that subgraphing does not work on nanowiremesh objects so we have to make a placeholder networkx graph
	# and subgraph this instead
	gg = nx.Graph(g)
	sg = gg.subgraph(nodelist)
	for edge in sg.edges:
		sg.edges[edge]['admittance'] = 1/sg.edges[edge]['resistance']

	# creating admittance matrix
	adj = nx.adjacency_matrix(sg, nodelist = nodelist, weight = 'admittance')
	diag_elements = [
		sum([sg.edges[e]['admittance'] for e in sg.edges(n)]) for n in nodelist
		]
	diag = sparse.diags(diag_elements)
	adm = diag - adj

	# creating current vector
	I = np.zeros(adm.shape[0])
	top_elec_idx = nodelist.index(g.topElectrode)
	bot_elec_idx = nodelist.index(g.bottomElectrode)
	I[top_elec_idx] = current
	I[bot_elec_idx] = -current
	
	# return A and b from Ax=b where A=admittance matrix and b = current vector
	return adm, I, nodelist