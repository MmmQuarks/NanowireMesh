import NanowireMesh as nwm
import pdb
import pandas as pd
import networkx as nx
import os
import random
import numpy as np
import ParallelCentrality as PC
import sys

# allowing for customization without changing this script or the submit script 
resultsFolder = sys.argv[1]
if resultsFolder[-1] != '/':
	resultsFolder = resultsFolder + '/'

taskID = sys.argv[2]

# this dict of a single dict is a little more complex than needed here, but it is just preserving existing structures that 
# work well (i.e. the version of this code that generates networks and calculates lots of centralities).
centralities = dict(lowBCUW = dict(name = 'unweighted_betweenness_centrality',
				color = 'dark-violet',
				func = lambda g : nx.betweenness_centrality(g, weight = None),
				clean = True,
				pointtype = 8
				)
		)

for networkFilePath in sorted(os.listdir(resultsFolder)):
	isPickle = networkFilePath[-2:] == '.p'
	isRightTask = '_'.join(['task' ,taskID.zfill(2)]) in networkFilePath # making sure that we are dividing calculations across nodes correctly
	if isPickle and isRightTask:
		# opening network to calculate unweighted BC on 
		h = nwm.NanowireMesh(inPickle = resultsFolder + networkFilePath)
		
		# calculating all centralities that will be needed
		for key in centralities.keys():
			print('Task', taskID, 'Calculating', centralities[key]['name'])
			try:
				c = centralities[key]['func'](h)
			except TypeError:
				pdb.set_trace()
			nx.set_node_attributes(h, c, centralities[key]['name'])
		
		# setting each node in each wire to have the same centrality as the max centrality node in that wire
		# finding all sets of wires
		contactEdges = [edge for edge in h.edges if h.edges[edge]['resistanceType'] == 'cont']
		poppedEdges = h.pop_edges(contactEdges)
		for key in centralities.keys():
			centralityName = centralities[key]['name']
			for wire in nx.connected_components(h):
				nodeCents = [h.nodes[node][centralityName] for node in wire]
				for node in wire:
					h.nodes[node][centralityName] = max(nodeCents)
		h.add_edges_from(poppedEdges)				
			
		# store graph, now with unweighted betweenness centralities calculated
		h.to_pickle(outFileName = resultsFolder + networkFilePath)
