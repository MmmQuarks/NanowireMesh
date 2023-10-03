import NanowireMesh as nwm
import pdb
import pandas as pd
import networkx as nx
import random
import numpy as np
import ParallelCentrality as PC
import sys

# allowing for customization without changing this script or the submit script 
resultsFolder = sys.argv[1]
if resultsFolder[-1] != '/':
	resultsFolder = resultsFolder + '/'
runtimeOptionsDf = pd.read_csv(resultsFolder + 'runtimeOptions.csv')
# columns 
runtimeOptions = dict(zip(runtimeOptionsDf.iloc[:,0], runtimeOptionsDf.iloc[:,1]))
for key,val in runtimeOptions.items(): # converting strings into other data types if possible
	if val in ['True', 'False']:
		runtimeOptions[key] = val == 'True'
	else:
		try:
			runtimeOptions[key] = float(val)
		except ValueError:
			pass # this is not a number

taskID = sys.argv[2]

#xycePath = 'xyce'
xycePath = '/home/amwt/TPV/2DNanowires/installations/build_Xyce-6.11_parallel_openMPI/bin/Xyce'


imgFormat = 'png'


centralities = dict(rand = dict(name = 'random',
			color = 'black',
			func = lambda g : {node : np.random.rand() for node in g},
			clean = False,
			pointtype = 4,
			),
		randCleaned = dict(name = 'random_cleaned',
			color = 'black',
			func = lambda g : {node : np.random.rand() for node in g},
			clean = True,
			pointtype = 6
			),
		lowBC = dict(name = 'betweenness_centrality',
				color = 'dark-violet',
				func = lambda g : nx.betweenness_centrality(g, weight = 'resistance'),
				clean = True,
				pointtype = 8
				),
		lowEC = dict(name = 'electrode_centrality', 
				color = 'dark-green', 
				func = lambda g : nwm.electrode_centrality(g, potential = 'voltage', weight = 'resistance'),
				clean = True,
				pointtype = 12
				),
		lowPC = dict(name = 'percolation_centrality',
				color = 'red', 
				func = lambda g : nx.percolation_centrality(g, attribute = 'voltage', weight = 'resistance'),
				clean = True,
				pointtype = 17
				),
		lowCW = dict(name = 'current_weighted_centrality', 
				color = 'dark-orange', 
				func = lambda g : PC.current_weighted_centrality(g),
				clean = True,
				pointtype = 32
				),
		lowPW = dict(name = 'power_weighted_centrality',
				color = 'dark-blue', 
				func = lambda g : PC.power_weighted_centrality(g),
				clean = True,
				pointtype = 10
				),
		eig = dict(name = 'eigenvector_centrality',
				color = 'gray50',
				func = lambda g : nx.eigenvector_centrality(g, weight = 'admittance'),
				clean = True,
				pointtype = 90
				)
			)

if not runtimeOptions['useEigenvectorCentrality']:
	del centralities['eig']

# values taken from the Bellew paper 
def gen_bellew_resistances(k = 1, exclude_outliers = False):
	contactResistances = 		[2.5, 10, 20, 30, 40, 50, 60, 220, 280]
	contactResistanceWeights =	[6,   15,  3,  1,  2,  1,  1,   1,   1]
	if exclude_outliers:
		del contactResistances[-2:]
		del contactResistanceWeights[-2:]
	return random.choices(contactResistances, weights = contactResistanceWeights, k = k)


params = {'width' : None, 
		'height' : None,
		'nwLength' : None,
		'rcMean' : None,
		'rcSD' : None,
		'buffer' : 1,
		'percMultiple' : None,
		'nwDiam' : None,
		'buffer' : 2,
		'initialTemp' : 298.15,
		'addInternalResistance' : None}

# setting the params that are specific in runtimeOptions.csv
for key in params.keys():
	if params[key] == None:
		params[key] = runtimeOptions[key]


counter = 0
while True:
	# making network for evolution
	h = nwm.NanowireMesh(**params)
	
	# modifying the contact resistances to have the same distribution as the Bellew paper
	# https://pubs.acs.org/doi/pdf/10.1021/acsnano.5b05469 
	contactEdges = [edge for edge in h.edges if h.edges[edge]['resistanceType'] == 'cont']
	contactResistances = gen_bellew_resistances(k = len(contactEdges),
							exclude_outliers = runtimeOptions['excludeOutliers'])
	nx.set_edge_attributes(h, 
				dict(zip(contactEdges, contactResistances)),
				'resistance')
	# adding admittances
	admittances = {edge: 1 / h.edges[edge]['resistance'] for edge in h.edges}
	nx.set_edge_attributes(h, admittances, 'admittance')
	print('about to solve')	
	h.solve_circuit_using_xyce(xycePath = xycePath,
					netlistName = resultsFolder + 'netlist_task_' + taskID.zfill(2),
					voltage = 1)
	# calculating all centralities that will be needed
	for key in centralities.keys():
		print('Task', taskID, 'iter', counter,  'Calculating', centralities[key]['name'])
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
		
	# store initial graph, now with centralities, in results folder
	pickleName = resultsFolder + '_'.join(['task', taskID.zfill(2), 'iter', str(counter).zfill(3), 'initial.p'])
	h.to_pickle(outFileName = pickleName)
	
	counter += 1
	
	
