import NanowireMesh as nwm
import re
import pdb
import pandas as pd
import networkx as nx
import random
import gnuplotlib as gp
import numpy as np
from scipy import stats, interpolate, optimize
from copy import deepcopy
import ParallelCentrality as PC
import csv
import itertools
from sortedcontainers import SortedDict
import os
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
srcNetworkFolder = runtimeOptions['srcNetworkFolder']

# making sure we have everything we need in runtime options



xycePath = '/home/amwt/TPV/2DNanowires/installations/build_Xyce-6.11_parallel_openMPI/bin/Xyce'
makeEvolutionData = True
makeGenData = True

simOptions = dict(removeEntireWires = True,
			addInternalResistance = True,
			saveLastDataPoint = True,
			trials = 10)

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
				),
		lowBCUW = dict(name = 'unweighted_betweenness_centrality',
				color = 'dark-violet',
				func = lambda g : nx.betweenness_centrality(g, weight = None),
				clean = True,
				pointtype = 8
				)
			)

if not runtimeOptions['useEigenvectorCentrality']:
	del centralities['eig']

if runtimeOptions['onlyUseUnweightedBC']:
	for key in list(centralities.keys()):
		if key != 'lowBCUW':
			del centralities[key]
	
# values taken from the Bellew paper 
def gen_bellew_resistances(k = 1, exclude_outliers = False):
	contactResistances = 		[2.5, 10, 20, 30, 40, 50, 60, 220, 280]
	contactResistanceWeights =	[6,   15,  3,  1,  2,  1,  1,   1,   1]
	if exclude_outliers:
		del contactResistances[-2:]
		del contactResistanceWeights[-2:]
	return random.choices(contactResistances, weights = contactResistanceWeights, k = k)

#keys = list(centralities.keys())
#for key in keys:
#	if 'rand' not in key:
#		del centralities[key]




def effective_perc_multiple_to_density(percMultiple, nwLength, height, width):
	assert height == width, "Height and width must be equal"
	nc_eff = 5.63726 / nwLength**2 + 1 / (nwLength * height) + 5.5 / (height**2)
	n_s = percMultiple * nc_eff
	return n_s

def fom(resistance, transparency):
	return 188.5 / (resistance * (transparency**(-1/2) - 1))


def get_wire_segments(g, startingNode):
	# private class to hold the nodes of a wire
	class _Segment:
		# Constructor
		def __init__(self):
			# set of nodes visited
			self.visited = set()

		# implementation of depth first search
		def _dfs_recursive(self, g, startingNode):
			# mark the current node as visited
			self.visited.add(startingNode)

			# recur for all nodes adjacent to this node
			# but only if they are connected by internal resistors
			for neighbor in g.neighbors(startingNode):
				if neighbor not in self.visited:
					if g.edges[startingNode, neighbor]['resistanceType'] == 'int':
						self._dfs_recursive(g, neighbor)
	segment = _Segment()
	segment._dfs_recursive(g, startingNode)
	return segment.visited

def clean_network(g):
	# removing non percolating nodes (also removes isolates)
	try:
		g.find_percolating_cluster()
	except KeyError:
		pdb.set_trace()
	nonPercolating = set(g.nodes) - set(g.percolatingCluster)
	g.remove_nodes_from(nonPercolating)


	# removing dangling ends - but must be sure to not remove
	# segments of participating wires because this introduces asymmetry between
	# random removal and random generation

	# getting all possibly dangling ends
	danglingEnds = {node for node in g.nodes if g.degree(node) == 1} 

	# filtering to remove entire dangling wires only if all segments of 
	# wire have degree less than three
	danglingEnds = {node for node in danglingEnds if all([g.degree(neighb) <= 2 for neighb in get_wire_segments(g, node)])}

	# making sure we don't remove the electrodes
	danglingEnds -= {g.topElectrode, g.bottomElectrode}

	# getting a set of all nodes to remove including the other segments connected to dangling ends
	nodesToRemove = {seg for node in danglingEnds for seg in get_wire_segments(g, node)}
	while nodesToRemove:
		g.remove_nodes_from(nodesToRemove)
		danglingEnds = {node for node in g.nodes if g.degree(node) == 1} 
		danglingEnds = {node for node in danglingEnds if all([g.degree(neighb) <= 2 for neighb in get_wire_segments(g, node)])}
		danglingEnds -= {g.topElectrode, g.bottomElectrode}
		nodesToRemove = {seg for node in danglingEnds for seg in get_wire_segments(g, node)}
#
#		nodesToRemove = {node for node in g.nodes if g.degree(node) in [0,1]}
#		nodesToRemove -= {g.topElectrode, g.bottomElectrode}

def sort_nodes_by_attribute(G, attribute, includeElectrodes = False):
	attrDict = nx.get_node_attributes(G, attribute)
	if not includeElectrodes:
		del attrDict[G.topElectrode]
		del attrDict[G.bottomElectrode]
	# returns attributes sorted by increasing key
	return SortedDict({val : key for key, val in attrDict.items()})
	
def get_perc_multiple(g):
	assert g.width == g.height, "Width and height are not the same"
	lTot = sum([g.nodes[node]['length'] for node in g if node not in [g.topElectrode, g.bottomElectrode]])
	l = g.nwLength
	Ls = g.width
	return lTot * l / (5.63726 * Ls**2 + l * Ls + 5.5 * l**2)

evolveData = pd.DataFrame(columns = ['task','trial', 'centrality', 'percMultiple', 'resistance'])

for origNetwork in os.listdir(srcNetworkFolder):
	correctTask = 'task_' + taskID.zfill(2) in origNetwork
	isPickle = origNetwork[-2:] == '.p'
	if correctTask and isPickle:
		# making network for evolution
		h = nwm.NanowireMesh(inPickle = srcNetworkFolder + origNetwork)

		# getting the trial number
		trial = re.search("trial_(\d{3})", origNetwork).group()
		trial = trial[-3:]
		# for each trial, iterating through all possible removal methods
		for key in centralities.keys():
			g = deepcopy(h)
			counter = 0
			isPercolating = nx.has_path(g, g.topElectrode, g.bottomElectrode)
			evolveData = evolveData.append(dict(task = taskID.zfill(2),
								trial = trial,
								centrality = centralities[key]['name'],
								percMultiple = get_perc_multiple(g),
								resistance = g.sheetResistance),
								ignore_index = True)

			while isPercolating:
#j				nonElectrodeNodes = [node for node in g if node not in [g.topElectrode, g.bottomElectrode]]
				# shuffling the order of the nodes
#				random.shuffle(nonElectrodeNodes)
				
				# return sorted dict of centrality : node pairs
				# sorted in increasing order of centrality
				sd = sort_nodes_by_attribute(g,
								attribute = centralities[key]['name'],
								includeElectrodes = False)
				c, node = sd.peekitem(index = 0)
				nodesToRemove = get_wire_segments(g, node) if simOptions['removeEntireWires'] else [node]
#				if simOptions['removeEntireWires'] == True:
#					nodesToRemove = get_wire_segments(g, nonElectrodeNodes[0])
#				else:
#					nodesToRemove = set(nonElectrodeNodes[:1])
				# if removeEntireWires is false, we just return a random segment from nonElectrodeNodes	
	
				# removing nodes from graph
				poppedNodes, poppedEdges = g.pop_nodes_with_edges(nodesToRemove)
				isPercolating = nx.has_path(g, g.topElectrode, g.bottomElectrode)
				
				# saving data periodically or if this is the last iteration
				if not isPercolating or (counter > 0 and counter % 40*5*5 == 0):
					# adding back in the last removed nodes so the circuit can be solved
					# if we have just destroyed conductivity
					if not isPercolating:
						# if we are not saving the last data point, there is no need to proceed.
						if not simOptions['saveLastDataPoint']:
							break
						g.add_nodes_from(poppedNodes)
						g.add_edges_from(poppedEdges)
					if centralities[key]['clean']:
						clean_network(g)
					# solving the circuit and adding the data to our list
					try:
						g.solve_circuit_using_xyce(xycePath = xycePath,
							netlistName = resultsFolder + '_'.join(['netlist',
												'task',
												taskID.zfill(2)]))
					except csv.Error as e:
						print(e)
						print(simOptions)
						pdb.set_trace()
					evolveData = evolveData.append(dict(trial = trial,
										centrality = centralities[key]['name'],
										percMultiple = get_perc_multiple(g),
										resistance = g.sheetResistance),
									ignore_index = True)

					# save to pickle
					pickleName = resultsFolder 
					pickleName += '_'.join([origNetwork.replace('_initial.p', ''), #remaining string looks like task_01_trial_0122
									centralities[key]['name'],
									'iter',
									str(counter).zfill(5) + '.p'])

					g.to_pickle(outFileName = pickleName)

					# saving data from evolutions to CSV
					# this should be done frequently so if it crashes we don't lose everything
					evolveData.to_csv(path_or_buf = resultsFolder + '_'.join(['evolveData',
													'task',
													taskID.zfill(2) + '.csv']))
	
					# breaking the while loop if the sample is not percolating
					if not isPercolating:
						break
			
				counter +=1
		

