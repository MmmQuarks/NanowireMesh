import NanowireMesh as nwm
import pdb
import pandas as pd
import networkx as nx
import random
import gnuplotlib as gp
import numpy as np
from scipy import stats
from copy import deepcopy
import time
import csv
import itertools
import ParallelCentrality as PC
import sys

resultsFolder = '/Users/adamtrebach/Documents/Research/TPV/Nanowires/Results/bc_paper_node_removal_data_v3_005/'
xycePath = 'xyce'
makeEvolutionData = False #toggle to evolve the networks
makeGenData = False # toggle to make "ground truth" resistance data


simOptions = dict(removeEntireWires = True,
			addInternalResistance = True,
			cleanNetwork = True,
			numTrials = 10,
			percMultiples = [float(val) for val in sys.argv[1:]],
			netlistSuffix = '_'.join([val for val in sys.argv[1:]]),
			nodesToRemovePerStep = 40,
			removalMethodsToCalculate = ['rand'],
#						'lowBC',
#						'lowEC',
#						'lowPC',
#						'lowCW',
#						'lowPW'],
			removalMethodParams = {'rand' : dict(color = 'black', centralityName = 'random'),
				'lowBC' : dict(color = 'light-red', centralityName = 'betweenness_centrality'),
				'highBC' : dict(color = 'dark-red', centralityName = 'betweenness_centrality'),
				'lowEC' : dict(color = 'web-green', centralityName = 'electrode_centrality'),
				'highEC' : dict(color = 'dark-green', centralityName = 'electrode_centrality'),
				'lowPC' : dict(color = 'web-blue', centralityName = 'percolation_centrality'),
				'highPC' : dict(color = 'dark-blue', centralityName = 'percolation_centrality'),
				'lowDownPotentialPC' : dict(color = 'sienna1', centralityName = 'down_potential_percolation_centrality'),
				'highDownPotentialPC' : dict(color = 'dark-orange', centralityName = 'down_potential_percolation_centrality'),
				'lowCW' : dict(color = 'orchid', centralityName = 'current_weighted_centrality'),
				'highCW' : dict(color = 'orchid4', centralityName = 'current_weighted_centrality'),
				'lowPW' : dict(color = 'khaki1', centralityName = 'power_weighted_centrality'),
				'highPW' : dict(color = 'gold', centralityName = 'power_weighted_centrality'),
				'lowAngle' : dict(color = 'light-pink', centralityName = 'absolute_angle'),
				'highAngle' : dict(color = 'dark-pink', centralityName = 'absolute_angle')}
					)
# setting network parameters	
params = {'width' : 100, 
		'height' : 100,
		'nwLength' : 10,
		'rcMean' : 10,
		'rcSD' : 0,
		'buffer' : 1,
		'percMultiple' : None,
		'addInternalResistance' : simOptions['addInternalResistance']}

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

#
#def clean_network(g):
#	# removing isolates and dangling ends
#	nodesToRemove = {node for node in g.nodes if g.degree(node) in [0,1]}
#	nodesToRemove -= {g.topElectrode, g.bottomElectrode}
#	while nodesToRemove:
#		g.remove_nodes_from(nodesToRemove)
#		nodesToRemove = {node for node in g.nodes if g.degree(node) in [0,1]}
#		nodesToRemove -= {g.topElectrode, g.bottomElectrode}
#	# removing non percolating nodes
#	try:
#		g.find_percolating_cluster()
#	except KeyError:
#		pdb.set_trace()
#	nonPercolating = set(g.nodes) - set(g.percolatingCluster)
#	g.remove_nodes_from(nonPercolating)

class ExceptionMsg(Exception):
	def __init__(self, msg):
		self.msg = msg

	def __str__(self):
		return(repr(self.msg))


def sort_nodes_by_attribute(G,
				attribute,
				includeElectrodes = False):
	# this function returns the non-electrode nodes from the graph G sorted in increasing  order of some attribute
	attrDict = nx.get_node_attributes(G, attribute)

	if attrDict == {}:
		sortedNodes = None
	else:
		# removing (or not removing) electrodes)
		if includeElectrodes:
			pass
		elif not includeElectrodes:
			del attrDict[G.topElectrode]
			del attrDict[G.bottomElectrode]
		# make the sorted list
		keys = np.array(list(attrDict.keys()))
		vals = list(attrDict.values())
		try:
			sortedInds = np.argsort(vals)
		except TypeError:
			pdb.set_trace()
		keys = keys[sortedInds]
		keys = list(keys)
		sortedNodes = keys

	#sortedNodes is the nodes sorted in increasing order of attribute
	# we confirm this in the below assert statement
	for n in range(len(sortedNodes) - 1):
		thisNode = sortedNodes[n]
		nextNode = sortedNodes[n+1]
		msg = ' '.join(['Error: attribute \"', 
					attribute,
					'\" of node',
					str(thisNode),
					'is not greater than that of node',
					str(nextNode)])
		assert G.nodes[thisNode][attribute] <= G.nodes[nextNode][attribute], msg
	# as long as we pass the assertion test, we return sortedNodes
	return sortedNodes


for percMultiple in simOptions['percMultiples']:
	# making data container
	evolveData = pd.DataFrame(columns = ['percMultipleGroup', 'trial', 'removalMethod', 'percMultiple', 'resistance'])
	for trial in range(simOptions['numTrials']):

		if makeEvolutionData:
			params['percMultiple'] = percMultiple
			# making network for evolution
			g = nwm.NanowireMesh(**params)
			g.solve_circuit_using_xyce(xycePath = xycePath, netlistName = 'netlist' + simOptions['netlistSuffix'], voltage = 1)

			# iterate through removalMethodsToCalculate and find the centralities 
			# we need to calculate 
			centralitiesToCalculate = set()
			for removalMethod in simOptions['removalMethodsToCalculate']:
				centralityName = simOptions['removalMethodParams'][removalMethod]['centralityName']
				centralitiesToCalculate.add(centralityName)

			print('Calculating Centralities')
			for centrality in centralitiesToCalculate:
				if centrality == 'random':
					# if random, assign random floats to all nodes 
					centDict = {node : np.random.rand() for node in g}
				elif centrality == 'betweenness_centrality':
					centDict = nx.betweenness_centrality(g, weight = 'resistance')
				elif centrality == 'electrode_centrality':
					centDict = nwm.electrode_centrality(g, potential = 'voltage', weight = 'resistance')
				elif centrality == 'percolation_centrality':
					centDict = nx.percolation_centrality(g, attribute = 'voltage', weight = 'resistance')
				elif centrality == 'current_weighted_centrality':
					centDict = PC.current_weighted_centrality(g)
				elif centrality == 'power_weighted_centrality':
					centDict = PC.power_weighted_centrality(g)
				else: # if centrality doesn't match 
					try:
						raise ExceptionMsg('Centrality \"' + centrality + '\" not associated with calculation function.')
					except ExceptionMsg as err:
						print(err)
						raise Exception('Ending Run')
				# storing centralities in graph
				nx.set_node_attributes(g, centDict, centrality)

			# store graph, now with centralities, in results folder
			pickleName = resultsFolder + '_'.join([str(params['width']) + 'x' + str(params['height']),
						'percMult' + str(percMultiple),
						'nwLen' + str(params['nwLength']),
						'Rc' + str(params['rcMean']),
						'trial' + str(trial),
						'initial.p'])
			g.to_pickle(outFileName = pickleName)

			for centrality in centralitiesToCalculate:
				counter = 0 # counts how many iterations of node removals
				h = deepcopy(g)
				isPercolating = nx.has_path(h, h.topElectrode, h.bottomElectrode)
				removedNodes = []
				onLastNode = False
	
				# iteratively removing nodes until percolation is destroyed
				while isPercolating:
					# this returns nodes sorted by INCREASING value of attribute
					# removing nodes one by one and ending the sim if we have destroyed percolation
					for n in range(simOptions['nodesToRemovePerStep']):
						sortedNodes = sort_nodes_by_attribute(h, attribute = centrality, includeElectrodes = False)
						node = sortedNodes[0] # sorted in decreasing order of centrality so lowest centrality is last node
						# making sure that we remove the entire wire of which this node is a segment
						allSegmentNodes = get_wire_segments(h, node)
						removedNodeData, removedEdgeData = h.pop_nodes_with_edges(allSegmentNodes)
						isPercolating = nx.has_path(h, h.topElectrode, h.bottomElectrode)
						if not isPercolating:
							# if not percolating, add back in removed nodes and break the removal loop
							h.add_nodes_from(removedNodeData)
							h.add_edges_from(removedEdgeData)
							break
					# calculating resistance
					if simOptions['cleanNetwork'] == True:
						clean_network(h)
					h.solve_circuit_using_xyce(xycePath = xycePath, netlistName = 'netlist' + simOptions['netlistSuffix'], voltage = 1)
					evolveData = evolveData.append({'percMultipleGroup': percMultiple,
									'trial' : trial,
									'removalMethod' : centrality,
									'percMultiple' : h.get_perc_multiple(),
									'resistance' : h.sheetResistance},
									ignore_index = True)

					# every 10 iterations and also on the last iteration, we store a snapshot of the network
					if counter % 5 == 0 or not isPercolating:
						pickleName = resultsFolder + '_'.join([str(params['width']) + 'x' + str(params['height']),
							'percMult' + str(percMultiple),
							'nwLen' + str(params['nwLength']),
							'Rc' + str(params['rcMean']),
							'trial' + str(trial),
							centrality,
							'iter' + str(counter).zfill(5) + '.p'])
						h.to_pickle(outFileName = pickleName)
					counter += 1

				# saving data from evolutions to CSV
				evolveData.to_csv(path_or_buf = resultsFolder + 'percMultiple' + str(percMultiple) + '_evolutionData.csv') 
			
		else:
			# if not making data, read the existing data
			evolveData = pd.read_csv(resultsFolder + 'percMultiple' + str(percMultiple) + '_evolutionData.csv')
	
	
# block to make ground truth

if makeGenData:
	# making generated data
	genData = pd.DataFrame(columns = ['percMultiple', 'resistance'])
	percMultipleList = list(np.linspace(1.0, 2.2, 75))
#	percMultipleList = percMultipleList + 2 * [pm for pm in percMultipleList if pm < 1.4]
	for percMultiple in percMultipleList:
		maxReps = 25
		minReps = 5
		minPerc, maxPerc = min(percMultipleList), max(percMultipleList)
		numReps = maxReps + (minReps - maxReps) / (maxPerc - minPerc) * (percMultiple - minPerc)
		numReps = int(np.ceil(numReps))
		for repetitions in range(numReps):
			params['percMultiple'] = percMultiple
			g = nwm.NanowireMesh(**params)
			
			# cleaning the network if requested
			if simOptions['cleanNetwork'] == True:
				clean_network(g)

			g.solve_circuit_using_xyce()
			genData = genData.append(dict(percMultiple = g.get_perc_multiple(),
						resistance = g.sheetResistance),
						ignore_index = True)
	
	genData.to_csv(path_or_buf = resultsFolder + 'genData.csv')
else:
	#if not making data then read the existing folder
	genData = pd.read_csv(resultsFolder + 'genData.csv')

# filtering genData to not include ugly data below percMultiple = 1
genData = genData[genData['percMultiple'] > 0.75]
#genData = genData[genData['percMultiple'] < evolveData['percMultiple'].max()]

for percMultiple in simOptions['percMultiples']:
	for removalMethod in simOptions['removalMethodsToCalculate']:
#		pdb.set_trace()
		evolveData = pd.read_csv(resultsFolder + 'percMultiple' + str(percMultiple) + '_evolutionData.csv')
		centralityName = simOptions['removalMethodParams'][removalMethod]['centralityName']
		evolveData = evolveData[evolveData['removalMethod'] == centralityName]
		evolveData = evolveData[evolveData['percMultiple'] > 0.75]
		# making binned statistics
		evolveResMeans, evolveResMeanBinEdges, evolveResMeanBinNumbers = stats.binned_statistic(evolveData['percMultiple'],
														evolveData['resistance'],
														statistic = 'mean',
														bins = 15)
		evolveResSems, evolveResSemBinEdges, evolveResSemBinNumbers = stats.binned_statistic(evolveData['percMultiple'],
														evolveData['resistance'],
														statistic = stats.sem,
														bins = 15)
		genResMeans, genResMeanBinEdges, genResMeanBinNumbers = stats.binned_statistic(genData['percMultiple'],
														genData['resistance'],
														statistic = 'mean',
														bins = 10)
		genResSems, genResSemBinEdges, genResSemBinNumbers = stats.binned_statistic(genData['percMultiple'],
														genData['resistance'],
														statistic = stats.sem,
														bins = 10)
		
		
		# formatting the evolution data
		evolveCurveOptions = {'with' : 'yerrorbars linestyle 7 lw 5 lc \"black\"',
					'tuplesize' : 3}
		evolveCurvePercMultiples = [np.average([evolveResMeanBinEdges[n], evolveResMeanBinEdges[n+1]]) for n in range(len(evolveResMeanBinEdges) - 1)]
		evolveCurveDf = pd.DataFrame(dict(percMultiple = evolveCurvePercMultiples,
						resistanceMean = evolveResMeans,
						resistanceSem = evolveResSems))
		evolveCurve = (evolveCurveDf['percMultiple'], 
					evolveCurveDf['resistanceMean'],
					evolveCurveDf['resistanceSem'],
					evolveCurveOptions)
		
		# formatting the gen data
		genCurveOptions = {'with' : 'filledcurves lc \"orange\"',
					'tuplesize' : 3}
		
		genCurvePercMultiples = [np.average([genResMeanBinEdges[n], genResMeanBinEdges[n+1]]) for n in range(len(genResMeanBinEdges) - 1)]
		genCurveDf = pd.DataFrame(dict(percMultiple = genCurvePercMultiples,
					resistanceMean = genResMeans,
					resistanceSem = genResSems))
		genCurve = (genCurveDf['percMultiple'], 
					genCurveDf['resistanceMean'] - genCurveDf['resistanceSem'],
					genCurveDf['resistanceMean'] + genCurveDf['resistanceSem'],
					genCurveOptions)
		
		genCurveMeanOptions = {'with' : 'lines lw 3 lc \"blue\"'}
		genCurveMean = (genCurveDf['percMultiple'],
					genCurveDf['resistanceMean'],
					genCurveMeanOptions)
		
		# formatting plot options
		plotOptions = {'xrange' : str(2.5) + ':0.4',
				'yrange' : '0:400',
#				'title' : 'Removing lowest ' + centralityName.replace('_',' ') + ' starting at pm ' + str(percMultiple),
				'hardcopy' : resultsFolder + 'percMultiple' + str(percMultiple) + removalMethod + 'comparison_plot.png',
				'unset' : ['grid', 'title'],
				'cmds' : ['set xtics scale 3 font \",40\" offset 0,-2',
						'set ytics scale 3 font \",40\" offset -2,0',
						'set tmargin 4',
						'set lmargin 15',
						'set bmargin 5']}
		
		gp.plot(genCurve, genCurveMean, evolveCurve,  **plotOptions)
