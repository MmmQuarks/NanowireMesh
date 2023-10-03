import NanowireMesh as nwm
import pdb
import pandas as pd
import networkx as nx
import random
import gnuplotlib as gp
import numpy as np
from scipy import stats
from copy import deepcopy
import csv
import itertools

resultsFolder = '/Users/adamtrebach/Documents/Research/TPV/Nanowires/Results/bc_paper_comparing_rand_to_rand_11/'

makeEvolutionData = True
makeGenData = True
simOptions = dict(removeEntireWires = True,
			addInternalResistance = True,
			cleanNetwork = True,
			saveLastDataPoint = True)

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

for boolOptions in itertools.product([True, False], repeat = len(simOptions)):
	# setting all of our switches
	for n, key in enumerate(simOptions.keys()):
		simOptions[key] = boolOptions[n]

	# skipping any permutations where we are either removing wire segments
	# or are not adding internal resistors
	if not simOptions['removeEntireWires'] or not simOptions['addInternalResistance'] or not simOptions['saveLastDataPoint']:
		continue
	startPercMultiple = 2
	params = {'width' : 100, 
			'height' : 100,
			'nwLength' : 10,
			'rcMean' : 10,
			'rcSD' : 0,
			'buffer' : 1,
			'percMultiple' : startPercMultiple,
			'addInternalResistance' : simOptions['addInternalResistance']}
	
	
	#making data containers
	if makeEvolutionData:
		evolveData = pd.DataFrame(columns = ['percMultiple', 'resistance'])
		
		
		for i in range(10):
			# making network for evolution
			g = nwm.NanowireMesh(**params)
			h = deepcopy(g)
			
			counter = 0
			isPercolating = nx.has_path(g, g.topElectrode, g.bottomElectrode)
			removedNodes = []
			while isPercolating:
				nonElectrodeNodes = [node for node in g if node not in [g.topElectrode, g.bottomElectrode]]
				# shuffling the order of the nodes
				random.shuffle(nonElectrodeNodes)
	
				# getting nodes from same segment
				if simOptions['removeEntireWires'] == True:
					nodesToRemove = get_wire_segments(g, nonElectrodeNodes[0])
				else:
					nodesToRemove = set(nonElectrodeNodes[:1])
				# if removeEntireWires is false, we just return a random segment from nonElectrodeNodes	
	
				# removing nodes from graph
				poppedNodes, poppedEdges = g.pop_nodes_with_edges(nodesToRemove)
				isPercolating = nx.has_path(g, g.topElectrode, g.bottomElectrode)
				
				# saving data periodically or if this is the last iteration
				if not isPercolating or counter % 40 == 0:
					# adding back in the last removed nodes so the circuit can be solved
					# if we have just destroyed conductivity
					if not isPercolating:
						# if we are not saving the last data point, there is no need to proceed.
						if not simOptions['saveLastDataPoint']:
							break
						g.add_nodes_from(poppedNodes)
						g.add_edges_from(poppedEdges)
					if simOptions['cleanNetwork'] == True:
						clean_network(g)
					# solving the circuit and adding the data to our list
					try:
						g.solve_circuit_using_xyce()
					except csv.Error as e:
						print(e)
						print(simOptions)
						pdb.set_trace()
					evolveData = evolveData.append(dict(percMultiple = g.get_perc_multiple(),
								resistance = g.sheetResistance),
								ignore_index = True)
	
					# breaking the while loop if the sample is not percolating
					if not isPercolating:
						break
			
				counter +=1
		
		# saving data from evolutions to CSV
		evolveData.to_csv(path_or_buf = resultsFolder + str(simOptions) + 'evolveData.csv')
		
		
	else:
		# if not making data, read the existing data
		evolveData = pd.read_csv(resultsFolder + str(simOptions) + 'evolveData.csv')
	
	if makeGenData:
		# making generated data
		genData = pd.DataFrame(columns = ['percMultiple', 'resistance'])
		percMultipleList = list(np.linspace(1.0, startPercMultiple, 100))
		percMultipleList = percMultipleList + 2 * [pm for pm in percMultipleList if pm < 1.4]
		for percMultiple in percMultipleList:
			for repetitions in range(2):
				params['percMultiple'] = percMultiple
				g = nwm.NanowireMesh(**params)
				
				# cleaning the network if requested
				if simOptions['cleanNetwork'] == True:
					clean_network(g)
	
				g.solve_circuit_using_xyce()
				genData = genData.append(dict(percMultiple = g.get_perc_multiple(),
							resistance = g.sheetResistance),
							ignore_index = True)
		
		genData.to_csv(path_or_buf = resultsFolder + str(simOptions) + 'genData.csv')
	else:
		#if not making data then read the existing folder
		genData = pd.read_csv(resultsFolder + str(simOptions) + 'genData.csv')
	
		
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
	evolveCurveOptions = {'with' : 'yerrorbars linestyle 7 lw 3 lc \"black\"',
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
	plotOptions = {'xrange' : str(startPercMultiple) + ':1',
			'yrange' : '0:400',
			'title' : str(simOptions),
			'hardcopy' : resultsFolder + str(simOptions) + 'comparison_plot.png'}
	
	gp.plot(genCurve, genCurveMean, evolveCurve,  **plotOptions)
