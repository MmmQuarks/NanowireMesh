import networkx as nx
import numpy as np
import pdb
import pandas as pd
import NanowireMesh as nwm
from copy import deepcopy
import subprocess
import random
from scipy import interpolate, optimize, stats

# parameters -----------

shouldGenerateNetworks = True
shouldRedrawNetworks = True
shouldMakePlot = True
normalizeResistances = True

resultsFolder = '/Users/adamtrebach/Documents/Research/TPV/Nanowires/Results/bc_paper_node_removal_grid_04/'
removalMethodsList = ['rand', 'lowBC', 'highBC', 'lowBCRecalculated', 'highBCRecalculated']
removalMethodsList = set(removalMethodsList) - {'lowBCRecalculated', 'highBCRecalculated'}
removalMethodsList = list(removalMethodsList)

pickleName = resultsFolder + 'grid_network.p'

# ------------

def sort_nodes_by_attribute(G,
				attribute,
				includeElectrodes = False):
	# this function returns the non-electrode nodes from the graph G sorted in decreasing order of some attribute
	attrDict = nx.get_node_attributes(G, attribute)
	if attrDict == {}:
		print('Attribute \"', attribute, '\" not found in any nodes of graph. Aborting sort.')
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
		sortedInds = np.argsort(vals)
		keys = keys[sortedInds]
		keys = list(keys)
		sortedNodes = keys

	#sortedNodes is the nodes sorted in descending order of betweenness centrality 
	return sortedNodes

# calculated the current percolation multiple of the network using the total lengths of wires (excluding electrodes) remaining on the sample
def get_perc_multiple(g):
	# NL = total length of wires on sample excluding electrodes
	# A = sample area
	# L = length of single wire
	nodesExcludingElectrodes = set(g.nodes) - {g.topElectrode, g.bottomElectrode}
	nodeLengths = [g.nodes[node]['length'] for node in nodesExcludingElectrodes]
	NL = sum(nodeLengths)
	A = g.width * g.height
	L = g.nwLength
	return NL * L / (5.63726 * A)

if shouldGenerateNetworks:
	# we want a 100 x 100 network where the wires make grids 
	height = width = 10
	
	# making horizontal wires
	wires = []
	for n in range(21):
		wire = dict(x = width / 2,
				y = 0.5 * n,
				angle = 0,
				length = width *  1.05,
				temp = 298.15)
		if n == 0:
			wire['isBottomElectrode'] = True
		if n == 10:
			wire['isTopElectrode'] = True
		wires.append(wire)
	
	# making vertical wires
	for n in range(21):
		wire = dict(x = 0.5 * n,
				y = height / 2,
				angle = np.pi/2,
				length = height * 1.05,
				temp = 298.15)
		if n == 0:
			wire['x'] = wire['x'] 
		wires.append(wire)
	
	
	wireDict = {n : wire for n, wire in enumerate(wires)}
	g = nwm.NanowireMesh(width = width,
				height = height,
				nwLength = width,
				rcMean = 10,
				rcSD = 0,
				wireDict = wireDict,
				buffer = 1,
				addInternalResistance = True)

	# calculate average nw length so the get_perc_multiple function will work correctly
	lengths = nx.get_node_attributes(g, 'length')
	meanLength = np.average([val for key, val in lengths.items() if key not in [g.topElectrode, g.bottomElectrode]])
	g.nwLength = meanLength
	
	g.solve_circuit_using_xyce()
	
	
	c = nx.betweenness_centrality(g, weight = 'resistance')
	nx.set_node_attributes(g, c, 'betweenness_centrality')
	
	for removalMethod in removalMethodsList:
		onLastNode = False
		nodesToRemovePerStep =  1
		h = deepcopy(g)
	
		if removalMethod == 'rand':
			pass
		else:
			centralitySortedNodes = sort_nodes_by_attribute(h, attribute = 'betweenness_centrality', includeElectrodes = False)
		
		# saving the sheet resistance data
		percMultipleList = [get_perc_multiple(h)]
		resistance = [h.sheetResistance]
		
		# saving the initial configuration
		newPickleEnding = removalMethod + '_iter' + '0'.zfill(5) + '.p'
		thisPickleName = pickleName.replace('.p', newPickleEnding)
		h.to_pickle(outFileName = thisPickleName)
		h.to_img(showElectrodes = True, showJunctions = True,
				outFile = thisPickleName.replace('.p','plot'),
				openImage = False,
				title = 'Removal Method: ' + removalMethod + ' Iteration 0')
	
		onLastNode = False
		n = 0
		while not onLastNode:
			n += 1
			print('method:', removalMethod, '\niteration:', n)
	
			if removalMethod == 'rand':
				nodesWithoutElectrodes = set(h.nodes) - {h.topElectrode, h.bottomElectrode}
				nodesToRemove = random.sample(nodesWithoutElectrodes, nodesToRemovePerStep)
			else:
				# recalculate if called for
				if 'Recalculated' in removalMethod:
					c = nx.betweenness_centrality(h, weight = 'resistance')
					nx.set_node_attributes(h, c, 'betweenness_centrality')
				# tried to DRY here but may be confusing
				centralitySortedNodes = sort_nodes_by_attribute(h,
										attribute = 'betweenness_centrality',
										includeElectrodes = False)
				if 'low' in removalMethod:
					nodesToRemove = centralitySortedNodes[0:nodesToRemovePerStep]
				elif 'high' in removalMethod:
					nodesToRemove = centralitySortedNodes[-nodesToRemovePerStep:]
					# we reverse this so the nodes are in order from highest centrality to lowest
					nodesToRemove.reverse()
	
	
			for node in nodesToRemove:
				removedNodeData, removedEdgeData = h.pop_nodes_with_edges(node)
				h.isPercolating = nx.has_path(h, h.topElectrode, h.bottomElectrode)
				# if the network is percolating, keep removing nodes
				# if the network is no longer percolating, set h = networkBeforeNodeRemoval and solve this last circuit
				# and save network and stop removing nodes
				if h.isPercolating:
					pass
				else:
					# add back in removed nodes and edges
					h.add_edges_from(removedEdgeData)
					h.add_nodes_from(removedNodeData)
					onLastNode = True
					h.lastNode = node
					# stop removing nodes
					break
			try:
				# remove dangling ends
				print('Removing dangling ends')
				danglingEnds = {node for node in h if h.degree(node) in [0,1]}
				danglingEnds -= {h.topElectrode, h.bottomElectrode}
				while danglingEnds:
					h.remove_nodes_from(danglingEnds)
					danglingEnds = {node for node in h if h.degree(node) in [0,1]}
					danglingEnds -= {h.topElectrode, h.bottomElectrode}
				
				print('Removing non percolating nodes')
				h.find_percolating_cluster()
				nonPercolating = [node for node in h if node not in h.percolatingCluster]
				h.remove_nodes_from(nonPercolating)
			                                                                                                                                          
				# solve circuit and update
				h.solve_circuit_using_xyce()
				# save the sheet resistance data
				percMultipleList.append(get_perc_multiple(h))
				resistance.append(h.sheetResistance)
			                                                                                                                                          
				# save the pickle of this network every 30 samples or if we're on the last node
				iterNumber = str(n)
				newPickleEnding = removalMethod + '_iter' + iterNumber.zfill(5) + '.p'
				thisPickleName = pickleName.replace('.p', newPickleEnding)
				h.to_pickle(outFileName = thisPickleName)
				h.to_img(showElectrodes = True, showJunctions = True,
					outFile = thisPickleName.replace('.p','plot'),
					openImage = False,
					title = 'Removal Method: ' + removalMethod + ' Iteration ' + iterNumber,
					lineWidth = 14,
					pointSize = 7)
	
	
	
				if onLastNode:
					# if the network is about to stop percolating, we are done. We want to record the last snapshot here and move on.
					print('Circuit is one node away from failure. Finished with simulations for method', removalMethod)
					# because we have set onLastNode = True earlier in the code, this is the last thing 
					# that will be executed in the while loop
			except subprocess.CalledProcessError as err:
				print('SHIT GOT FUCKED')
				print(err.output)
			
		# saving the data as csv
		data = [(percMultipleList[i] , resistance[i]) for i in range(len(percMultipleList))]
		data = np.array(data)
		fileEnding = removalMethod + '_data'
		fileName = pickleName.replace('.p', fileEnding)
		# so fileName will be resultsFolder/grid_networkremovalMethod_data.csv
		np.savetxt(fileName, data, fmt = '%f', delimiter = ' ')

		if normalizeResistances:
			data[:,1] = data[:,1] / data[0,1]
			np.savetxt(fileName + '_normalized', data, fmt = '%f', delimiter = ' ')

	


if shouldMakePlot:
	plotScriptFileName = resultsFolder + 'plot_script.txt'

	colorDict = {'rand' : 'black',
			'highBC' : 'red',
			'lowBC' : 'blue',
			'lowBCRecalculated': 'dark-green',
			'highBCRecalculated' : 'orange'}
	
	plotScript = ['set term png size 8000, 6000 linewidth 15',
			'set output \"' + plotScriptFileName.replace('plot_script.txt', 'evolution_plot.png\"'),
			'set title \"Network Evolution On Square Grid Under Node Removal \" font \",200\" offset 35,15',
			'set tmargin 50',
			'set lmargin 60',
			'set rmargin 230',
			'set bmargin 40',
			'set xlabel \"Network Density as Multiple of Percolation Threshold\" font \",170\" offset 0,-30',
			'set ylabel \"Normalized Series Resistance of Network\" font \",170\" offset -27,0',
			'set xrange [0.4:0]',
			'set yrange [0:20]',
			'set xtics font \",140\" offset 0,-12',
			'set ytics font \",140\" offset -2,0',
			'set key outside right top font \",170\"']
	
	#plotScript += ['set key off']
	# calculating the densities at 50% to failure
	halfwayData = [['removalMethod', 'halfwayRelativeDensityDecrease', 'halfwayRelativeResistanceIncrease']]

	#making container for the summary data
	summaryDf = pd.DataFrame(columns = ['removalMethod', 'RDoublePercMult', 'EndPercMult'])

	subplotCommandsList = []
	for n, removalMethod in enumerate(removalMethodsList):
		fileEnding = removalMethod + '_data'
		fileName = pickleName.replace('.p', fileEnding) + ('_normalized' if normalizeResistances else '')
		thisHalfwayData = np.loadtxt(fileName)
		halfwayDensity = np.average( thisHalfwayData[[0,-1],0])
		initialDensity = thisHalfwayData[0,0]
		# we use np.flip because np.interp requires x data to be in increasing order
		halfwayResistance = np.interp(halfwayDensity, np.flip(thisHalfwayData[:,0]), np.flip(thisHalfwayData[:,1]))
		halfwayData.append([removalMethod, halfwayDensity/initialDensity, halfwayResistance])
	
		command = ('plot \"' if n == 0 else '\"') + fileName + '\" using 1:2 with lines title \"' + removalMethod + '\" lc \"' + colorDict[removalMethod] + '\" lw 8'
		subplotCommandsList.append(command)
	
		# adding the markers for the final density at the end of the plot
		df = pd.read_csv(fileName, sep = ' ', names = ['percMultiple', 'resistance'])
		minPercMultiple = df['percMultiple'].min()	
		maxResistance = df['resistance'].max()
		minMarker = np.array([(minPercMultiple, 10)])
		minMarker = np.array([(minPercMultiple, 0), (minPercMultiple, maxResistance)])
		minMarkerFileName = fileName + 'min_percMultiple'
		np.savetxt(minMarkerFileName, minMarker, fmt = '%f', delimiter = ' ')
		command = '\"' + minMarkerFileName + '\" using 1:2 with lines lw 5 lc \"' + colorDict[removalMethod] + '\" notitle'
		subplotCommandsList.append(command)

		# calculating the density at which resistance doubles
		percToResFunc = interpolate.interp1d(df['percMultiple'], df['resistance'])
		funcToOptimize = lambda x : percToResFunc(x) - 2 * df['resistance'].min()
		rootResults = optimize.root_scalar(funcToOptimize,
							bracket = (df['percMultiple'].min(), df['percMultiple'].max()))
		if rootResults.converged:
			summaryDf = summaryDf.append(
							dict(RDoublePercMult = rootResults.root, 
								EndPercMult = df['percMultiple'].min(),
								removalMethod = removalMethod),
							ignore_index = True)

	summaryDf.to_csv(resultsFolder + 'summary.csv', index = False)

	subplotCommand = ', '.join(subplotCommandsList)
	plotScript.append(subplotCommand)
	plotScriptString = '\n'.join(plotScript)
	plotScriptFile = open(plotScriptFileName , 'w')
	plotScriptFile.write(plotScriptString)
	plotScriptFile.close()
	
	try:
		subprocess.check_call(['gnuplot', plotScriptFileName])
	except subprocess.CalledProcessError as err:
		print('Error in Gnuplot Script')
		print(err.output)
	
	for l in halfwayData:
		print(l)
