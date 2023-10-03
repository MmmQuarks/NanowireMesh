import random
from copy import deepcopy
import subprocess
import NanowireMesh as nwm
import os
import networkx as nx
import numpy as np
import itertools
from scipy.special import factorial as factorial
from scipy import stats, optimize, interpolate
import pandas as pd
import time

#setting network properties
networkSize = 100
width = height = networkSize
nwLength = 10
nwLengthSD = 0
rcMean = 10
rcSD = 0

# setting simulation properties
#removalMethodsList = ['rand', 'highBC', 'lowBC']
#removalMethodsList = ['highBC', 'highBCRecalculated']
removalMethodsList = ['lowBC', 'lowBCRecalculated']
#removalMethodsList = ['lowPC', 'lowPCRecalculated']

# choosing which bits to execute
shouldGenerateNetworks = False
shouldMakePlot = True

numSamples = 6

#making sure that the iterations get saved in the correct folder
# also giving them descriptive names
prefix = '/Users/adamtrebach/Documents/Research/TPV/Nanowires/results/bc_paper_node_removal_test_bc_recalculation_during_evolution_03/'

# defining function to get a list of non electrode nodes sorted by an attribute in INCREASING order of that attribute
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

def clean_network(G):
	nodesToRemove = [node for node in G if G.degree[node] in [0,1]]
	while nodesToRemove:
		G.remove_nodes_from(nodesToRemove)
		nodesToRemove = [node for node in G if G.degree[node] in [0,1]]

netlistName = 'bc_paper_node_removal_test_bc_recalculation_during_evolution_netlist'


percMultiplesToSample = [1.2, 1.5, 1.8]
percMultiplesToSample = [1.5]
for percMultiple in percMultiplesToSample:
	if shouldGenerateNetworks:
		# defining function to get fraction of the percolation threshold the network is currently at
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
		
		
		for trialNum in range(numSamples):
			print('Beginning trial', trialNum)
			pickleName = prefix + '_'.join([str(width) + 'x' + str(height),
							'percMult' + str(percMultiple),
							'nwLen' + str(nwLength),
							'Rc' + str(rcMean),
							'trial' + str(trialNum) + '.p'])
			#when going through the iterations, we will add suffixes to indicate which removal method and which trial
			g = nwm.NanowireMesh(width = width,
					height = height,
					percMultiple = percMultiple,
					nwLength = nwLength,
					nwLengthSD = 0,
					rcMean = rcMean,
					rcSD = rcSD,
					buffer = 1)
			
			#get ordered list of betweenness centralities
			print('Calculating betweenness centralities')
			bc = nx.betweenness_centrality(g, weight = 'resistance')
			nx.set_node_attributes(g, bc, 'betweenness_centrality')
			sortedNodes = sort_nodes_by_attribute(g, attribute = 'betweenness_centrality', includeElectrodes = False)
# commented out but left for possible future reference
#			del bc[g.topElectrode]
#			del bc[g.bottomElectrode]
#			# make the sorted list
#			bcKeys = np.array(list(bc.keys()))
#			bcVals = list(bc.values())
#			sortedInds = np.argsort(bcVals)
#			bcKeys = bcKeys[sortedInds]
#			bcKeys = list(bcKeys)
#			sortedNodes = deepcopy(bcKeys)
			#sortedNodes is the nodes sorted in descending order of betweenness centrality 
		
			# calculating the initial resistance of the network
			g.to_netlist(netlistName = netlistName, voltage = 1)
			try:
				#solve circuit and update
				command = ['xyce', netlistName]
				subprocess.check_call(command,
							stdout = subprocess.DEVNULL)
				g.update_with_xyce_output(inFile = netlistName + '.csv', disableTQDM = True)
			except subprocess.CalledProcessError as err:
				print('SHIT GOT FUCKED')
				print(err.output)
				
			# calculate percolation centralities which can only be done after calculating network voltages
			print('Calculating percolation centralities')
			pc = nx.percolation_centrality(g, attribute = 'voltage', weight = 'resistance')
			nx.set_node_attributes(g, pc, 'percolation_centrality')
	
			for removalMethod in removalMethodsList:
	
				# make sure that we have a complete list of nodes sorted by betweenness centrality
				#( this list gets modified every time we use it to remove nodes so we have to reset it)
				if removalMethod in ['lowBC', 'highBC']:
					sortedNodes =  sort_nodes_by_attribute(g, attribute = 'betweenness_centrality', includeElectrodes = False)
				elif removalMethod in ['lowPC', 'lowPCRecalculated']:
					sortedNodes =  sort_nodes_by_attribute(g, attribute = 'percolation_centrality', includeElectrodes = False)


				h = deepcopy(g)
				
				# setting nodes to remove per step to be different for the different approaches
				nodesToRemovePerStep = 30 if 'low' in removalMethod else 5

				percMultipleList =[get_perc_multiple(h)]
				#
				resistance = [h.sheetResistance]
			
				newPickleEnding = removalMethod + '_iter' + '0'.zfill(5) + '.p'
				thisPickleName = pickleName.replace('.p', newPickleEnding)
				h.to_pickle(outFileName = thisPickleName)
	
				# flag indicating we are not on the last node to be removed yet?
				onLastNode = False
				n = 0
				while not onLastNode:
					n += 1
					print('method:', removalMethod, '\niteration:', n + 1)
					
					# selecting the node to remove based on the method
					if removalMethod == 'rand':
						nodesWithoutElectrodes = set(h.nodes) - {h.topElectrode, h.bottomElectrode}
						nodesToRemove = random.sample(nodesWithoutElectrodes, nodesToRemovePerStep)
					elif removalMethod == 'highBC':
						nodesToRemove = sortedNodes[-nodesToRemovePerStep:]
						del sortedNodes[-nodesToRemovePerStep:]
					elif removalMethod == 'highBCRecalculated':
						#recalculate BC every 20 nodes
						if n % 20 == 0:
							print('Recalculating BC for method', removalMethod)
							bc = nx.betweenness_centrality(h, weight = 'resistance')
							nx.set_node_attributes(h, bc, 'betweenness_centrality')
						# get a list of nodes sorted by BC
						sortedNodes = sort_nodes_by_attribute(h, attribute = 'betweenness_centrality', includeElectrodes = False)
						nodesToRemove = sortedNodes[-nodesToRemovePerStep:]
						del sortedNodes[-nodesToRemovePerStep:]
					elif removalMethod == 'lowBC':
						nodesToRemove = sortedNodes[:nodesToRemovePerStep]
						del sortedNodes[:nodesToRemovePerStep]
					elif removalMethod == 'lowBCRecalculated':
						#recalculate BC every 20 nodes
						if n % 20 == 0:
							print('Recalculating BC for method', removalMethod)
							bc = nx.betweenness_centrality(h, weight = 'resistance')
							nx.set_node_attributes(h, bc, 'betweenness_centrality')
						# get a list of nodes sorted by BC
						sortedNodes = sort_nodes_by_attribute(h, attribute = 'betweenness_centrality', includeElectrodes = False)
						nodesToRemove = sortedNodes[:nodesToRemovePerStep]
						del sortedNodes[:nodesToRemovePerStep]
					elif 'lowPC' in removalMethod:
						if 'Recalculated' in removalMethod:
							if n % 20 == 0:
								print('Recalculating PC for method', removalMethod)
								pc = nx.percolation_centrality(h, attribute = 'voltage', weight = 'resistance')
								nx.set_node_attributes(h, pc, 'percolation_centrality')
								sortedNodes = sort_nodes_by_attribute(h, attribute = 'percolation_centrality',
													includeElectrodes = False)
						nodesToRemove = sortedNodes[:nodesToRemovePerStep]
						del sortedNodes[:nodesToRemovePerStep]
					# removing the nodes one by one and making sure network is percolating after each removal
					for node in nodesToRemove:
						removedNodeData, removedEdgeData = h.pop_nodes_with_edges(node)
#						removedNodeData = {node : h.nodes[node]}
#						removedEdgeData = [(edge[0], edge[1], h.edges[edge]) for edge in h.edges(node)]
#						h.remove_node(node)
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
					clean_network(g)
					h.find_percolating_cluster()
					# writing to circuit and solving
					h.to_netlist(netlistName = netlistName, voltage = 1)
					command = ['xyce', netlistName]
					try:
						# solve circuit and update
						subprocess.check_call(command, stdout = subprocess.DEVNULL)
						h.update_with_xyce_output(inFile = netlistName + '.csv', disableTQDM = True)
						
						# save the sheet resistance data
#						clean_network(h)
						percMultipleList.append(get_perc_multiple(h))
						resistance.append(h.sheetResistance)
		
						# save the pickle of this network every 30 samples or if we're on the last node
						iterNumber = str(n + 1)
						if (n + 1) % 30 == 0 or onLastNode:
							newPickleEnding = removalMethod + '_iter' + iterNumber.zfill(5) + '.p'
							thisPickleName = pickleName.replace('.p', newPickleEnding)
							h.to_pickle(outFileName = thisPickleName)
							if onLastNode:
								# if the network is about to stop percolating, we are done. We want to record the last snapshot here and move on.
								print('Circuit is one node away from failure. Finished with simulations for method', removalMethod)
								# because we have set onLastNode = True earlier in the code, this is the last thing 
								# that will be executed in the while loop
					except subprocess.CalledProcessError as err:
						print('SHIT GOT FUCKED')
						print(err.output)
			
				# saving the data as a csv
				data = [(percMultipleList[i] , resistance[i]) for i in range(len(percMultipleList))]
				data = np.array(data)
				fileEnding = removalMethod + '_data'
				fileName = pickleName.replace('.p', fileEnding)
				np.savetxt(fileName, data, fmt = '%f', delimiter = ' ')
			
	if shouldMakePlot:
		prefix = prefix + '_'.join([str(width) + 'x' + str(height),
							'percMult' + str(percMultiple),
							'nwLen' + str(nwLength),
							'Rc' + str(rcMean)])
	
		plotScriptFileName = '/Users/adamtrebach/Documents/Research/TPV/Nanowires/results/bc_paper_node_removal_test_bc_recalculation_during_evolution_04/evolution_under_node_removal_test_bc_recalculation_multiplot_' + str(width) + 'x' + str(height) + 'percMult' + str(percMultiple) + 'plot_script.txt'
		plotScript = ['set term png size 8000, 5000 linewidth 7',
				'set output \"' + plotScriptFileName.replace('plot_script.txt', '.png\"'),
				'set multiplot layout 3,2 scale 1,1 title \"Testing PC Recalculation with Network Cleaning\" font \",100\" offset 0,-15',
				'set tmargin 20',
				'set lmargin 30',
				'set rmargin 30',
				'set xlabel \"Network Density as Multiple of Percolation Threshold\" font \",50\" offset 0,-6',
				'set ylabel \"Series resistance of network in ohms\" font \",50\" offset -9,0',
				'set xrange [2:0]',
				'set yrange [0:600]',
				'set xtics font \",40\" offset 0,-3',
				'set ytics font \",40\"',
				'set key font \",50\"']
	
		# setting line colors so I can set matching point colors later on 
		lineColorDict = {'rand' : '\"red\"',
				'highBC' : '\"forest-green\"',
				'lowBC' : '\"blue\"',
				'highBCRecalculated' : '\"blue\"',
				'lowBCRecalculated' : '\"orange\"',
				'lowPC' : '\"blue\"',
				'lowPCRecalculated' : '\"orange\"'}

		runDataDf = pd.DataFrame(columns = ['removalMethod', 'RDoublePercMult', 'EndPercMult'])
	
		for trialNumber in range(numSamples):
			plotScript.append('set title \"Trial ' + str(trialNumber) + '\" font \",75\" offset 0,-3')
			subplotCommandsList = []
			for n, removalMethod in enumerate(removalMethodsList):
				fileName = prefix + '_trial' + str(int(trialNumber)) + removalMethod + '_data'
				command = ('plot \"' if n == 0 else '\"') + fileName + '\" using 1:2 with lines title \"' + removalMethod + ' nodes removed first\" linecolor ' + lineColorDict[removalMethod]
				subplotCommandsList.append(command)
	
				# adding the markers for the final density at the end of the plot
				df = pd.read_csv(fileName, sep = ' ', names = ['percMultiple', 'resistance'])
				minPercMultiple = df['percMultiple'].min()	
				minMarker = np.array([(minPercMultiple, 10)])
				minMarkerFileName = fileName + 'min_percMultiple'
				np.savetxt(minMarkerFileName, minMarker, fmt = '%f', delimiter = ' ')
	
				command = '\"' + minMarkerFileName + '\" using 1:2 with points pointsize 10 pointtype 2 linecolor ' + lineColorDict[removalMethod] + ' notitle'
				subplotCommandsList.append(command)

				# calculating the density at which resistance doubles
				percToResFunc = interpolate.interp1d(df['percMultiple'], df['resistance'])
				funcToOptimize = lambda x : percToResFunc(x) - 2 * df['resistance'].min()
				rootResults = optimize.root_scalar(funcToOptimize,
									bracket = (df['percMultiple'].min(), df['percMultiple'].max()))
				if rootResults.converged:
					runDataDf = runDataDf.append(
									dict(RDoublePercMult = rootResults.root, 
										EndPercMult = df['percMultiple'].min(),
										removalMethod = removalMethod),
									ignore_index = True)
		
		
			subplotCommand = ', '.join(subplotCommandsList)
			plotScript.append(subplotCommand)
	#	plotScript = plotScript + ['unset multiplot', 'set bmargin at screen 0.2']
		plotScriptString = '\n'.join(plotScript)
		plotScriptFile = open(plotScriptFileName , 'w')
		plotScriptFile.write(plotScriptString)
		plotScriptFile.close()
	
		try:
			subprocess.check_call(['gnuplot', plotScriptFileName])
		except subprocess.CalledProcessError as err:
			print('Error in Gnuplot Script')
			print(err.output)
			


		# getting average and sterr from the summary Df
		sem = stats.sem
		meanDf = runDataDf.groupby(by = ['removalMethod']).mean()
		stdDf = runDataDf.groupby(by = ['removalMethod']).sem()
		stdNameMapping = {name : name + 'SEM' for name in stdDf.columns}
		stdDf = stdDf.rename(columns = stdNameMapping)
		resultsDf = pd.concat([meanDf, stdDf], axis = 1, join = 'inner')
		resultsDf.to_csv(prefix + 'summary_by_removal_method.csv', index = True)

