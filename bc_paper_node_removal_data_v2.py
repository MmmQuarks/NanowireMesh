import random
from copy import deepcopy
import subprocess
import NanowireMesh as nwm
import os
import networkx as nx
import numpy as np
import itertools
import ParallelCentrality as PC
from scipy.special import factorial as factorial
from scipy import interpolate, optimize, stats
import pandas as pd
import time
import pdb
import sys
#import gnuplotlib as gp

# simulation options
shouldGenerateNetworks = True
shouldGenerateSimulatedResistances = True
shouldEvolveNetworks = True
shouldCalculateInitialCentralities = True
shouldMakePlot = False # makes plots of evolution of each trial using multiplot
shouldMakeAverageFailurePlot = False
shouldNormalizeResistances = False
numSamples = 10 #note: if testing, numSamples must be >= 2 because we calculat std at the end
xycePath = '/home/amwt/TPV/2DNanowires/installations/build_Xyce-6.11_parallel_openMPI/bin/Xyce'
netlistSuffix = '_'.join(sys.argv[1:]) # making a netlist suffix so different instances of this script can run concurrently
existingNetworkFolder = '/home/amwt/TPV/2DNanowires/Data/evolution_under_node_removal/' # name directory where existing networks can be found
resultsFolder = '/home/amwt/TPV/2DNanowires/Data/evolution_under_node_removal/bc_paper_node_removal_figure_v2_0001/'
simulatedResistancesFolder = '/Users/adamtrebach/Documents/Research/TPV/Nanowires/Results/bc_paper_simulated_resistances_04/'
removalMethodsToCalculate = ['rand',
				'lowBC',
				'lowEC',
				'lowPC',
				'lowCW',
				'lowPW',
				'lowAngle'
				]

#setting network properties
params = {'width' : 100,
		'height' : 100,
		'nwLength' : 10,
		'nwLengthSD' : 0,
		'rcMean' : 10,
		'rcSD' : 0,
		'buffer' : 1,
		'percMultiple' : None}

# setting simulation properties
#removalMethodsList = ['rand'] 
#colorsList = ['black']
#
#removalMethodsList += ['lowBC', 'highBC']
#colorsList += ['light-red', 'dark-red']
#
#removalMethodsList += ['lowEC', 'highEC']
#colorsList += ['web-green', 'dark-green']
#
#removalMethodsList += ['lowPC', 'highPC'] 
#colorsList += ['web-blue', 'dark-blue']
#
#removalMethodsList += ['lowDownPotentialPC',  'highDownPotentialPC']
#colorsList += ['sienna1', 'dark-orange']
#
#removalMethodsList += ['lowCW', 'highCW']
#colorsList += ['orchid', 'orchid4']
#
#removalMethodsList += ['lowPW', 'highPW']
#colorsList += ['khaki1', 'gold']
#
#removalMethodsList += ['lowAngle', 'highAngle']
#colorsList += ['light-pink', 'dark-pink']
#
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
			'lowPW' : dict(color = 'khak:i1', centralityName = 'power_weighted_centrality'),
			'highPW' : dict(color = 'gold', centralityName = 'power_weighted_centrality'),
			'lowAngle' : dict(color = 'light-pink', centralityName = 'absolute_angle'),
			'highAngle' : dict(color = 'dark-pink', centralityName = 'absolute_angle')}
pdb.set_trace()
		


#mapFromRemovalMethodToCentralityName = {'highBC' : 'betweenness_centrality',
#					'lowBC' : 'betweenness_centrality',
#					'highEC' : 'electrode_centrality',
#					'lowEC' : 'electrode_centrality',
#					'lowPC' : 'percolation_centrality',
#					'highPC' : 'percolation_centrality',
#					'lowDownPotentialPC' : 'down_potential_percolation_centrality',
#					'highDownPotentialPC' : 'down_potential_percolation_centrality',
#					'lowCW' : 'current_weighted_centrality',
#					'highCW' : 'current_weighted_centrality',
#					'lowPW' : 'power_weighted_centrality',
#					'highPW' : 'power_weighted_centrality',
#					'lowAngle' : 'absolute_angle',
#					'highAngle' : 'absolute_angle'}
#
#colorsList = ['black', 'green', 'red', 'blue', 'gold', 'dark-magenta',
#		 'dark-green', 'gray-40', 'sienna1', 'royalblue',
#		'olive', 'aquamarine']
#if len(colorsList) == len(removalMethodsList):
#	colorsDict = { removalMethodsList[n] : colorsList[n] for n in range(len(removalMethodsList))}
#else:
#	raise Exception('Number of colors does not match number of removal methods. Plots will be fucky.')	
#		

#------------------------------------end settings -------------------------------

try:
	# try to make the folder
	os.mkdir(resultsFolder)
	print('Created folder', resultsFolder)
except FileExistsError:
	# if the folder already exists, check to make sure it's empty
	resultsFolderContents =	os.listdir(resultsFolder)
	if len(resultsFolderContents) > 0:
		# if results folder has contents, raise exception
#		msg = 'Results folder ' + resultsFolder + ' is not empty. Either change the results folder or empty the current one.'
		print(' Results folder is not empty. Proceed with caution.')			
#		raise Exception(msg)
#	else:
#		print('Using empty folder', resultsFolder)
#	



def sort_nodes_by_attribute(G,
				attribute,
				includeElectrodes = False):
	# this function returns the non-electrode nodes from the graph G sorted in decreasing order of some attribute
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


percMultiplesToSample = [float(val) for val in sys.argv if val != sys.argv[0]]
for percMultiple in percMultiplesToSample:
	if shouldEvolveNetworks:
		for trialNum in range(numSamples):
			print('Beginning trial', trialNum)
			#when going through the iterations, we will add suffixes to indicate which removal method and which trial
	
			# this sets the pickle name for the networks we are going to be saving
			pickleName = resultsFolder + '_'.join([str(params['width']) + 'x' + str(params['height']),
							'percMult' + str(percMultiple),
							'nwLen' + str(params['nwLength']),
							'Rc' + str(params['rcMean']),
							'trial' + str(trialNum) + '.p'])
			
			if shouldGenerateNetworks:
				# generating a new network if we're doing so
				params['percMultiple'] = percMultiple
				g = nwm.NanowireMesh(**params)
				g.to_pickle(outFileName = pickleName)
				 
			else:
				#opening the existing networks that we're already using
				# assuming everything is named in the style of pickleName
				g = nwm.NanowireMesh(inPickle = pickleName)

				# this should not longer be necessary because the centralities will be saved with the first pickle
				# now
				# copying calculated centralities from  other pickles (iter00000.p)
#				pickleWithCentralities = pickleName.replace('.p', 'rand_iter00000.p')
#				g2 = nwm.NanowireMesh(inPickle = pickleWithCentralities)
#				if set(g.nodes) == set(g2.nodes):
#					centralities = ['betweenness_centrality',
#							'electrode_centrality',
#							'percolation_centrality',
#							'down_potential_percolation_centrality',
#							'current_weighted_centrality',
#							'power_weighted_centrality',
#							'absolute_angle']
#					for elem in centralities:
#						cDict = nx.get_node_attributes(g2, elem)
#						nx.set_node_attributes(g, cDict, elem)

	
			# calculating the initial resistance of the network
			netlistName = 'netlist' + netlistSuffix
			# cleaning network
			#clean_network(g)
			g.solve_circuit_using_xyce(xycePath = xycePath, netlistName = netlistName, voltage = 1)
			
			if shouldCalculateInitialCentralities:
				print('Calculating initial centralities')
				# betweenness centrality
				print('Starting initial betweenness centrality calc.')
				c = nx.betweenness_centrality(g, weight = 'resistance')
				nx.set_node_attributes(g, c, 'betweenness_centrality')
				
				# electrode centrality
				print('Starting initial electrode centrality calc.')
				c = nwm.electrode_centrality(g, potential = 'voltage', weight = 'resistance', normalized = True)
				nx.set_node_attributes(g, c, 'electrode_centrality')	
		
				# percolation centrality
				print('Starting initial percolation centrality calc.')
				c = nx.percolation_centrality(g, attribute = 'voltage', weight = 'resistance')
				nx.set_node_attributes(g, c, 'percolation_centrality')
		
				# down potential  percolation centrality
#				print('Starting initial down potential percolation centrality calc.')
#				c = nwm.percolation_centrality(g, attribute = 'voltage', weight = 'resistance', downPotentialOnly = True)
#				nx.set_node_attributes(g, c, 'down_potential_percolation_centrality')
				
				# current weighted centrality
				print('Starting initial current weighted centrality calculation')
				c = PC.current_weighted_centrality(g)
				nx.set_node_attributes(g, c, 'current_weighted_centrality')
	
				# power weighted centrality
				print('Starting initial power weighted centrality calculation')
				c = PC.power_weighted_centrality(g)
				nx.set_node_attributes(g, c, 'power_weighted_centrality')
	
				# angle centrality
				print('Calculating absolute value of angles')
				c = {node : abs(g.nodes[node]['angle']) for node in g}
				nx.set_node_attributes(g, c, 'absolute_angle')

				# save network with centralities calculated.
				g.to_pickle(outFileName = pickleName)
	
	
			for removalMethod in removalMethodsToCalculate:
	
				# make sure that we have a complete list of nodes sorted by betweenness centrality
				#( this list gets modified every time we use it to remove nodes so we have to reset it)
				if removalMethod == 'rand':
					pass
				else:
					centralityName = removalMethodParams[removalMethod]['centralityName']
					centralitySortedNodes = sort_nodes_by_attribute(g,
											attribute = centralityName,
											includeElectrodes = False)
									
				h = deepcopy(g)
				
				# setting nodes to remove per step to be different for the different approaches
				# if we're starting with low centrality, we want to remove more nodes at a time
				nodesToRemovePerStep = 20 if 'low' in removalMethod  else 5
	
			
				percMultipleList =[get_perc_multiple(h)]
				resistance = [h.sheetResistance]
			
				newPickleEnding = removalMethod + '_iter' + '0'.zfill(5) + '.p'
				thisPickleName = pickleName.replace('.p', newPickleEnding)
				h.to_pickle(outFileName = thisPickleName)
	
				# flag indicating we are not on the last node to be removed yet
				onLastNode = False
				n = 0
				while not onLastNode:
					n += 1
					print('method:', removalMethod, '\niteration:', n)
	
					#re-calculate centrality as network evolves. This is unnecessary on the first iteration because 
					# we calculated all centralities above before the removalMethod loop
					# we want to skip this for now so we are not recalculating
					if n != 1 and removalMethod != 'rand' and False:
						centralityName = mapFromRemovalMethodToCentralityName[removalMethod]
						if centralityName == 'betweenness_centrality':	
							c = nx.betweenness_centrality(h, weight = 'resistance')
						elif centralityName == 'electrode_centrality':	
							c = nwm.electrode_centrality(h, potential = 'voltage', weight = 'resistance', normalized = True)
						elif centralityName == 'percolation_centrality':
							c = nwm.percolation_centrality(h, attribute = 'voltage', weight = 'resistance', downPotentialOnly = False)
						elif centralityName == 'down_potential_percolation_centrality':
							c = nwm.percolation_centrality(h, attribute = 'voltage', weight = 'resistance', downPotentialOnly = True)
						# assigning centralities to nodes in graph
						nx.set_node_attributes(g, c, centralityName)
					
					# selecting the nodes to remove based on the method
					if removalMethod == 'rand':
						nodesWithoutElectrodes = set(h.nodes) - {h.topElectrode, h.bottomElectrode}
						nodesToRemove = random.sample(nodesWithoutElectrodes, nodesToRemovePerStep)
					else:
						# tried to DRY here but may be confusing
						centralityName = removalMethodParams[removalMethod]['centralityName']
						centralitySortedNodes = sort_nodes_by_attribute(h,
												attribute = centralityName,
												includeElectrodes = False)
						if 'low' in removalMethod:
							nodesToRemove = centralitySortedNodes[0:nodesToRemovePerStep]
						elif 'high' in removalMethod:
							nodesToRemove = centralitySortedNodes[-nodesToRemovePerStep:]
							# we reverse this so the nodes are in order from highest centrality to lowest
							nodesToRemove.reverse()
						
	
	
					# removing the nodes one by one and making sure network is percolating after each removal
					while nodesToRemove:
						node = nodesToRemove.pop(0)

						# making sure that we remove the entire wire of which this node is a segment
						allSegmentNodes = get_wire_segments(h, node)

						# removing the segment nodes from the nodesToRemove list (if they are in there)
						# this makes sure that we never attempt to remove the same node twice
						for n in allSegmentNodes - {node}: 
							# we don't need to remove node from nodesToRemove because it has already been popped
							try:
								nodesToRemove.remove(n)
							except ValueError:
								# if we get a value error it means that this segment node=
								# was not in the list of nodesToRemove in this iteration
								# we should ignore this and check the next node
								pass
						removedNodeData, removedEdgeData = h.pop_nodes_with_edges(allSegmentNodes)
						try:
							assert h.topElectrode in h
						except AssertionError:
							print('Assetion error')
							pdb.set_trace()
						h.isPercolating = nx.has_path(h, h.topElectrode, h.bottomElectrode)
						# if the network is percolating, keep removing nodes
						# if the network is no longer percolating, we should
							# 1. set h = networkBeforeNodeRemoval 
							# 2. clean the network one last time
							# 3. solve this last circuit
							# 4. save network
							# 5. stop removing nodes
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
						# solve circuit and update
						h.solve_circuit_using_xyce(xycePath = xycePath, netlistName = netlistName, voltage = 1)

						# save the sheet resistance data
						percMultipleList.append(get_perc_multiple(h))
						resistance.append(h.sheetResistance)
		
						# save the pickle of this network every 10 samples or if we're on the last node
						iterNumber = str(n)
						if n % 10 == 0 or onLastNode:
							newPickleEnding = removalMethod + '_iter' + iterNumber.zfill(5) + '.p'
							thisPickleName = pickleName.replace('.p', newPickleEnding)
							h.to_pickle(outFileName = thisPickleName)
							if onLastNode:
								# if the network is about to stop percolating, we are done. 
								# We want to record the last snapshot here and move on.
								print('Circuit is one node away from failure. Finished with simulations for method', removalMethod)
								# because we have set onLastNode = True earlier in the code, this is the last thing 
								# that will be executed in the while loop
					except subprocess.CalledProcessError as err:
						print('SHIT GOT FUCKED')
						print(err.output)
			
				# saving the data as a csv
				data = [(percMultipleList[i] , resistance[i]) for i in range(len(percMultipleList))]
				data = np.array(data)
				fileEnding = removalMethod + '_data.txt'
				fileName = pickleName.replace('.p', fileEnding)
				np.savetxt(fileName, data, fmt = '%f', delimiter = ' ')
				
	# making plots of the evolution of each trial using multiplot
	if shouldMakePlot:
		plotScriptFileName = resultsFolder  + str(params['width']) + 'x' + str(params['height']) + 'percMult' + str(percMultiple) + 'plot_script.txt'

		plotScript = ['set term png size 8000, 6000 linewidth 7',
				'set output \"' + plotScriptFileName.replace('plot_script.txt', '.png\"'),
				'set multiplot layout 4,3 scale 1,1 title \"Network Evolution Under Node Removal Starting at Percolation Multiple ' + str(percMultiple) + '\" font \",100\" offset 0,-15',
				'unset key',
				'set tmargin 20',
				'set lmargin 30',
				'set rmargin 60',
				'set xlabel \"Network Density as Multiple of Percolation Threshold\" font \",50\" offset 0,-6',
				'set ylabel \"Series resistance of network in ohms\" font \",50\" offset -9,0',
				'set xrange [' + str(1.15 * percMultiple) + ':0]',
				'set yrange [0:30]',
				'set xtics font \",40\" offset 0,-3',
				'set ytics font \",40\"',
				'#set key inside center top horizontal font \",50\"']
	
	
		for trialNumber in range(numSamples):
			plotScript.append('set title \"Trial ' + str(trialNumber) + '\" font \",75\" offset 0,-3')
			subplotCommandsList = []
			# if we are on the last trial, add the legend
			if trialNumber == numSamples - 1:
				plotScript.append('set key at graph 1.5,screen 0.1 center bottom font \",50\"')


			for n, removalMethod in enumerate(removalMethodsToCalculate):
				fileName = ''.join([resultsFolder,
							str(params['width']) + 'x' + str(params['height']),
							'_percMult' + str(percMultiple),
							'_nwLen' + str(params['nwLength']),
							'_Rc' + str(params['rcMean']),
							'_trial' + str(int(trialNumber)),
							removalMethod + '_data.txt'])
				if shouldNormalizeResistances:
					nData = np.loadtxt(fileName)
					nData[:,1] = nData[:,1] / nData[0,1]
					fileName = fileName.replace('.txt.', '_normalized.txt')
					np.savetxt(fileName, nData)
				command = ''.join(['plot ' if n == 0 else '',
							'\"' + fileName + '\"',
							'using 1:2 with lines ',
							'title \"' + removalMethod + ' nodes removed first\" ',
							'lc \"' + colorsDict[removalMethod] + '\"'])
				subplotCommandsList.append(command)
				# adding a key only if we're on the last one
				# adding the markers for the final density at the end of the plot
				df = pd.read_csv(fileName, sep = ' ', names = ['percMultiple', 'resistance'])
				minPercMultiple = df['percMultiple'].min()	
				maxResistance = df['resistance'].max()
#				minMarker = np.array([(minPercMultiple, 10)])
				minMarker = np.array([(minPercMultiple, 0), (minPercMultiple, maxResistance)])
				minMarkerFileName = fileName + 'min_percMultiple.txt'
				np.savetxt(minMarkerFileName, minMarker, fmt = '%f', delimiter = ' ')
	
				command = ''.join(['\"' + minMarkerFileName + '\"',
							'using 1:2 with lines lw 0.7 ',
							'lc \"' + colorsDict[removalMethod] + '\" ',
							'dt \"   .   \" notitle'])
				subplotCommandsList.append(command)
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
			


# make plot of average densities at failure for different removal methods and different starting densities
if shouldMakeAverageFailurePlot:
	dfColumns = ['percMultipleGroup', 'startPercMultiple', 'trial', 'removalMethod', 'endPercMultiple', 
			'halfwayPercMultiple', 'halfwayNormalizedPercMultiple', 'halfwayResistance', 'halfwayNormalizedResistance']
	df = pd.DataFrame(data = None, index = None, columns = dfColumns)
	relDf = pd.DataFrame(data = None, index = None, columns = dfColumns)

	# making container for the high density data to be compared to 
	# simulated networks at various densities
	hdDf = pd.DataFrame(columns = ['percMultiple', 'resistance', 'removalMethod'])

	# iterating through data to make data for average failure
	# and high density comparison to simulated resistances
	for percMultiple in percMultiplesToSample:
		for trialNumber in range(numSamples):
			for removalMethod in removalMethodsToCalculate:
				# calculating summary data like doubling point, triplign point, etc
				fileName = ''.join([resultsFolder,
							str(params['width']) + 'x' + str(params['height']),
							'_percMult' + str(percMultiple),
							'_nwLen' + str(params['nwLength']),
							'_Rc' + str(params['rcMean']),
							'_trial' + str(int(trialNumber)),
							removalMethod + '_data.txt'])
				
				#opening the file in pandas
				fileDf = pd.read_csv(fileName, sep = ' ', names = ['percMultiple', 'resistance'])
				thisData = {'percMultipleGroup' : percMultiple,
						'startPercMultiple' : fileDf['percMultiple'].max(),
						'trial' : trialNumber, 
						'removalMethod' : removalMethod,
						'endPercMultiple' : fileDf['percMultiple'].min()}

				thisData['halfwayPercMultiple'] = np.average([thisData['startPercMultiple'], thisData['endPercMultiple']])
				thisData['halfwayNormalizedPercMultiple'] = thisData['halfwayPercMultiple'] / thisData['startPercMultiple']
				# making interpolting function to get f(percMult) = resistance
				percToResFunc = interpolate.interp1d(fileDf['percMultiple'], fileDf['resistance'])
				# have to flip the data b/c interp requires arguments to be ascending
				halfwayResistance = percToResFunc(thisData['halfwayPercMultiple'])
			#	np.interp(thisData['halfwayPercMultiple'], 
			#					np.flip(fileDf['percMultiple']),
			#					np.flip(fileDf['resistance']))
				thisData['halfwayResistance'] = halfwayResistance
				thisData['halfwayNormalizedResistance'] = halfwayResistance  / fileDf['resistance'].min()

				# using interpolated function to find doubling point, tripling point, etc
				for factor in range(2,6):
					func = lambda x : percToResFunc(x) - factor * fileDf['resistance'].min()
					try:
						rootResults = optimize.root_scalar(func,
									bracket = (fileDf['percMultiple'].min(), fileDf['percMultiple'].max()))
						if rootResults.converged:
							thisData['ResistanceX' + str(factor) + 'percMult'] = rootResults.root
					except ValueError:
						# sometimes the resistance never increase by eg 5x so we get a value error
						# this is not a problem
						pass
				df = df.append(thisData, ignore_index = True)
				
				# making the version with the data I need 

				#filling out the high density data
				if percMultiple == 2.1:
#					if removalMethod == 'lowBC':
					tempDf = fileDf.loc[fileDf['percMultiple'] > 7 / 5.63726, :]
					tempDf.insert(len(tempDf.columns), 'removalMethod', removalMethod) 
					hdDf = hdDf.append(tempDf)

	# making summary data (tripling point, etc)	
	#for percMultiple in percMultiplesToSample:
	#thisDfRows = df['percMultipleGroup'] == percMultiple
	#thisDf = df[thisDfRows]
	#thisDf = df[['removalMethod', 'startPercMultiple', 'endPercMultiple']]
	meanDf = df.groupby(by = ['removalMethod', 'percMultipleGroup']).mean()
	
	# get std
	sem = stats.sem
	stdDf = df.groupby(by = ['removalMethod', 'percMultipleGroup']).sem()
	nameMapping = {name : name + 'SEM' for name in stdDf.columns}
	stdDf = stdDf.rename(columns = nameMapping)

	
	resultsDf = pd.concat([meanDf, stdDf], axis = 1, join = 'inner')
	pdb.set_trace()
	cols = resultsDf.columns.tolist()
	cols = sorted(cols)
	resultsDf = resultsDf[cols]
	resultsDf.to_csv(resultsFolder + 'summary_by_method_and_starting_density.csv', index = True, sep = ',')
	

	# now just grouping by method and ignoring the different starting densities
	meanDfMethods = df.groupby(by = ['removalMethod']).mean()
	stdDfMethods = df.groupby(by = ['removalMethod']).sem()
	nameMapping = {name : name + 'SEM' for name in stdDfMethods.columns}
	stdDfMethods = stdDfMethods.rename(columns = nameMapping)
	resultsDfMethods = pd.concat([meanDfMethods, stdDfMethods], axis = 1, join = 'inner')
	cols = resultsDfMethods.columns.tolist()
	cols = sorted(cols)
	resultsDfMethods = resultsDfMethods[cols]
	print(resultsDfMethods)
	resultsDfMethods.to_csv(resultsFolder + 'summary_by_method.csv', index = True, sep = ',')


	# generating simulated networks to probe "experimental" values of resistance
	# at given densities
	print("Generating 'experimental' data to compare to network evolution.")
	percMultiplesToSimulate = np.linspace(1.0, 2.4, 20)
	trialsPerPercMultiple = 30

	if shouldGenerateSimulatedResistances:
		df = pd.DataFrame(columns = ['percMultiple' , 'resistance'])
		for percMultiple in percMultiplesToSimulate:
			for trialNum in range(trialsPerPercMultiple):
				params['percMultiple'] = percMultiple
				g = nwm.NanowireMesh(**params)
				#clean_network(g)
				pickleName = ''.join([resultsFolder,
							'percMultiple' + str(round(percMultiple, 3)),
							'_trial' + str(trialNum).zfill(2)])
				g.to_pickle(outFileName = pickleName)

				# solving circuit
				try:
					g.solve_circuit_using_xyce(xycePath = xycePath,
									netlistName = 'netlist',
									voltage = 1)
					data = dict(percMultiple = g.get_perc_multiple(),
							resistance = g.sheetResistance)
					df = df.append(data, ignore_index = True)
				except:
					print('Error in solving circuit')
					pdb.set_trace()
		df.to_csv(path_or_buf = simulatedResistancesFolder + 'simulated_resistances.csv')
	else:
		df = pd.read_csv(filepath_or_buffer = simulatedResistancesFolder + 'simulated_resistances.csv')

	# calculating binned statistics
	percMultipleMeans, percMultipleBinEdges, percMultipleBinNumbers = stats.binned_statistic(df['percMultiple'],
													df['percMultiple'],
													statistic = 'mean',
													bins = 25)
	percMultipleStds = stats.binned_statistic(df['percMultiple'],
							df['percMultiple'],
							statistic = sem,
							bins = 25)[0]
	resistanceMeans = stats.binned_statistic(df['percMultiple'],
							df['resistance'],
							statistic = 'mean',
							bins = 25)[0]
	resistanceStds = stats.binned_statistic(df['percMultiple'],
							df['resistance'],
							statistic = sem,
							bins = 25)[0]
	resultsDict = dict(percMultipleMean = percMultipleMeans,
				percMultipleStd = percMultipleStds,
				resistanceMean = resistanceMeans,
				resistanceStd = resistanceStds)
	resultsDf = pd.DataFrame(resultsDict)
	resultsDf.to_csv(simulatedResistancesFolder + 'summary.csv', index = True)


	# calculating binned statistics for evolution of high density (pm = 2.1) networks
	# to compare to simulated networks that start at given densities
	for n, removalMethod in enumerate(removalMethodsToCalculate):
		tempDf = hdDf.loc[hdDf['removalMethod'] == removalMethod, :]
		hdMeans, hdMeanBinEdges, hdMeanBinNumbers = stats.binned_statistic(tempDf['percMultiple'],
									tempDf['resistance'],
									statistic = 'mean',
									bins = 15)
		hdSems, hdSemBinEdges, hdSemBinNumbers  =  stats.binned_statistic(tempDf['percMultiple'], 
									tempDf['resistance'],
									statistic = sem,
									bins = 15)
		percMultiples = [np.average([hdMeanBinEdges[n], hdMeanBinEdges[n+1]]) for n in range(len(hdMeanBinEdges) - 1)]
		methodDf = pd.DataFrame()
		methodDf['percMultiple'] = percMultiples
		methodDf['resistanceMean'] = hdMeans
		methodDf['resistanceSEM'] = hdSems
		csvName = '_'.join([resultsFolder + 'highDensity', removalMethod + '.csv'])
		methodDf.to_csv(path_or_buf = csvName,
				index = False)
		
		if n == 0:
			# calculating the forro data
			# function signature is forro_sheet_resistance(percMultiple, Rw, Rc)
			# first must calc_internal_resistivity(self, diam, temp)
#			nwDiam = 0.15
#			nwLength = 10
#			temp = 298.15
#			# only instantiating a network so we can use a class method
#			g = nwm.NanowireMesh(width = 50, height = 50)
#			rho = g.calc_internal_resistivity(diam = nwDiam, temp = temp)
#			Rw = rho * nwLength / (np.pi * (nwDiam / 2)**2)
#			Rc = 10
#			forroPercMults = np.linspace(min(hdMeanBinEdges), max(hdMeanBinEdges),50)
#			forroResistances = np.array([nwm.forro_sheet_resistance(percMultiple = pm, Rw = Rw, Rc = Rc) for pm in forroPercMults])
#		
#			forroDf = pd.DataFrame()
#			forroDf['percMultiple'] = forroPercMults
#			forroDf['resistance'] = forroResistances
#			forroDf.to_csv(path_or_buf = resultsFolder + 'highDensity_forro.csv')
#		
#			# making gnuplot parameters
#			forroCurveOptions = {'with' : 'lines lw 4',
#						'legend' : 'Theoretical'}
#			forroCurve = (forroDf['percMultiple'],
#					forroDf['resistance'],
#					forroCurveOptions)
			simResDf = pd.read_csv(filepath_or_buffer = '/Users/adamtrebach/Documents/Research/TPV/Nanowires/Results/bc_paper_simulated_resistances_02/summary.csv')
			simResDf = simResDf.loc[simResDf['percMultipleMean'] < max(hdMeanBinEdges), :]
			simResRegionOptions = {'with' : 'filledcurves lc \"orange\"',
						'tuplesize' : 3}
			simResRegion = (simResDf['percMultipleMean'],
					(simResDf['resistanceMean'] - simResDf['resistanceStd']) / simResDf['resistanceMean'].min(),
					(simResDf['resistanceMean'] + simResDf['resistanceStd']) / simResDf['resistanceMean'].min(),
					simResRegionOptions)

		
		#making the plots using gnuplotlib
		curveOptions = {'with' : 'yerrorbars linestyle 7 lw 3 lc \"black\"',
					'tuplesize' : 3}
		#maxResistance = max(max(forroDf['resistance']), max(methodDf['resistanceMean']))
		plotOptions = dict(xrange = '2.4:1.2',
#					yrange = '5:' + str(maxResistance + 1),
					hardcopy = csvName.replace('.csv', '_ForroComparisonPlot.png'),
					cmds = ['set lmargin 18',
						'set rmargin 10',
						'set bmargin 10',
						'set tmargin 3',
						'set yrange [0:8]',
						'set key tmargin font \",20\"',
						'set ylabel \"Resistance Relative to Start\" font \",30\" offset -4,0',
						'set xlabel \"Multiple of Percolation Threshold\" font \",30\" offset 0,-5',
						'set ytics font \",30"',
						'set xtics font \",30\" offset 0,-2'],
					terminal = 'png size 1000,1000')
		curve = (methodDf['percMultiple'],
				methodDf['resistanceMean'] / methodDf['resistanceMean'].min(),
				methodDf['resistanceSEM'] / methodDf['resistanceMean'].min(),
				curveOptions)
		gp.plot(simResRegion, curve,  **plotOptions)







	# making a different data file for each removal method
	removalMethodFileNameDict = {}
	for removalMethod in removalMethodsToCalculate:
		try:
			dfToPrint = resultsDf.loc[(removalMethod)]
		except KeyError:
			pdb.set_trace()
		fileName = ''.join([resultsFolder,
					str(params['width']) + 'x' + str(params['height']),
					'_nwLen' + str(params['nwLength']),
					'_Rc' + str(params['rcMean']),
					'_removalMethod=' + removalMethod + '.txt'])
		removalMethodFileNameDict[removalMethod] = fileName
		dfToPrint.to_csv(fileName, index = False, sep = ' ')
	
	#making plot script

	plotScriptFileName = ''.join([resultsFolder,
					str(params['width']) + 'x' + str(params['height']),
					'_nwLen' + str(params['nwLength']),
					'_Rc' + str(params['rcMean']),
					'comparing_failure_density_plot_script.txt'])

	plotScript = ['set term png size 7000, 5000 linewidth 7',
			'set output \"' + plotScriptFileName.replace('plot_script.txt', '.png\"'),
			'set  title \"Network Evolution Under Node Removal\" font \",100\" offset 0,-2',
			'set tmargin 20',
			'set lmargin 50',
			'set rmargin 270',
			'set bmargin 30',
			'set xrange [2:0]',
			'set xlabel \"Network Density as Multiple of Percolation Threshold when Percolation is Destroyed\" font \",80\" offset 0,-14',
			'set ylabel \"Network Density as Multiple of Percolation Threshold when Simulations Began\" font \",80\" offset -14,0',
			'set xtics font \",60\" offset 0,-5',
			'set ytics font \",60\" offset -1,0',
			'set key outside right center font \",70\"']

	plotCommands = []	
	for n, removalMethod in enumerate(removalMethodsToCalculate):
		fileName = removalMethodFileNameDict[removalMethod]
		command = ('plot \"' if n == 0 else '\"') + fileName + '\" using 2:1:4:3  with xyerrorbars  title \"' + removalMethod + ' nodes removed first\" lc \"' + colorsDict[removalMethod] + '\"'
		plotCommands.append(command)	

	plotCommands = ', '.join(plotCommands)
	plotScript.append(plotCommands)
	plotScriptString = '\n'.join(plotScript)
	plotScriptFile = open(plotScriptFileName , 'w')
	plotScriptFile.write(plotScriptString)
	plotScriptFile.close()

	try:
		subprocess.check_call(['gnuplot', plotScriptFileName])
	except subprocess.CalledProcessError as err:
		print('Error in Gnuplot Script')
		print(err.output)

