import NanowireMesh as nwm
import subprocess
import networkx as nx
import numpy as np
import itertools
from scipy.special import factorial as factorial
from copy import deepcopy
import time


# designed to make a real life network example of betweenness centrality and its importance
width = height = 10
nwLength = 10
percMultiple = 1.1
rcMean = 10
rcSD = 0

conds = {'connected' : False,
		'highestBCCausesBiggestChange' : False}
condsValues = list(conds.values())

iterCounter = 1
while all(condsValues) != True:
	print('Attempt', iterCounter)
	iterCounter += 1
	print('Creating new network to try')
	g = nwm.NanowireMesh(width = width,
				height = height,
				percMultiple = percMultiple,
				rcMean = rcMean,
				rcSD = rcSD)

	# resetting all the conditions to be false so nothing gets erroneously preserved from a previous run
	for key in conds.keys():
		conds[key] = False

	#calculating series resistance of network
	print('Calculating series resistance of network')
	netlistName = 'BC_example_network'
	g.to_netlist(netlistName = netlistName)
	command = ['xyce', netlistName]
	subprocess.check_call(command, stdout = subprocess.DEVNULL)
	g.update_with_xyce_output(netlistName + '.csv')


	# finding node with top betweenness centrality

	#finding bc values
	print('Calculating betweenness centralities')
	bc = nx.betweenness_centrality(g, weight = 'resistance')

	# removing electrodes from the dict
	del bc[g.topElectrode]
	del bc[g.bottomElectrode]

	# finding highest BC node
	bcVals = list(bc.values())
	bcKeys = list(bc.keys())
	sortedInds = np.argsort(bcVals)
	bcKeys = np.array(bcKeys)
	bcKeys = bcKeys[sortedInds]
	bcKeys = list(bcKeys)
	bcKeys.reverse()
	topBCNode = bcKeys[0]

	#iterate through nodes to find which one gives the greatest voltage change when removed
	# this list is in decreasing order of betweenness centrality
	# first set condition highestBCCausesBiggestChange to True. It will only remain true if 
	# the highest BC node does indeed cause the biggest resistnace change
	conds['highestBCCausesBiggestChange'] = True
	print('Calculating resistance changes when nodes are removed')
	for n, node in enumerate(bcKeys):
		h = deepcopy(g)
		h.remove_node(node)
		h.find_percolating_cluster()
		connected = h.isPercolating
		if not connected and node != topBCNode:
			print('Circuit disconnected when a single node was removed that was not the top BC node. Starting over.')
			print('Beginning attempt', iterCounter)
			time.sleep(1)
			break
		else:
			#calculating series resistance of network
			netlistName = 'BC_example_network_1_node_removed'
			h.to_netlist(netlistName = netlistName)
			command = ['xyce', netlistName]
			subprocess.check_call(command, stdout = subprocess.DEVNULL)
			h.update_with_xyce_output(netlistName + '.csv')

			# compare to reference voltage
			if node == topBCNode:
				referenceResistanceChange = h.sheetResistance - g.sheetResistance
			else:
				resistanceChange = h.sheetResistance - g.sheetResistance
				if resistanceChange > referenceResistanceChange:
					conds['highestBCCausesBiggestChange'] = False
					print('Found node which causes greater resistance change than node with highest betweenness centrality. Starting over')
					print('Beginning attempt', iterCounter)
					time.sleep(1)
					break

	shouldWrite = all(list(conds.values()))
	if shouldWrite:
		h.to_pickle(outFileName = 'BC_example_network')
