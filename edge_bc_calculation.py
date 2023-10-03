import sys
import NanowireMesh as nwm
import re
import pdb
import pandas as pd
import os
import networkx as nx

# program outline
#1. find original network files for a given taskID
#2. calculate edge betweenness centrality on original networks
#3. save new networks into destination folder

# reading the results folder from the submit script
resultsFolder = sys.argv[1]
if resultsFolder[-1] != '/':
	resultsFolder = resultsFolder + '/'
runtimeOptionsDf = pd.read_csv(resultsFolder + 'runtimeOptions.csv')
# columns 
runtimeOptions = dict(zip(runtimeOptionsDf.iloc[:,0], runtimeOptionsDf.iloc[:,1]))
for key,val in runtimeOptions.items(): # converting strings into other data types if possible
	if val in ['True', 'False']:
		runtimeOptions[key] = val == 'True' # converting string to bool
	else:
		try:
			runtimeOptions[key] = float(val)
		except ValueError:
			pass # this is not a number

taskID = sys.argv[2]
srcNetworkFolder = runtimeOptions['srcNetworkFolder']

def dangling_ends(graph):
	return [node for node in graph if graph.degree(node) == 1]

#1. starting by iterating through original set of networks (and making sure that we are only considering this task)
for origNetworkName in os.listdir(srcNetworkFolder):
	correctTask = 'task_' + taskID.zfill(2) in origNetworkName
	isPickle = origNetworkName[-2:] == '.p'
	if correctTask and isPickle:
		# open network before it has evolved
		origNetwork = nwm.NanowireMesh(inPickle = srcNetworkFolder + origNetworkName,
			makeEmpty = False)

		# this run we are calculating current flow betweenness centrality so we
		# have to remove disconnected nodes
		origNetwork.remove_nodes_from(set(origNetwork.nodes) - set(origNetwork.percolatingCluster))
		# also have to remove dangling ends
		while dangling_ends(origNetwork):
			origNetwork.remove_nodes_from(
				dangling_ends(origNetwork)
				)
		

#2 calculating edge betweenness centrality
		ebc = nx.edge_current_flow_betweenness_centrality(origNetwork,
			normalized = True,
			solver = 'full',
			weight = 'admittance')

		nx.set_edge_attributes(origNetwork, # writing ebc values to graph
			ebc,
			name = 'edge_current_flow_betweenness_centrality')

		taskTrialString = origNetworkName[:17]

		origNetwork.to_pickle(outFileName = resultsFolder + origNetworkName)
