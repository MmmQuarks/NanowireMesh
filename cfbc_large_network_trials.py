import sys
import NanowireMesh as nwm
import re
import pdb
import pandas as pd
import os
import networkx as nx
import approximate_current_flow_betweenness_centrality as acfbc

# program outline
#1. read from runtimeOptions.csv in resultsFolder as specified from sys.argv[1]
#2. generate network with properties set in runtimeOptions.csv
	#2a clean network so it has only the percolating cluster and no dangling ends
	# insert admittance values
#3. calculate approximate current flow betweenness centrality on nodes
#4. write cfbc values to graph and store as pickle in resultsFolder

# reading the results folder from the submit script

#1 reading from runtime options
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

# defining function for step 2a
def prep(graph):
	def _dangling_ends(graph):
		return [node for node in graph if graph.degree(node) == 1]

	# remove non percolating
	graph.remove_nodes_from(set(graph.nodes) - set(graph.percolatingCluster))

	# remove dangling ends
	while _dangling_ends(graph):
		graph.remove_nodes_from(
			_dangling_ends(graph)
			)

	admittanceDict = {edge : 1/g.edges[edge]['resistance'] for edge in g.edges}
	nx.set_edge_attributes(
		g,
		admittanceDict,
		'admittance'
		)

# 2 generating graph
g = nwm.NanowireMesh(
	makeEmpty = False,
	width = runtimeOptions['width'],
	height = runtimeOptions['height'],
	nwLength = runtimeOptions['nwLength'],
	nwDiam = runtimeOptions['nwDiam'],
	percMultiple = runtimeOptions['percMultiple']
	)
# 2a - prepping network
prep(g)

# 3 - calculating cfbc
k = acfbc.min_k(g,
	epsilon = 0.05,
	failureProbability = 0.01
	)
print(k)
cfbc = acfbc.approximate_current_flow_betweenness_centrality(
	g,
	normalized = True,
	weight = 'admittance',
	solver = 'lu',
	k = acfbc.min_k(
		g,
		epsilon = 0.05,
		failureProbability = 0.01
		)
	)

nx.set_node_attributes(g,
	cfbc,
	'cfbc_approx'
	)

g.to_pickle(
	outFileName = resultsFolder + '_'.join([
		str(int(runtimeOptions['width'])) + 'x' + str(int(runtimeOptions['height'])),
		'percMultiple' + str(runtimeOptions['percMultiple'])
		]
		)
	)
print('calculated and saved approx cfbc')

cfbcExact = nx.current_flow_betweenness_centrality(
	g,
	normalized = True,
	weight = 'admittance',
	solver = 'lu',
	)
nx.set_node_attributes(g,
	cfbcExact,
	'cfbc_exact'
	)

g.to_pickle(
	outFileName = resultsFolder + '_'.join([
		str(int(runtimeOptions['width'])) + 'x' + str(int(runtimeOptions['height'])),
		'percMultiple' + str(runtimeOptions['percMultiple'])
		]
		)
	)

print('calculated and saved exact cfbc')


print('program exited successfully')
