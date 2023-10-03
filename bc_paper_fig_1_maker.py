import NanowireMesh as nwm
import os
import networkx as nx
import numpy as np
import itertools
from scipy.special import factorial as factorial
import DrawFromDistribution as dfd

# choosing which bits of code to execute
shouldMakeNewNetwork = False
solveCircuit = True
calculateCentralities = True
loadElectricalProperties = True
findShortestPaths = True
makeCurrentPlot = True
makePowerPlot = True
makeBetweennessCentralityPlot = True
makePercolationCentralityPlot = True
randomizeJunctionResistances = True
makeShortestPathsPlots = True
prefix = '/Users/adamtrebach/Documents/Research/TPV/Nanowires/'

if not shouldMakeNewNetwork:
	pickleName = prefix + '200x200_percMult1.5_nwLen10_Rc10.p'
	g = nwm.NanowireMesh(inPickle = pickleName)
else:
	width = 200
	height = 200
	nwLength = 10
	nwLengthSD = 0
	nwDiam = .15,
	percMultiple = 1.5
	rcMean = 10
	rcSD = 0
	removeIsolates = True
	g = nwm.NanowireMesh( width = width,
				height = height,
				nwLength = nwLength,
				nwLengthSD = nwLengthSD,
				nwDiam = nwDiam,
				percMultiple = percMultiple,
				rcMean = 10,
				rcSD = 0,
				removeIsolates = removeIsolates)
	pickleName = prefix + '_'.join([str(int(width)) + 'x' + str(int(height)),
					'percMult' + str(percMultiple),
					'nwLen' + str(nwLength),
					'Rc' + str(rcMean)])
	g.to_pickle(outFileName = pickleName)



# randomizing junction resistances
if randomizeJunctionResistances:
	pickleName = pickleName.replace('.p', '_randomized_junction_resistances.p')

	mean1 = 10
	mean2 = 300
	sd1 = 0.5
	sd2 = 150
	print('randomizing junction resistances')
	def resistance_distributuion(x,y):
		pdf1 = 1/np.sqrt(2 * np.pi * sd1**2) * np.exp( - (x - mean1)**2 / (2 * sd1**2))
		pdf2 = 1/np.sqrt(2 * np.pi * sd2**2) * np.exp( - (x - mean2)**2 / (2 * sd2**2))
		return 1/2 * pdf1 + 1/2 * pdf2

	def rescaled_resistance_distribution(x,y):
		maxVal = resistance_distributuion(mean1, None)
		return 1 / maxVal * resistance_distributuion(x, y)

	contactResistors = [edge  for edge in g.edges if g.edges[edge]['resistanceType'] == 'cont']
	numDraws = len(contactResistors)
	xBounds = [0.000001, 10000]
	yBounds = [1,2]
	resistanceValues, trash = dfd.draw(rescaled_resistance_distribution,
					xBounds = xBounds,
					yBounds = yBounds,
					numDraws = numDraws)

	resistanceDict = {contactResistors[n] : resistanceValues[n] for n in range(len(contactResistors))}
	nx.set_edge_attributes(g, resistanceDict, 'resistance')
	print('finished randomizing junction resistances')
	g.to_pickle(pickleName)


# choosing whether or not to solve the circuit
netlistName = pickleName.replace('.p','') + '_netlist'
if solveCircuit:
	g.to_netlist(netlistName = netlistName)
	os.system('xyce ' + netlistName)

# reading from the solved circuit file
if loadElectricalProperties:
	g.update_with_xyce_output(netlistName + '.csv')
	print(g.nodes[g.topElectrode]['voltage'])
	g.to_pickle(pickleName)

# calculate betweenness centralities
if calculateCentralities:
	print('Calculating Betweenness Centrality')
	btwnCentrality = nx.betweenness_centrality(g, k = 300, weight = 'resistance')
	nx.set_node_attributes(g, btwnCentrality, 'betweenness_centrality')
	g.to_pickle(pickleName)
	
	print('Calculating Percolation Centrality')
	percCentrality = nx.percolation_centrality(g, attribute = 'voltage', weight = 'resistance')
	nx.set_node_attributes(g, percCentrality, 'percolation_centrality') 
	g.to_pickle(pickleName)

# calculating currents through edges from powers
edgeCurrents = {edge : np.sqrt(g.edges[edge]['power'] / g.edges[edge]['resistance']) for edge in g.edges}
nx.set_edge_attributes(g, edgeCurrents, 'current')
nodeCurrents = dict()
nodePowers = dict()
print('Averaging edge attributes to apply to nodes')
for node in g.nodes:
	edges = g.edges(node)
	if len(edges) > 0:
		thisNodeCurrents = [g.edges[edge]['current'] for edge in edges]
		avgCurrent = np.average(thisNodeCurrents)
		nodeCurrents[node] = avgCurrent

		thisNodePowers = [g.edges[edge]['power'] for edge in edges]
		avgPower = np.average(thisNodePowers)
		nodePowers[node] = avgPower
	else:
		nodeCurrents[node] = 0
nx.set_node_attributes(g, nodeCurrents, 'current')
nx.set_node_attributes(g, nodePowers, 'power')
g.to_pickle(pickleName)
# make image of currents (not informative it turns out)
if makeCurrentPlot:
	g.to_img(z = 'current',
		lineWidth = 4,
		xLabel = 'X [um]',
		yLabel = 'Y [um]',
		title = 'Color Coded by Current Carried', 
		outFile = pickleName.replace('.p', '_current_plot'),
		zLabel = 'Current [Amps]')

if makePowerPlot:
	g.to_img(z = 'power',
	lineWidth = 4,
	xLabel = 'X [um]',
	yLabel = 'Y [um]',
	title = 'Color Coded by Power Dissipated',
	zLabel = 'Power [Watts]',
	outFile = pickleName.replace('.p', '_power_plot'))

# find the k shortest simple paths
#def k_shortest_paths(G, weight, k):
#	return list(itertools.islice(nx.shortest_simple_paths(G, G.topElectrode, G.bottomElectrode, weight = 'resistance'), k))
#
#def nodes_in_k_shortest_paths(G, weight, k):
#	paths = k_shortest_paths(G, weight, k)
#	nodes = set()
#	for path in paths:
#		nodes = nodes.union(set(path))
#	return nodes

# instead of finding the k shortest paths, we're going to find shortest paths until the paths include 5% of all nodes 
if findShortestPaths:
	nodesInShortestPaths = set()
	shortestPaths = list()
	print('Beginning to find shortest paths')
	shortestPathsIterator = nx.shortest_simple_paths(g, g.topElectrode, g.bottomElectrode, weight = 'resistance')
	while len(nodesInShortestPaths) / len(g.nodes) < 0.02:
		thisShortestPath = next(shortestPathsIterator)
		shortestPaths.append(thisShortestPath)
		nodesInShortestPaths = nodesInShortestPaths | set(thisShortestPath)
		#for path in shortestPaths:
#			nodesInShortestPaths = nodesInShortestPaths.union(set(path))
#		
#		g.shortestPaths = shortestPaths
#		g.nodesInShortestPaths = nodesInShortestPaths
		print('Found a path. List of shortest paths includes', str(len(nodesInShortestPaths)/len(g.nodes) * 100), '% of all nodes')

	g.shortestPaths = shortestPaths
	g.nodesInShortestPaths = nodesInShortestPaths
	for node in g.nodes:
		g.nodes[node]['in_shortest_paths'] = 1 if node in nodesInShortestPaths else 0
	
	g.to_pickle(pickleName)

# assigning powers to nodes 
nodePowers = dict()
for node in g.nodes:
	edges = g.edges(node)
	if len(edges) > 0:
		thisNodePowers = [g.edges[edge]['power'] for edge in edges]
		totalPower = sum(thisNodePowers) / 2
		nodePowers[node] = totalPower
	else:
		nodePowers[node] = 0
nx.set_node_attributes(g, nodePowers, 'power')
g.to_pickle(pickleName)

# make shortest path plot
if makeShortestPathsPlots:
	# make plot of just shortest paths
		
	extraCode = ['set label 1 \"Shortest paths are calculated until at least 2% of nodes are included in the shortest paths\" at -30,-25',
			'set label 1 font \",120\"',
			'set label 1 front',
			'show label 1',
			'unset colorbox']
	
	g.to_img(z = 'in_shortest_paths',
		xLabel = 'X [um]',
		yLabel = 'Y [um]',
		zLabel = '',
		title = 'Shortest Paths',
		outFile = pickleName.replace('.p', 'shortest_paths_no_junctions'),
		extraCode = extraCode)

	# making the same plot as above but with the junctions shown
	# finding hottest 2% of points
	numPoints = len(g.nodesInShortestPaths) / len(g.nodes) * len(g.edges)
	numPoints = int(numPoints)


	# getting powers from graph
	powersDict = nx.get_edge_attributes(g, 'power')
	powers = list(powersDict.values())
	sortedInds = np.argsort(powers)
	edges =list(powersDict.keys())
	edges = np.array(edges)
	edgesSorted = edges[sortedInds,:]
	edgesToShow = edgesSorted[-1 - numPoints: -1,:]
	#converting to list of tuples
	edgesToShow = [(row[0], row[1]) for row in edgesToShow]
	
	# finding which points from the above lie on the shortest paths
	hottestEdgesInShortestPaths = [edge for edge in edgesToShow if edge[0] in g.nodesInShortestPaths and edge[1] in g.nodesInShortestPaths]

	# calculating the probability of this happening by chance
	# m = number of hot spots on shortest path
	# N = number of nodes in the entire network
	# n = number of nodes in shortest paths
	m = len(hottestEdgesInShortestPaths)
	N = len(g.nodes)
	n = len(g.nodesInShortestPaths)
	print('m, N, n', m, N, n)
#	numerator = np.sqrt(n * (3 * n /2 - m - 1)) * (3 * n / 2)**(3 * n/2) * (3 * N/2 - m - 1)**(3 * N / 2 - m - 1)
#	denominator = np.sqrt(N * (3 * N / 2 - m - 1)) * (3 * N / 2)**(3 * N / 2) * (3 * n / 2 - m - 1)**(3 * n / 2 - m - 1)
#	probability = numerator / denominator

	percentHottestEdgesInShortestPaths = round(len(hottestEdgesInShortestPaths) / len(edgesToShow) * 100,2)
	# making the plot
	extraCode = extraCode + ['set label 2 \"The same percentage of hottest junctions are included. ' + str(percentHottestEdgesInShortestPaths) + '% of these junctions lie on the shortest paths.\" at -30,-35',
			'set label 2 font \",120\"',
			'set label 2 front',
			'show label 2']
	
	g.to_img(z = 'in_shortest_paths',
		xLabel = 'X [um]',
		yLabel = 'Y [um]',
		zLabel = '',
		showJunctions = edgesToShow,
		title = 'Shortest Paths and Hot Spots',
		outFile = pickleName.replace('.p', 'shortest_paths_with_junctions'),
		extraCode = extraCode)


