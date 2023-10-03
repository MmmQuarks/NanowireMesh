import networkx as nx
import NanowireMesh as nwm
from tqdm import tqdm
import pdb
from copy import deepcopy
import Timer
import itertools

def path_length(g, path, weight):
	# make list of edges
	edges = list( zip( path, path[1:]))
	lengths = [g.edges[e][weight] for e in edges]
	return sum(lengths)

def _parallel_all_pairs_shortest_path(g, weight = None):
#	for node in g:
#		try:
#			int(node)
#		except ValueError:
#			msg = ' '.join(['Node', str(node), 'is not int'])
#			raise Exception(msg)

	path = nx.single_source_dijkstra_path
	sortedNodes = sorted(list(g.nodes))
	b = db.from_sequence(sortedNodes, npartitions = 4)
	results = b.map(lambda n : path(g, n, cutoff = None, weight = weight)).compute()
	resultsDict = {sortedNodes[int(n)] : results[int(n)] for n in sortedNodes}
	return resultsDict


def betweenness_centrality(g, 
			weight = None, 
			parallel = False, 
			normalized = True):

	# note that this algorithm returns only one shortest path for each pair of nodes.
	# if there are multiple shortest paths, algorithm will not yield correct answers
	shortestPaths = nx.shortest_path(g, weight = weight)

	# make centrality container
	c = dict.fromkeys(g.nodes, 0)


	# in an undirected graph, the path with endpoints A&B only counts once
	# aka A->B is the same path as B->A. We will only count one of these in the below
	# and then if the graph is directed we will correct this with the normalization
	# at the end
	while shortestPaths:
		source = list(shortestPaths.keys())[0]
		targetsPathsDict = shortestPaths.pop(source)
		for target in targetsPathsDict.keys():
			path = targetsPathsDict[target]
			for node in path[1:-1]:
				c[node] += 1

			# removing the reversed version of this 
			# pair of endpoints from shortestPaths
			try:
				del shortestPaths[target][source]
			except KeyError:
				# handling the issues in case this pair has 
				# already been removed
				pass

	if normalized:
		N = len(g)
		# scaling the normalization correctly for directed and undirected
		numerator = 1 if g.is_directed() else 2
		for node in c.keys():
			c[node] *= numerator  / ((N-1) * (N-2))
	return c

			

# deprecated versions of functions. probably should delete but testing first.

#def current_weighted_centrality(g, endpoints = False, parallel = False):
#	weight = 'resistance'
#	def pathWeight(G, path, node):
#		source = path[0]
#		target = path[-1]
#		deltaV = G.nodes[source]['voltage'] - G.nodes[target]['voltage']
#		R = path_length(g, path = path, weight = weight)
#		return deltaV / R
#
#	def pathCondition(G, path):
#		source = path[0]
#		target = path[-1]
#		return g.nodes[source]['voltage'] > g.nodes[target]['voltage']
#
#	systemVoltage = g.nodes[g.topElectrode]['voltage'] - g.nodes[g.bottomElectrode]['voltage']
#	systemCurrent = systemVoltage / g.sheetResistance
#	normalizationCoefficient = 1/systemCurrent
#
#	return centrality(g,
#			weight = weight,
#			endpoints = endpoints,
#			parallel = parallel,
#			pathWeight = pathWeight,
#			pathCondition = pathCondition)
#
#def power_weighted_centrality(g, endpoints = False, parallel = False):
#	# same as current_weighted_centrality but with deltaV**2 at the end
#	weight = 'resistance'
#	def pathWeight(G, path, node):
#		source = path[0]
#		target = path[-1]
#		deltaV = G.nodes[source]['voltage'] - G.nodes[target]['voltage']
#		R = path_length(g, path = path, weight = weight)
#		return deltaV**2 / R
#
#	def pathCondition(G, path):
#		source = path[0]
#		target = path[-1]
#		return g.nodes[source]['voltage'] > g.nodes[target]['voltage']
#
#	systemVoltage = g.nodes[g.topElectrode]['voltage'] - g.nodes[g.bottomElectrode]['voltage']
#	systemPower = systemVoltage**2 / g.sheetResistance
#	normalizationCoefficient = 1/systemPower
#
#	return centrality(g,
#			weight = weight,
#			endpoints = endpoints,
#			parallel = parallel,
#			pathWeight = pathWeight,
#			pathCondition = pathCondition)

def percolation_centrality(g, 
			weight = 'resistance',
			attribute = 'voltage',
			targets = False):
	# targets determines whether we are using source and target percolation
	# or just source (aka is the weight V(s) or V(s) - V(t))
	
	if targets:
		return _percolation_centrality_with_targets(g, 
							weight = 'resistance',
							attribute = 'voltage')
	else:
		return _percolation_centrality_without_targets(g, 
							weight = 'resistance',
							attribute = 'voltage')


def _percolation_centrality_without_targets(g,
					weight = 'resistance',
					attribute = 'voltage'):

	shortestPaths = nx.shortest_path(g, weight = weight)

	# make centrality container
	c = dict.fromkeys(g.nodes, 0)

	# in an undirected graph, the path with endpoints A&B only counts once
	# aka A->B is the same path as B->A. We will only count one of these in the below
	# and then if the graph is directed we will correct this with the normalization
	# at the end
	attributeSum = sum([g.nodes[node][attribute] for node in g])
	for source in shortestPaths.keys():
		targetsPathsDict = shortestPaths[source]
		for target in targetsPathsDict.keys():
			path = targetsPathsDict[target]
			sourceVal = g.nodes[source][attribute]
			for node in path[1:-1]:
				nodeVal = g.nodes[node][attribute]
				c[node] += sourceVal / (attributeSum - nodeVal)


	N = len(g)
	# scaling the normalization correctly for directed and undirected
	for node in c.keys():
		try:
			assert c[node] is not None
		except AssertionError:
			print('Node with None percolation centrality.')
			pdb.set_trace()
		c[node] *= 1 / (N-2)
	return c



def _percolation_centrality_with_targets(g, 
					weight = 'resistance',
					attribute = 'voltage'):

	shortestPaths = nx.shortest_path(g, weight = weight)

	# make centrality container
	c = dict.fromkeys(g.nodes, 0)

	# PC(v) = sum_{s != v != t} sigma_st (v) / sigma_st * w_vst
	# w_vst = ramp(x_s - x_t) / sum_{s != v != t} ramp(x_s - x_t)
	# this sum in the weight is summed over all s and t with s!= t
	# and neither s nor t equal to v
	# note that we can include the terms with s = t and not change
	# the sum because these terms contribute zero then we can rewrite the 
	# denom of w_vst as 
	# denom = sum_{s != v, t != v} ramp(x_s - x_t)
	# denom = sum_{st} R(x_s - x_t) - sum_{t} R(x_v - x_t) - sum_{s} R(x_s - x_v)
	ramp = lambda x : x if x >= 0 else 0
	attributes = nx.get_node_attributes(g, attribute)
	pairSum = 0
	for s in g:
		for t in g:
			pairSum += ramp(attributes[s] - attributes[t])

	def sourceSum(node, attributes = attributes):
		return sum([ramp( attributes[source] - attributes[node]) for source in attributes.keys()])
	
	def targetSum(node, attributes = attributes):
		return sum([ramp( attributes[node] - attributes[target]) for target in attributes.keys()])

	for source in shortestPaths.keys():
		targetsPathsDict = shortestPaths[source]
		for target in targetsPathsDict.keys():
			path = targetsPathsDict[target]
			for node in path[1:-1]:
				numerator = ramp(g.nodes[source][attribute] - g.nodes[target][attribute])
				denominator = pairSum - sourceSum(node) - targetSum(node)
				c[node] += numerator / denominator


	return c

def current_weighted_centrality(g):
	shortestPaths = nx.shortest_path(g, weight = 'resistance')

	# make centrality container
	c = dict.fromkeys(g.nodes, 0)

	ramp = lambda x : x if x >= 0 else 0
	voltages = nx.get_node_attributes(g, 'voltage')

	for source in shortestPaths.keys():
		targetsPathsDict = shortestPaths[source]
		for target in targetsPathsDict.keys():
			# making sure we skip path from node to itself
			if target == source:
				continue
			path = targetsPathsDict[target]
			pathResistance = path_length(g, path = path, weight = 'resistance')
			pathVoltage = ramp(voltages[source] - voltages[target])
			try:
				pathWeight = pathVoltage / pathResistance
			except ZeroDivisionError:
				pdb.set_trace()
			
			for node in path[1:-1]:
				c[node] += pathWeight

	# scaling the normalization correctly for directed and undirected
	totalVoltage = voltages[g.topElectrode] - voltages[g.bottomElectrode]
	totalCurrent = totalVoltage / g.sheetResistance
	N = len(g)
	for node in c.keys():
		c[node] *= 1 / (N-2) / totalCurrent
	return c

def power_weighted_centrality(g):
	shortestPaths = nx.shortest_path(g, weight = 'resistance')

	# make centrality container
	c = dict.fromkeys(g.nodes, 0)

	ramp = lambda x : x if x >= 0 else 0
	voltages = nx.get_node_attributes(g, 'voltage')

	for source in shortestPaths.keys():
		targetsPathsDict = shortestPaths[source]
		for target in targetsPathsDict.keys():
			# making sure we skip path from node to itself
			if target == source:
				continue
			path = targetsPathsDict[target]
			pathResistance = path_length(g, path = path, weight = 'resistance')
			pathVoltage = ramp(voltages[source] - voltages[target])
			try:
				pathWeight = pathVoltage**2 / pathResistance
			except ZeroDivisionError:
				pdb.set_trace()
			
			for node in path[1:-1]:
				c[node] += pathWeight

	# scaling the normalization correctly for directed and undirected
	totalVoltage = voltages[g.topElectrode] - voltages[g.bottomElectrode]
	totalPower = totalVoltage**2 / g.sheetResistance
	N = len(g)
	for node in c.keys():
		c[node] *= 1 / (N-2) / totalPower
	return c



def main():
	g = nwm.NanowireMesh(width = 30, height = 30)
	results = dict()
	print('Testing betweenness centrality')
	for normalized in [True, False]:
		nxc = nx.betweenness_centrality(g,
						weight = 'resistance',
						normalized = normalized,
						endpoints = False)
		c = betweenness_centrality(g,
						weight = 'resistance',
						normalized = normalized)

		name = 'btwn centrality' + (' normalized' if normalized else '')
		results[name] = nxc == c

	g.solve_circuit_using_xyce()
	print('Testing percolation centrality without targets')
	nxc = nx.percolation_centrality(g,
					weight = 'resistance',
					attribute = 'voltage')
	c = percolation_centrality(g, 
					weight = 'resistance',
					attribute = 'voltage')
	# rounding nxc and c to 8 decimal places
	for key in nxc.keys():
		nxc[key] = round(nxc[key], 8)
	for key in c.keys():
		c[key] = round(c[key], 8)
	results['percolation centrality'] = nxc == c

	print('testing current weighted centrality')
	c = current_weighted_centrality(g) 
	print('----')	
	for key, val in results.items():
		print(key, ':',  val)





if __name__ == '__main__':
	main()
