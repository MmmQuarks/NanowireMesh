import NanowireMesh as nwm
import networkx as nx
from sortedcontainers import SortedDict


def get_wire_segments(g, startingNode):
	"""
	returns all nodes in the same nanowire as startingNode for the graph g
	"""
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
	g.find_percolating_cluster()
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

def get_perc_multiple(g):
	"""
	returns the effective percolation multiple of the graph g
	"""
	assert g.width == g.height, "Width and height are not the same"
	lTot = sum([g.nodes[node]['length'] for node in g if node not in [g.topElectrode, g.bottomElectrode]])
	l = g.nwLength
	Ls = g.width
	return lTot * l / (5.63726 * Ls**2 + l * Ls + 5.5 * l**2)


def effective_perc_multiple_to_density(percMultiple, nwLength, height, width):
	"""
	args
		percMultiple - size-normalized percolation multiple of the system
		nwLength - length of average nanowire IN METERS
		height - height of sample IN METERS
		width - width of sample IN METERS
	returns
		number density in wires per METER^2
	"""
	assert height == width, "Height and width must be equal"
	nc_eff = 5.63726 / nwLength**2 + 1 / (nwLength * height) + 5.5 / (height**2)
	n_s = percMultiple * nc_eff
	return n_s


def fom(resistance, transparency):
	return 188.5 / (resistance * (transparency**(-1/2) - 1))


def sort_nodes_by_attribute(G, attribute, includeElectrodes = False):
	attrDict = nx.get_node_attributes(G, attribute)
	if not includeElectrodes:
		del attrDict[G.topElectrode]
		del attrDict[G.bottomElectrode]
	# returns attributes sorted by increasing key
	return SortedDict({val : key for key, val in attrDict.items()})

