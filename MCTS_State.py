import NanowireMesh as nwm
import pdb
import networkx as nx
import numpy as np
from copy import copy, deepcopy
import sys
import pandas as pd

# utility functions
# the below are copied from bc_paper_node_removal_data_v4.py commit 39cebcd
# (this is when the data was generated) # docstrings are new though
xycePath = '/home/amwt/TPV/2DNanowires/installations/build_Xyce-6.11_parallel_openMPI/bin/Xyce'

def valid(g, msg):
	assert {g.topElectrode, g.bottomElectrode}.issubset(set(g.nodes)), 'Electrodes removed {}'.format(msg)

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
	try:
		g.find_percolating_cluster()
	except KeyError:
		pdb.set_trace()
	nonPercolating = set(g.nodes) - set(g.percolatingCluster)
	g.remove_nodes_from(nonPercolating)
	valid(g, 'pt1')


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
		valid(g, 'pt2')
		danglingEnds = {node for node in g.nodes if g.degree(node) == 1} 
		danglingEnds = {node for node in danglingEnds if all([g.degree(neighb) <= 2 for neighb in get_wire_segments(g, node)])}
		danglingEnds -= {g.topElectrode, g.bottomElectrode}
		nodesToRemove = {seg for node in danglingEnds for seg in get_wire_segments(g, node)}
#
#		nodesToRemove = {node for node in g.nodes if g.degree(node) in [0,1]}
#		nodesToRemove -= {g.topElectrode, g.bottomElectrode}

# the below functions are copied from bc_paper_node_removal_analysis_v3.py
# on commit 86294e8

transparencyParams = dict(wavelength = 550E-9,
				mr = 0.055 + 3.32j)


# these below parameters are all copied exactly from the runtimeOptions file in bc_paper_network_generator_001
# full path below
#/home/amwt/TPV/2DNanowires/Data/evolution_under_node_removal/bc_paper_network_generator_001/runtimeOptions.csv
params = {'width' : 100,  # in um
		'height' : 100, # in um
		'nwLength' : 10, # in um
		'nwDiam' : 0.15 # in um
}

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
	converts the effective percolation multiple to a density in units of \
	# wires per um^2
	"""
	assert height == width, "Height and width must be equal"
	nc_eff = 5.63726 / nwLength**2 + 1 / (nwLength * height) + 5.5 / (height**2)
	n_s = percMultiple * nc_eff
	return n_s


def fom(resistance, transparency):
	return 188.5 / (resistance * (transparency**(-1/2) - 1))

# beginning of truly new code
class State():
	def __init__(self, g, parent = None):
		self.g = g
		self.player = 1
		self.parent = parent
		# adding random suffix to netlist so that concurrently running xyce instances don't conflict
		if parent is not None:
			self.netlist_suffix = parent.netlist_suffix
		else:
			self.netlist_suffix = str(
				int(np.random.randint(low = 0, high = 1e10))
				).zfill(11)
		
	def getCurrentPlayer(self):
		return self.player
	
	def getPossibleActions(self):
		"""
		assert player == 1,	then return a list of frozensets of nodes like [ {n1, n2, n3}, {n4, n5}, ...] \
		where n1,n2,n3 are all the nodes of a single nanowire, and n4,4 are all the nodes of\
		another nanowire
		"""
		
		assert self.getCurrentPlayer() == 1
		g = self.g
		possibleNodes = set(g.nodes) - {g.topElectrode, g.bottomElectrode}
		actions = []
		while possibleNodes:
			startingNode = possibleNodes.pop()
			singleWireNodes = get_wire_segments(g, startingNode)
			possibleNodes -= set(singleWireNodes)
			actions.append(frozenset(singleWireNodes))
		return actions



	def takeAction(self, action):
		"""
		returns a new instance of state containing a copy of graph with \
			nodes in action removed
		sets parent or returned state to self
		"""
		h = nwm.copy(self.g)
		h.topElectrode = self.g.topElectrode
		h.bottomElectrode = self.g.bottomElectrode
		h.remove_nodes_from(action)
		# if the state is non terminal, clean the network
		# otherwise just return the new state without cleaning
		if nx.has_path(h, h.topElectrode, h.bottomElectrode):
			clean_network(h)
		else:
			print('terminal!!!')
		return State(
			h,
			parent = self
		)


	def isTerminal(self):
		"""
		returns False is g.topElectrode has a path to g.bottomElectrode
		returns True if electrodes are disconnected
		"""
		return not nx.has_path(self.g, self.g.topElectrode, self.g.bottomElectrode)

	def getReward(self):
		"""
		return TCE figure of merit of parent state
		note that FOMs are all divided by 200 to ensure that values fall on the interval [0,1]
		for which UCB for trees holds
		need to confirm that nothing will have FOM over 200 though
		"""
		g = copy(self.parent.g)
		if not nx.has_path(g, g.topElectrode, g.bottomElectrode):
			pdb.set_trace()
		g.solve_circuit_using_xyce(
			xycePath = xycePath,
			netlistName = 'netlist_mcts_{}'.format(self.netlist_suffix)
		)
		resistance = g.sheetResistance
		percMultiple = get_perc_multiple(g)
		density = effective_perc_multiple_to_density(
			percMultiple, 
			nwLength = params['nwLength'] * 1E-6,  # because nw Len given in um in params
			height = params['height'] * 1E-6, # because height given in um in params
			width = params['width'] * 1E-6 # because width given in um in params
		)
		transparency = nwm.transparency(radius = 0.15E-6 / 2,
				n_s = density,
				nwLength = params['nwLength'] * 1E-6,
				wavelength = transparencyParams['wavelength'],
				mr = transparencyParams['mr']
		)
		return fom(resistance, transparency)/200
