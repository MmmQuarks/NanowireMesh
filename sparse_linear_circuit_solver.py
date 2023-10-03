import networkx as nx
import pdb
from scipy import sparse, linalg
import numpy as np
from copy import deepcopy
import numba


def _setup_admittance_matrices(g, currentInjections = None, isSparse = False):

	# make current vector assuming current injections is a dict mapping each node to the 
	# current injected there (positive sign) or withdrawn from there (negative sign)
	# if currentInjections = None we put 1 amp in the top electrode and withdraw
	# 1 amp from the bottom electrode
	if currentInjections == None:
		currentInjections = {g.topElectrode : 1,
					g.bottomElectrode : -1}


	# confirm that there exists a pathway from current injection site to current withdrawal site
	currentNodes = list(currentInjections.keys())
	assert nx.has_path(g, currentNodes[0], currentNodes[1]), 'There is no percolating path from current injection to current withdrawal'


	# confirm that both current injection/withdrawal nodes are in the percolating cluster
	g.find_percolating_cluster()
	assert currentNodes[0] in g.percolatingCluster, 'Current source/sink not in percolating cluster'
	assert currentNodes[1] in g.percolatingCluster, 'Current source/sink not in percolating cluster'
	
	# temporarily removing nodes not in the percolating cluster
	nonPercolatingNodes = set(g.nodes) - set(g.percolatingCluster)
	poppedNodes, poppedEdges = g.pop_nodes_with_edges(nonPercolatingNodes)

	
	nodeList = sorted(g.nodes)

	# off diagonal elements should be negative reciprocals of resistances
	# note the offDiags.power(-1) maps 0 -> 0 so don't need to deal with infinities
	offDiagMatrix = nx.adjacency_matrix(g, nodelist = nodeList, weight = 'resistance')
	offDiagMatrix = -1 * offDiagMatrix.power(-1)

	# diagonal elements in row n should be y_1n + y_2n + ... y_Nn
	# this is the negative sum of all elements in column (or row) n of offDiagMatrix
	diagElements = [-offDiagMatrix[:,n].sum() for n in range(offDiagMatrix.shape[0])]
	diagMatrix = sparse.diags(diagElements, format = 'csr')
	admittanceMatrix = offDiagMatrix + diagMatrix


	# using numpy to confirm that the choice of current injections is sensible	
	currentsArray = np.array(list(currentInjections.values()))
	# current injected must equal current withdrawn
	assert  currentsArray.sum() == 0, 'Current Injections violate conservation of current'

	# must be exactly two current injections
	assert currentsArray.shape[0] == 2, 'Must be exactly two current injection sites.'

	# current injections must be nonzero
	assert np.alltrue(currentsArray != 0), 'Current injections must all be nonzero'

	# make current array in the order of nodeList
	I = np.zeros(len(g.nodes))
	for n, node in enumerate(nodeList):
		try:
			I[n] = currentInjections[node]
		except KeyError:
			pass	# if the node isn't listed in the current injections the value should remain 0

	# scale rows of admittanceMatrix and I
	for row in range(admittanceMatrix.shape[0]):
		if nodeList[row] not in currentInjections.keys():
			scale = 1/sparse.linalg.norm(admittanceMatrix[row,:])
			admittanceMatrix[row,:] *= scale
			I[row] *= scale


	return admittanceMatrix if isSparse else admittanceMatrix.todense(), I, nodeList, poppedNodes, poppedEdges

@numba.jit(nopython = True)
def lstsq_numba(A, b):
	return np.linalg.lstsq(A, b, rcond = -1)

@numba.jit(nopython = True)
def solve_numba(A, b):
	return np.linalg.solve(A, b)
	
def solve_with_admittances(g, 
			currentInjections = None,
			nodeList = None, 
			showInfo = False,
			solver = 'spsolve',
			voltageName = 'voltage',
			currentName = 'current'):

	
	try:
		assert solver in ['solve', 'lstsq_numba', 'lstsq', 'solve_numba', 'solve_scipy', 'spsolve']
	except AssertionError as e:
		e.args += ('Unknown solver \"' + solver + '\" requested',)
		raise

	# saving some time if we are working with sparse matrices by not converting from sparse and then back to sparse
	# inside spsolve
	if solver in ['spsolve']:
		isSparse = True
	else:
		isSparse = False
	admittanceMatrix, I, nodeList, poppedNodes, poppedEdges = _setup_admittance_matrices(g,
												currentInjections = currentInjections,
												isSparse = isSparse)


	if solver == 'solve':
		x = np.linalg.solve(admittanceMatrix, I)
		info = 'Solved using numpy.linalg.solve. No extra info given.'
	elif solver == 'solve_numba':
		x = solve_numba(admittanceMatrix, I)
		info = 'Solved using numpy.linalg.solve. No extra info given.'
	elif solver == 'lstsq_numba':
		info = lstsq_numba(admittanceMatrix, I)
		x = info[0]
	elif solver == 'lstsq':
		info = np.linalg.lstsq(admittanceMatrix ,I, rcond = -1) 
		x = info[0]
	elif solver == 'solve_scipy':
		x = linalg.solve(admittanceMatrix, I, assume_a = 'sym')
		info = 'Solved using scipy.linalg.solve. No extra info given'
	elif solver == 'spsolve':
		x = sparse.linalg.spsolve(admittanceMatrix, I)
		info = 'Solved using scipy.sparse.linalg.spsolve. No extra info given'

	
	# sometimes the numerics give us a solution where a handful of nodes have lower potential than ground.
	# we fix that manually here
	withdrawalIndex = np.argmin(I)
	groundVoltage = x[withdrawalIndex]
	x[x < groundVoltage] = groundVoltage
	minVoltageIndex = np.argmin(x)
	try:
		assert np.alltrue(x >= groundVoltage), 'Current withdrawal site is not at lowest voltage'
	except AssertionError:
		pdb.set_trace()

	# resetting the zero of the voltages so that withdrawal site is at zero
	x = x - x.min()

	voltagesDict = {nodeList[n] : x[n] for n in range(len(nodeList))}
	nx.set_node_attributes(g, voltagesDict, voltageName)

	currentsDict = {edge : abs(voltagesDict[edge[0]] - voltagesDict[edge[1]]) / g.edges[edge]['resistance'] for edge in g.edges}
	nx.set_edge_attributes(g, currentsDict, currentName)

	# setting the correct attributes for the non percolating nodes and edges
	# that were popped earlier
	# note that the elements of these lists are [(node, {keys : vals}), ...]
	# or [(node1, node1, {keys, vals}), ...]
	for nodeTuple in poppedNodes:
		nodeTuple[-1][voltageName] = x.min()
	for edgeTuple in poppedEdges:
		edgeTuple[-1][currentName] = 0

	# adding these non percolating nodes back in
	g.add_nodes_from(poppedNodes)
	g.add_edges_from(poppedEdges)
	
	# recording sheet resistance assuming it is equal to series resistance
	

	if showInfo:
		return info

	

def main():
	import NanowireMesh as nwm
	dim = 75
	g = nwm.NanowireMesh(width = dim, height = dim, percMultiple = 1.5)
	# making example from textbook
#	g.remove_nodes_from(list(g.nodes))
#	ebunch = [(4,1, {'id' : 1, 'resistance' : 1, 'resistanceType' : 'cont'}),
#			(1,0, {'id' : 2, 'resistance' : 1, 'resistanceType' : 'cont'}),
#			(1,2, {'id' : 3, 'resistance' : 1, 'resistanceType' : 'cont'}),
#			(2,0, {'id' : 4, 'resistance' : 1, 'resistanceType' : 'cont'}),
#			(2,3, {'id' : 5, 'resistance' : 1, 'resistanceType' : 'cont'}),
#			(3,0, {'id' : 6, 'resistance' : 1, 'resistanceType' : 'cont'})]
#	g.add_edges_from(ebunch)
#	g.topElectrode = 4
#	g.bottomElectrode = 0

	solve(g)
if __name__ == "__main__":
	main()
