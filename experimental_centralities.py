import sys
import networkx as nx
import NanowireMesh as nwm
import ParallelCentrality as PC
from pathlib import Path
import PardisoSolver
from SolverUtils import prune_network
import pandas as pd


def main(**kwargs):
	# opening the network
	path = Path(kwargs['path'])
	g = nwm.NanowireMesh(makeEmpty = False, inPickle = path)


	# defining save function
	def save(c, label):
		new_path = Path( path.parent, path.stem + label + '.csv')
		pd.DataFrame(
			[dict(key = key, centrality = value) for key,value in c.items()]
			).to_csv(new_path)


	# calculating the vanilla bc for nodes and edges

	print('Calculating Node Betweenness Centrality')
	c = nx.betweenness_centrality(g, weight = 'resistance')
	save(c, 'node_betweenness_centrality')

	print('Calculating Edge Betweenness Centrality')
	c = nx.edge_betweenness_centrality(g, weight = 'resistance')
	save(c, 'edge_betweenness_centrality')


	# cleaning network so no isolates and no dangling ends (which all have undefined voltages)
	g, popped_nodes, popped_edges = prune_network(g)

	# calculating centralities
	print('Calculating Electrode Centrality')
	c = nwm.electrode_centrality(
		g,
		potential = 'voltage',
		weight = 'resistance'
		)
	save(c, 'electrode_centrality')

	print('Calculating Percolation Centrality')
	c = nx.percolation_centrality(
		g,
		attribute = 'voltage',
		weight = 'resistance'
		)
	save(c, 'percolation_centrality')

	print('Calculating Current Weighted Centrality')
	c = PC.current_weighted_centrality(g)
	save(c, 'current_weighted_centrality')

	print('Calculating Power Weighted Centrality')
	c = PC.power_weighted_centrality(g)
	save(c, 'power_weighted_centrality')
	
	print('Calculating Eigenvector Centrality')
	for edge in g.edges:
		g.edges[edge]['admittance'] = 1 / g.edges[edge]['resistance']
	c = nx.eigenvector_centrality(
		g,
		weight = 'admittance'
		)
	save(c, 'eigenvector_centrality')

	print('All Centralities Calculated Successfully')
	
if __name__ == '__main__':
	main(path = sys.argv[1])
	
