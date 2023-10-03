
def prune_network(g, pop = True):
	# removing non_percolating wires
	g.find_percolating_cluster(verbose = False)
	non_percolating = set(g.nodes) - set(g.percolatingCluster)
	popped_nodes, popped_edges = g.pop_nodes_with_edges(non_percolating)

	# removing nodes with degree 1
	degree_one_nodes = [node for node in g if g.degree(node) == 1]
	while degree_one_nodes:
		pn, pe = g.pop_nodes_with_edges(degree_one_nodes)
		popped_nodes += pn
		popped_edges == pe
		degree_one_nodes = [node for node in g if g.degree(node) == 1]

	if pop:
		return g, popped_nodes, popped_edges
	else:
		return g
		

