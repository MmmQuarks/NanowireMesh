import NanowireMesh as nwm
nx = nwm.nx
import pdb

def add_shunts(g):
	nodesToShunt = set(g.nodes) - set([g.bottomElectrode])

def prepare(g,
	removeDanglingEnds = True):
	# reomve everything not in the percolating cluster
	g.remove_nodes_from(
		set(g.nodes) - set(g.percolatingCluster)
		)
	# remove dangling ends
	if removeDanglingEnds:
		dangling_ends = lambda x : [node for node in x if x.degree(node) == 1]
		while dangling_ends(g):
			g.remove_nodes_from(
				dangling_ends(g)
				)

	# add shunts
	g.add_shunt_resistors()

	# add admittances
	nx.set_edge_attributes(
		g,
		{edge : 1/g.edges[edge]['resistance'] for edge in g.edges},
		name = 'admittance'
		)

	


g = nwm.NanowireMesh(makeEmpty = False,
	width = 20,
	height = 20)


#g.to_img(outFile = 'p1',
#	showJunctions = True)
prepare(g,
	removeDanglingEnds = False)
#g.to_img(outFile = 'p2',
#	showJunctions = True)


ebc = nx.edge_current_flow_betweenness_centrality(g,
	weight = 'admittance'
	)

pdb.set_trace()
