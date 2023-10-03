import os
import networkx as nx
from matplotlib import pyplot as plt
import numpy as np
import itertools 
# making the polygonal graph
g = nx.Graph()
numNodes = 11

# making exterior nodes and all edges
for node in range(1, numNodes + 1):
	edge = (node, node + 1) if node < numNodes else (node, 1)
	g.add_edge(*edge)

	# making the coordinates from the unit circle
	angle = (node - 1) * 2 * np.pi / numNodes
	radius = 300
	g.nodes[node]['x'] = radius * np.cos(angle)
	g.nodes[node]['y'] = radius * np.sin(angle)

	# adding edge to central node
	g.add_edge(node, 0)

# making properties of central node
g.nodes[0]['x'] = g.nodes[0]['y'] = 0
pos = {node : (g.nodes[node]['x'], g.nodes[node]['y']) for node in g.nodes}


# making matplotlib axes with settings
ax = plt.axes()
ax.set_xlim(xmin = -1.5 * radius, xmax = 1.5 * radius)
ax.set_ylim(ymin = -1.5 * radius, ymax = 1.5 * radius)
ax.set_aspect('equal')

# first we just plot the darn thing
nx.draw_networkx(g, with_labels = False, ax = ax, pos = pos)
plt.savefig('2A.png', dpi = 300)
os.system('open 2A.png')
plt.clf()

# making matplotlib axes with settings
ax = plt.axes()
ax.set_xlim(xmin = -1.5 * radius, xmax = 1.5 * radius)
ax.set_ylim(ymin = -1.5 * radius, ymax = 1.5 * radius)
ax.set_aspect('equal')


# calculating betweenness centrality
bc = nx.betweenness_centrality(g)
bc = {node : round(val,2) for node, val in bc.items()}

# making node sizes proportional to betweenness centrality
nodeSizeList = [ 3000* bc[node] for node in g.nodes]

# calculating average internode distance
distances = [nx.shortest_path_length(g, source = edge[0], target = edge[1]) for edge in itertools.combinations(list(g.nodes), 2)]
avgDistance = np.average(distances)

# making text box with path length data
ax.text(-radius, -1.3 * radius, 'Average Distance Between Nodes = ' + str(round(avgDistance, 2)), verticalalignment = 'top')

# then we plot with sizes proportional to betweenness centrality
nx.draw_networkx(g, with_labels = False, ax = ax, pos = pos, node_size = nodeSizeList, nodelist = list(g.nodes))
plt.savefig('2B.png', dpi = 300)
os.system('open 2B.png')
plt.clf()

# making matplotlib axes with settings
ax = plt.axes()
ax.set_xlim(xmin = -1.5 * radius, xmax = 1.5 * radius)
ax.set_ylim(ymin = -1.5 * radius, ymax = 1.5 * radius)
ax.set_aspect('equal')


# remove an exterior node and see how the average internode distance changes
def filter_node(node):
	return node != 1
h = nx.subgraph_view(g, filter_node = filter_node)

# calculating betweenness centrality
bc = nx.betweenness_centrality(h)
bc = {node : round(val,2) for node, val in bc.items()}

# making node sizes proportional to betweenness centrality
nodeSizeList = [ 3000* bc[node] for node in h.nodes]

# calculating average internode distance
distances = [nx.shortest_path_length(h, source = edge[0], target = edge[1]) for edge in itertools.combinations(list(h.nodes), 2)]
avgDistance = np.average(distances)

# making text box with path length data
ax.text(-radius, -1.3 * radius, 'Average Distance Between Nodes = ' + str(round(avgDistance, 2)), verticalalignment = 'top')

# then we plot with sizes proportional to betweenness centrality
nx.draw_networkx(h, with_labels = False, ax = ax, pos = pos )
plt.savefig('2C.png', dpi = 300)
os.system('open 2C.png')
plt.clf()

# making matplotlib axes with settings
ax = plt.axes()
ax.set_xlim(xmin = -1.5 * radius, xmax = 1.5 * radius)
ax.set_ylim(ymin = -1.5 * radius, ymax = 1.5 * radius)
ax.set_aspect('equal')


# remove the interior node and see how the average internode distance changes
def filter_node(node):
	return node != 0
h = nx.subgraph_view(g, filter_node = filter_node)

# calculating betweenness centrality
bc = nx.betweenness_centrality(h)
bc = {node : round(val,2) for node, val in bc.items()}

# making node sizes proportional to betweenness centrality
nodeSizeList = [ 3000* bc[node] for node in h.nodes]

# calculating average internode distance
distances = [nx.shortest_path_length(h, source = edge[0], target = edge[1]) for edge in itertools.combinations(list(h.nodes), 2)]
avgDistance = np.average(distances)

# making text box with path length data
ax.text(-radius, -1.3 * radius, 'Average Distance Between Nodes = ' + str(round(avgDistance, 2)), verticalalignment = 'top')

# then we plot with sizes proportional to betweenness centrality
nx.draw_networkx(h, with_labels = False, ax = ax, pos = pos)
plt.savefig('2D.png', dpi = 300)
os.system('open 2D.png')
plt.clf()
