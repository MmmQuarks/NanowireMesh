import NanowireMesh as nwm
import numpy as np
from matplotlib import collections as mc
import pylab as pl
import matplotlib.pyplot as plt
import networkx as nx
import os

g = nwm.NanowireMesh(height = 250, width = 250, nwLength = 7, nwLengthSD = 3, percMultiple = 1.5, removeNonPercolating=False, addInternalRes= True)

lines = []
colors = []

# plot wires and color them black for percolating, blue for non percolating
for node in g.nodes:
	#p1 = g.nodes[node]['endpoints'][0]
	#p2 = g.nodes[node]['endpoints'][1]
	lines.append(g.nodes[node]['endpoints'])

	if node in g.percolatingCluster:
		colors.append('black')
	else:
		colors.append('red')

fig, ax = pl.subplots()
ax.set_aspect(aspect = 1)
# plot junctions
# plt.scatter( nx.get_edge_attributes(g, 'x').values(), nx.get_edge_attributes(g, 'y').values(), c = 'blue', s = 1)

#plot voltages (but first must solve with xyce)
#g.toNetlist('netlist')
#os.system('xyce netlist')
#g.updateNetworkWithXyceOutput('netlist.csv')
# plot voltages
# plt.scatter( 
# 	list(nx.get_node_attributes(g,'x').values()), 
# 	list(nx.get_node_attributes(g,'y').values()), 
# 	s = 1,
# 	c = list(nx.get_node_attributes(g,'voltage').values())
# )

#plot voltages ONLY of nodes at voltage 1
# x = [g.nodes[node]['x'] for node in g.nodes if g.nodes[node]['voltage']==1]
# y = [g.nodes[node]['y'] for node in g.nodes if g.nodes[node]['voltage']==1]
# plt.scatter(x,y, s = 10)

lc = mc.LineCollection(lines, alpha = 0.3, colors = colors, linewidths = .5)
ax.add_collection(lc)
ax.autoscale()

ax.set_xlim(0, g.width)
ax.set_ylim(0, g.height)
plt.show()