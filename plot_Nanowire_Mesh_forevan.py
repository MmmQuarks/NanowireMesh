import NanowireMesh as nwm
import numpy as np
from matplotlib import collections as mc
import pylab as pl
import matplotlib.pyplot as plt
import networkx as nx
import os
import pandas as pd


for n in range(467):
	lines = []
	colors = []
	folderName = '/home/amwt/TPV/2DNanowires/Frames3'
	# fileName = folderName + '/' + str(int(n)) + '_plotData.csv'
	# df = pd.read_csv(fileName)
	g = nwm.NanowireMesh(inPickle = folderName + '/' + str(int(n)) + '.p')

	endpoints = [g.nodes[node]['endpoints'] for node in g.nodes]
	# endpoints = [[(df.x1[i], df.y1[i]), (df.x2[i], df.y2[i])] for i in range(len(df))]

	lines = endpoints
	isPerc = np.isin(g.nodes, g.percolatingCluster)
	nodesList = list(g.nodes)
	percDict = {nodesList[i] : isPerc[i] for i in range(len(nodesList))}

	colors = [ ( (g.nodes[node]['temp'] - 298.15)/(498.15-298.15), 0, 0)  if percDict[node] else (0,0,1) for node in g.nodes]

	colors = [(0,0,0) if color[0] < 0 else color for color in colors]


	# colors = [( (df.temp[i] - 298.15)/(498.15 - 298.15) , 0 , 0) for i in range(len(df))]

	lineWidths = [color[0] *5 + .1 for color in colors]
	# lineWidths = [ color[0] * 5 + .5 for color in colors]


	# plot wires and color them black for percolating, blue for non percolating
	# for node in g.nodes:
	# 	#p1 = g.nodes[node]['endpoints'][0]
	# 	#p2 = g.nodes[node]['endpoints'][1]
	# 	lines.append(g.nodes[node]['endpoints'])

	# 	if node in g.percolatingCluster:
	# 		colors.append('black')
	# 	else:
	# 		colors.append('red')
	print('Making Plot')
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

	lc = mc.LineCollection(lines, alpha = .4, colors = colors, linewidths = lineWidths)
	ax.add_collection(lc)
	ax.autoscale()

	ax.set_xlim(0, 1000)
	ax.set_ylim(0, 1000)
	# plt.show()
	plotName = folderName + '/plots/' + str(int(n)) + '_plot.png'
	plt.savefig(plotName, dpi = 300)
	plt.close(fig)
