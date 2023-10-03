import NanowireMesh as nwm
import numpy as np
import networkx as nx
import os
import seaborn as sns
from matplotlib import pyplot as plt

width = height =1000 
f = 2/3 # the amount of vertical space to be filled by punctures
n = 5 # cutouts will have size of roughly nl x nl where l is wire length
nwLength = 10
nSquares = round(f * width / (n * nwLength) + 1)
radius = width / 6
# below we make a function which zeroes out the pdf along a slanted line in the middle of the sample. 
# the zeroed out region extends nwLength to the left and right of the central line for a total width of 
# 2 * nwLength
pdfList = [ lambda x,y : 1 - np.heaviside(n * nwLength / 2 - np.abs(x - width/3 + 1/3*(y-height)),0) * np.heaviside(-np.sin(2 * np.pi * nSquares * y / height),0),
	lambda x,y : 0.5*(1 +  np.heaviside(-np.sin(2 * np.pi * x/(2/5 * width)),0)),
	lambda x,y : 1 - 0.3 * np.heaviside(-np.sin(2 * np.pi * x/(2/5 * width)),0),
	lambda x,y : 0.5*(1 + np.heaviside(-np.sin(2 * np.pi * y/(2/5 * height)),0)),
	lambda x,y : 1 - 0.3 * np.heaviside(-np.sin(2 * np.pi * y/(2/5 * height)),0),
	lambda x,y : 0.5* (1 + np.heaviside(radius**2 - (x - width/2)**2 - (y - height/2)**2,0)),
	lambda x,y : 1 - 0.3 *  np.heaviside(radius**2 - (x - width/2)**2 - (y - height/2)**2,0)
	]

pdfNames = ['zipper',
		'horiz_stripes_1',
		'horiz_stripes_2',
		'vert_stripes_1',
		'vert_stripes_2',
		'circle_extra',
		'circle_hole'
		]
for n in range(len(pdfList)):
	g = nwm.NanowireMesh(width = width, height = height, pdf = pdfList[n], percMultiple = 2)
	# removing isolates
	isolates = list(nx.isolates(g))
	g.remove_nodes_from(isolates)
	
	# making betweenness centrality calculation and adding to graph
	bc = nx.betweenness_centrality(g, k = int(5 * np.log(len(g.nodes))), normalized = True, weight = 'resistance')
	attrs = {node : {'bc' : bc[node]} for node in g.nodes}
	nx.set_node_attributes(g, attrs)

	# making electrical calculations
	g.to_netlist(netlistName = 'netlist')
	os.system('xyce netlist')
	g.update_with_xyce_output('netlist.csv')
	# saving network so I can refer to it later
	g.to_pickle('bc_patt_' + pdfNames[n] + '.p')

	numBins = 50
	xBins = np.linspace(0, g.width, numBins + 1)
	yBins = np.linspace(0, g.height, numBins + 1)

	# making numpy arrays of node data
	# note that we are using np.digitize to get the bin into which each x and y coordinate will fall
	data = [(node, np.digitize(g.nodes[node]['x'], xBins), np.digitize(g.nodes[node]['y'], yBins), g.nodes[node]['bc']) for node in g.nodes]
	dtype = [('node', float), ('xInd', float), ('yInd', float), ('bc', float)]
	data = np.array(data, dtype = dtype)

	# sorting the data into its bins
	# note that it is in the same order as data 
#	inds = [(np.digitize(data['x'][n], xBins), np.digitize(data['y'][n], yBins)) for n in range(len(data))]
#	inds = np.array(inds, dtype = [('xInd', int), ('yInd', int)])
#	data = [(data['node'][n], data['x'][n], data['y'][n], data['bc'][n], inds['xInd'][n], inds['yInd'][n]) for n in range(len(data))]
#	dtype = dtype + [('xInd', int), ('yInd', int)]
#	data = np.array(data, dtype = dtype)

	# making the heatmap for betweenness centrality
	hmpData = np.zeros((numBins, numBins))
	for xInd in range(1,numBins + 1):
		for yInd in range(1, numBins + 1):
			binMatches = data[np.logical_and(data['xInd'] == xInd, data['yInd'] == yInd)]
			if len(binMatches) > 0:
				hmpData[xInd - 1][yInd - 1] = np.average(binMatches['bc'])
	plot = sns.heatmap(hmpData, cbar = True, square = True )
	plot.set(xlabel = '20 um', ylabel = '20 um', title = pdfNames[n] + ' bc')
	plot.figure.savefig('_'.join(['bc_patt', pdfNames[n], 'bc_plot.png']), dpi = 300)
	plt.clf()

	# creating heatmap data for power dissipated 
	data = [(np.digitize(edge[2]['x'], xBins), np.digitize(edge[2]['y'], yBins), edge[2]['power']) for edge in g.edges(data = True)]
	dtype = [('xInd', int), ('yInd', int), ('power', float)]
	data = np.array(data, dtype = dtype)

	hmpData = np.zeros((numBins, numBins))
	for xInd in range(1,numBins + 1):
		for yInd in range(1, numBins + 1):
			binMatches = data[np.logical_and(data['xInd'] == xInd, data['yInd'] == yInd)]
			if len(binMatches) > 0:
				hmpData[xInd - 1][yInd - 1] = np.log(np.average(binMatches['power']))
	plot = sns.heatmap(hmpData, cbar = True, square = True )
	plot.set(xlabel = '20 um', ylabel = '20 um', title = pdfNames[n]+ ' ln(power)')
	plot.figure.savefig('_'.join(['bc_patt', pdfNames[n], 'power_plot.png']), dpi = 300)
	plt.clf()


# zipper pattern
#pdf = lambda x,y : 1 - np.heaviside(n * nwLength / 2 - np.abs(x - width/3 + 1/3*(y-height)),0) * np.heaviside(-np.sin(2 * np.pi * nSquares * y / height),0)
#g = nwm.NanowireMesh(width = width, height = height, pdf = pdf)
#bc = nx.betweenness_centrality(g, k = int(5 * np.log(len(g.nodes))), normalized = True, weight = 'resistance')
#attrs = {node : {'bc' : bc[node]} for node in g.nodes}
#nx.set_node_attributes(g, attrs)
#g.to_netlist(netlistName = 'netlist')
#os.system('xyce netlist')
#g.update_with_xyce_output('netlist.csv')
#g.to_pickle('bc_patt_zipper.p')
#
## vertical stripes (light heavy light heavy light, then there is heavy light heavy light heavy)
##pdf = lambda x,y : 0.5*(1 +  np.heaviside(-np.sin(2 * np.pi * x/(2/5 * width)),0))
#g = nwm.NanowireMesh(width = width, height = height, pdf = pdf)
#bc = nx.betweenness_centrality(g, k = int(5 * np.log(len(g.nodes))), normalized = True, weight = 'resistance')
#attrs = {node : {'bc' : bc[node]} for node in g.nodes}
#nx.set_node_attributes(g, attrs)
#g.to_netlist(netlistName = 'netlist')
#os.system('xyce netlist')
#g.update_with_xyce_output('netlist.csv')
#g.to_pickle('bc_patt_vert_stripes_1.p')
#
##pdf = lambda x,y : 1 - 0.5 * np.heaviside(-np.sin(2 * np.pi * x/(2/5 * width)),0)
#g = nwm.NanowireMesh(width = width, height = height, pdf = pdf)
#bc = nx.betweenness_centrality(g, k = int(5 * np.log(len(g.nodes))), normalized = True, weight = 'resistance')
#attrs = {node : {'bc' : bc[node]} for node in g.nodes}
#g.to_netlist(netlistName = 'netlist')
#os.system('xyce netlist')
#g.update_with_xyce_output('netlist.csv')
#nx.set_node_attributes(g, attrs)
#
#g.to_pickle('bc_patt_vert_stripes_2.p')
#
## horizontal stripes (light heavy light heavy light, then there is heavy light heavy light heavy)
##pdf = lambda x,y : 0.5*(1 + np.heaviside(-np.sin(2 * np.pi * y/(2/5 * height)),0))
#g = nwm.NanowireMesh(width = width, height = height, pdf = pdf)
#bc = nx.betweenness_centrality(g, k = int(5 * np.log(len(g.nodes))), normalized = True, weight = 'resistance')
#attrs = {node : {'bc' : bc[node]} for node in g.nodes}
#nx.set_node_attributes(g, attrs)
#g.to_pickle('bc_patt_horiz_stripes_1.p')
#
##pdf = lambda x,y : 1 - 0.5 * np.heaviside(-np.sin(2 * np.pi * y/(2/5 * height)),0)
#g = nwm.NanowireMesh(width = width, height = height, pdf = pdf)
#bc = nx.betweenness_centrality(g, k = int(5 * np.log(len(g.nodes))), normalized = True, weight = 'resistance')
#attrs = {node : {'bc' : bc[node]} for node in g.nodes}
#nx.set_node_attributes(g, attrs)
#g.to_pickle('bc_patt_horiz_stripes_2.p')
#
## circles of high density and low density
#radius = width / 5
##pdf = lambda x,y : 0.5* (1 + np.heaviside(radius**2 - (x - width/2)**2 - (y - height/2)**2,0))
#g = nwm.NanowireMesh(width = width, height = height, pdf = pdf)
#bc = nx.betweenness_centrality(g, k = int(5 * np.log(len(g.nodes))), normalized = True, weight = 'resistance')
#attrs = {node : {'bc' : bc[node]} for node in g.nodes}
#nx.set_node_attributes(g, attrs)
#g.to_pickle('bc_patt_circle_extra.p')
#
##pdf = lambda x,y : 1 - 0.5 *  np.heaviside(radius**2 - (x - width/2)**2 - (y - height/2)**2,0)
#g = nwm.NanowireMesh(width = width, height = height, pdf = pdf)
#bc = nx.betweenness_centrality(g, k = int(5 * np.log(len(g.nodes))), normalized = True, weight = 'resistance')
#attrs = {node : {'bc' : bc[node]} for node in g.nodes}
#nx.set_node_attributes(g, attrs)
#g.to_pickle('bc_patt_circle_hole.p')
#
