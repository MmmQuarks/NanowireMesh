import NanowireMesh as nwm
import plotly.graph_objects as go
import numpy as np
import networkx as nx
from scipy import interpolate as interp
from tqdm import tqdm


#interpolating manually
xBins = [[n, n+1] for n in range(50)]
yBins = [[n, n+1] for n in range(50)]

z = [ [ [] for m in range(len(xBins))] for n in range(len(yBins))]

for trial in tqdm(range(100)):
	g = nwm.NanowireMesh(width = 100, height = 100, percMultiple= 1.5)

	# find all non electrode nodes
	availNodes = [node for node in g.nodes if node not in [g.topElectrode, g.bottomElectrode]]

	print('calculating betweenness centrality')
	bc = nx.betweenness_centrality(g, k = int(0.1 * len(g.nodes)), weight = 'resistance')

	data = [(node, g.nodes[node]['x'], g.nodes[node]['y'], bc[node]) for node in g.nodes]
	data = np.array(data, dtype = [('node', 'f8'), ('x', 'f8'), ('y', 'f8'), ('bc', 'f8')])

	for col, xbin in enumerate(tqdm(xBins)):
		for row, ybin in enumerate(tqdm(yBins)):
			binData = data[ xbin[0] <= data['x']] 
			binData = binData[ binData['x'] < xbin[1]]
			binData = binData[ ybin[0] <= binData['y']]
			binData = binData[ binData['y'] < ybin[1]]
			if len(binData) > 0:
				z[row][col] = z[row][col] + list(binData['bc'])

for row in range(len(z)):
	for col in range(len(z[row])):
		if len(z[row][col]) > 0:
			z[row][col] = sum(z[row][col])
		else:
			z[row][col] = 0
print(z)
print('begin plotting')

fig = go.Figure(data = go.Contour(z = z))
fig.show()