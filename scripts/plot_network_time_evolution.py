import NanowireMesh as nwm
import numpy as np
from matplotlib import collections as mc
import pylab as pl
import matplotlib.pyplot as plt
import networkx as nx
import os
from tqdm import tqdm
import time

def make_plot(g, outFile):
	lines = [g.nodes[node]['endpoints'] for node in g.nodes]
	colors = ['black' if node in g.percolatingCluster else 'blue' for node in g.nodes]



	# plotting line segments
	fig, ax = pl.subplots()
	ax.set_aspect(aspect = 1)

	# plotting failed junctions if they exist
	if 'failedJunctions' in dir(g):
		x = [g.failedJunctions[key]['x'] for key in g.failedJunctions.keys()]
		y = [g.failedJunctions[key]['y'] for key in g.failedJunctions.keys()]
		plt.scatter(x, y, c = 'red', s = 5)

	lc = mc.LineCollection(lines, colors = colors, linewidths = 0.25, alpha = 0.5)
	ax.add_collection(lc)
	ax.autoscale()



	ax.set_xlim(0, g.width)
	ax.set_ylim(0, g.height)
	plt.draw()
	fig.savefig(outFile + '.svg') #   savefig(outFile + '.svg', fig)




g = nwm.NanowireMesh(height = 500, width = 500, nwLength = 7, nwLengthSD = 3, percMultiple = 1.5, removeNonPercolating=False, addInternalRes= True)
folder = 'frames1'


pickleFolder = 'frames1/source_data'
g.to_pickle(pickleFolder + '/p0.p')
numFailedJunctions = 0
voltage = 1

for n in range(1,100):
	g.to_netlist('netlist', voltage = voltage)
	os.system('xyce netlist')
	g.update_with_xyce_output('netlist.csv')
	g.remove_failing_components()
	if len(g.failedJunctions) == numFailedJunctions:
		voltage = voltage * 1.05
	numFailedJunctions = len(g.failedJunctions)
	g.to_pickle(pickleFolder + '/p' + str(n) + '.p')
	print('is percolating:', g.isPercolating)
	if g.isPercolating == False:
		break

	# make_plot(g,'/'.join([folder, 'p' + str(n)]))

for n in range(100):
	g = nwm.NanowireMesh(inPickle= 'frames1/source_data/p' + str(n))
	make_plot(g, 'frames1/p' + str(n))
