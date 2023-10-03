import NanowireMesh as nwm
import networkx as nx
import gnuplotlib as gp
from scipy import stats
import numpy as np

g = nwm.NanowireMesh(width = 100, height = 100, buffer = 2)
g.solve_circuit_using_xyce()
bc = nx.edge_betweenness_centrality(g, weight = 'resistance')
x = np.array([bc[edge] for edge in g.edges])
y = np.array([g.edges[edge]['power'] for edge in g.edges])
curve0Options = {'with' : 'points pointtype 7 pointsize 4'}
curve0 = (x,y, curve0Options)

slope, yIntercept, r, p, stdErr = stats.linregress(x,y)
fitString = str(yIntercept) + '+' + str(slope) + '*x with lines lw 8 title \"r^2 = ' + str(round(r**2,3)) + '\"'

plotOptions = {'title' : 'Edge Power vs Edge Betweenness Centrality',
		'terminal' : 'png font \",50\" size 5000,3000',
		'output' : 'bc_paper_power_vs_bc_figure.png',
		'equation' : fitString,
		'xlabel' : 'Betweenness Centrality',
		'ylabel' : 'Power [W]'}
gp.plot(curve0,**plotOptions)

