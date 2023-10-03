from copy import deepcopy
import NanowireMesh as nwm
import networkx as nx
import numpy as np
import itertools
import pandas as pd
import pdb
import gnuplotlib as gp
from scipy import stats
import dask
from dask.distributed import Client

def draw_random_angle():
	pi = np.pi
	anglePossibilities = [-pi/2, -pi/4, 0, pi/4, pi/2]
	angleWeights = [1/6, 1/6, 1/3, 1/6, 1/6]
	return np.random.choice(anglePossibilities, p = angleWeights)

def main(g = None):
	def _average_edge_attribute_to_node(g, node, attribute):
		if nx.degree(g, node) > 0:
			edges = g.edges(node)
			return np.mean([g.edges[edge][attribute] for edge in edges]) / 2

		# if the degree is zero, we return nan
		return np.nan

	resultsFolder = '/Users/adamtrebach/Documents/Research/TPV/Nanowires/Results/bc_paper_angle_correlations_003/'
	# attributes to plot vs angle
	attributesToPlot = ['angle',
				'betweenness_centrality',
				'percolation_centrality',
				'power']
	
	
	#setting network properties
	
	pdf = lambda x,y : 1
	
	if g == None:
		g = nwm.NanowireMesh(width = 100,
					height = 100,
					buffer = 1,
					angleFunc = draw_random_angle,
					addInternalResistance = True)

		
		
		# ------ end network initialization
		
		# saving initial snapshot of network
#		g.to_pickle(resultsFolder + 'network.p')

		# solving electrical properties
		g.solve_circuit_using_xyce(netlistName = resultsFolder + 'netlist')
		
	
		print('Calculating Centralities')
		bc, pc = dask.compute( dask.delayed(nx.betweenness_centrality)(g, weight = 'resistance', normalized = True),
					dask.delayed(nx.percolation_centrality)(g, attribute = 'voltage', weight = 'resistance'))
		
		nx.set_node_attributes(g, bc, 'betweenness_centrality')
		nx.set_node_attributes(g, pc, 'percolation_centrality')
	
		# assign powers to nodes
		nodesWithoutElectrodes = set(g.nodes) - {g.topElectrode, g.bottomElectrode}
		powersDict = {node : np.average([g.edges[edge]['power']/2 for edge in g.edges(node)]) for node in nodesWithoutElectrodes}
		nx.set_node_attributes(g, powersDict, 'power')

	else:
		# if we are passed a non-empty g it is assumed it has all the properties calculated
		pass

	# storing the network as a pickle for later analysis if necessary
	g.to_pickle(resultsFolder + 'network.p')

	# making pandas dataframe
	nodesWithoutElectrodes = set(g.nodes) - {g.topElectrode, g.bottomElectrode}
	data = dict.fromkeys(attributesToPlot, 0)
	for attr in data.keys():
		data[attr] = np.array([g.nodes[node][attr] for node in nodesWithoutElectrodes])
	# taking absolute value of angles and rounding angles so grouping is done correctly
	for n in range(len(data['angle'])):
		val = data['angle'][n]
		val = round( abs(val), 5)
		data['angle'][n] = val
		
	df = pd.DataFrame(data = data)
	df['angle'] = df['angle'].apply(np.abs)
	df.to_csv(resultsFolder + 'dataframe.csv')


	for attribute in attributesToPlot:
		if attribute != 'angle':
			print(attribute)
			plotData = df[['angle', attribute]].dropna()
			# calculating the quantiles of each angle so we have multiple lines on plot
			plotData = plotData.groupby('angle').agg(['mean', 'std', 'count'])
			plotOptions = {'title' : 'Mean ' + attribute.replace('_', ' ') + ' as function of Angle',
				'xlabel' : 'Abs(Angle)',
				'ylabel' : attribute.replace('_', ' '),
				'terminal' : 'png font \",80\" size 5000,3000',
				'output' : resultsFolder + 'bc_paper_angle_correlations_' + attribute + '.png'}

			x = np.array(plotData.index)
			y = np.array(plotData[attribute]['mean'])
			yError = np.array(plotData[attribute]['std']) / np.sqrt(plotData[attribute]['count'])
			curveOptions = {'tuplesize' : 3,
				'with' : 'yerrorbars lw 10 pointtype 7 pointsize 10'
				}
			# making a linear fit
			slope, yIntercept, r, p, stdErr = stats.linregress(x,y)
			fitString = str(yIntercept) + ' + ' + str(slope) + ' * x lw 10 title \"r^2 = ' + str(round(r**2, 3)) + '\"'
			plotOptions['equation'] = fitString
			gp.plot((x,y,yError, curveOptions),**plotOptions)
		
if __name__ == '__main__':
	client = Client(n_workers = 2)
	main()
