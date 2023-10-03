import NanowireMesh as nwm
import numpy as np
import seaborn as sns
import networkx as nx
import itertools
from Timer import Timer
from tqdm import tqdm
from matplotlib import pyplot as plt
import os
system = os.system

from dask.distributed import Client, progress

if __name__ == '__main__':

	# file path settings 
	imgFolder = 'betweenness_centrality_plots/100 50umx50um nets, percMultiple=2, 3um wires, no electrodes in bc'
	
	basePath = '/Users/adamtrebach/Documents/Research/TPV/Nanowires/centrality_networks/variable_wire_lengths'
	folders = os.listdir(basePath)
	folders = set(folders) -set(['.DS_Store'])
	folders = sorted(list(folders))

	# settings for which actions the script should perform (sometimes I will be making plots for data that has already been computed)
	shouldCalculateElectrical = False 
	shouldCalculateBC = False
	shouldMakePlots = True

	# making the parallel workers, but only if we are calculating BCs
	if shouldCalculateBC:
		client = Client(threads_per_worker = 2, n_workers = 2)

	# settings for plots
	width = 100
	height = 100
	
	numBins = 50 # number of bins per side
	xBins = np.linspace(0-.01,width + .01, numBins + 1)
	yBins = np.linspace(0 -.01, height + .01, numBins + 1)


	
		
	for folder in folders:
		# making data containers for plots	
		nodesX = []
		nodesY = []
		edgesX = []
		edgesY = []
		bc_unweighted_with_electrodes_list = []
		bc_weighted_with_electrodes_list = []
		bc_unweighted_no_electrodes_list = []
		bc_weighted_no_electrodes_list = []
		powerList = []
		networks = sorted(os.listdir( '/'.join([basePath, folder])))
		for n, fileName in enumerate(networks):
			#initializing the network
			print('Opening', fileName, 'in folder', folder)
			filePath = '/'.join([basePath, folder, fileName])
			g = nwm.NanowireMesh(inPickle = filePath)
			
			# sending to xyce and updating with results if told to do so
			if shouldCalculateElectrical:
				g.to_netlist(voltage = 1, netlistName = 'netlist')
				system('xyce netlist')
				g.update_with_xyce_output('netlist.csv')
	
			# calculating graph theoretic properties
			if shouldCalculateBC:
				print('BC Calculation')	
				bc_unweighted_with_electrodes = client.submit(nx.betweenness_centrality, g, k = 100, normalized = True)
				bc_weighted_with_electrodes = client.submit(nx.betweenness_centrality, g, k = 100, weight = 'resistance', normalized = True)
	
				# making a subgraph that excludes the top and bottom electrodes
				# and calculating centralities for it
				filter_node_fxn = lambda node : node not in [g.topElectrode, g.bottomElectrode]
				gNoElectrodes = nx.subgraph_view(g, filter_node = filter_node_fxn)
				bc_unweighted_no_electrodes = client.submit(nx.betweenness_centrality, gNoElectrodes, k = 100, normalized = True)
				bc_weighted_no_electrodes = client.submit(nx.betweenness_centrality, gNoElectrodes, k = 100, weight = 'resistance', normalized = True)
	
				# making sure that all computations have completed before we add them to the graph
				allDone = False
				futureList = [bc_unweighted_with_electrodes,
						bc_weighted_with_electrodes,
						bc_unweighted_no_electrodes,
						bc_weighted_no_electrodes]
				while not allDone:
					isDone = [future.status == 'finished' for future in futureList]
					allDone = all(isDone)
	
				# storing centralities in nodes	
				nx.set_node_attributes(g, bc_unweighted_with_electrodes.result(), 'bc_unweighted_with_electrodes')
				nx.set_node_attributes(g, bc_weighted_with_electrodes.result(), 'bc_weighted_with_electrodes') 
				nx.set_node_attributes(g, bc_unweighted_no_electrodes.result(), 'bc_unweighted_no_electrodes')
				nx.set_node_attributes(g, bc_weighted_no_electrodes.result(), 'bc_weighted_no_electrodes')
	
				# saving pickle of graph
				g.to_pickle(filePath)
				print('saved to file', filePath)
	
			# block to make (or not make) heatmap data
			if shouldMakePlots:

				# making lists of datapoints
				# making list of nodes that lie on the sample so no external nodes are included
				# note that we also filter out the electrodes since they are really not the same as other wires
				nodeList = [node for node in g.nodes if 0.01 <= g.nodes[node]['x'] <= g.width + 0.01  and 0.01 <= g.nodes[node]['y'] <= g.height + 0.01]
				nodeList = set(nodeList) - {g.topElectrode, g.bottomElectrode}
				nodesX = nodesX + [np.digitize(g.nodes[node]['x'], xBins) - 1 for node in nodeList]
				nodesY = nodesY + [np.digitize(g.nodes[node]['y'], yBins) - 1 for node in nodeList]
				bc_unweighted_with_electrodes_list = bc_unweighted_with_electrodes_list + [g.nodes[node]['bc_unweighted_with_electrodes'] for node in nodeList]
				bc_weighted_with_electrodes_list = bc_weighted_with_electrodes_list + [g.nodes[node]['bc_weighted_with_electrodes'] for node in nodeList]
				bc_unweighted_no_electrodes_list = bc_unweighted_no_electrodes_list + [g.nodes[node]['bc_unweighted_no_electrodes'] for node in nodeList]
				bc_weighted_no_electrodes_list = bc_weighted_no_electrodes_list + [g.nodes[node]['bc_weighted_no_electrodes'] for node in nodeList]

				#making list of edges that lie on the sample no external edges are included
				edgeList = [edge for edge in g.edges if 0.01 <= g.edges[edge]['x'] <= g.width + 0.01 and 0.01 <= g.edges[edge]['y'] <= g.height + 0.01]
				edgesX = edgesX + [np.digitize(g.edges[edge]['x'], xBins) - 1 for edge in edgeList]
				edgesY = edgesY + [np.digitize(g.edges[edge]['y'], yBins) - 1 for edge in edgeList]
				powerList = powerList + [g.edges[edge]['power'] for edge in edgeList]

	
		# making data containers for heatmaps (this happens in the folders block)	
		if shouldMakePlots:
		 	# making a dictionary of the datasets so we can iterate over them more easily
			bcDict = {'bc_unweighted_with_electrodes' : bc_unweighted_with_electrodes_list,
			 		'bc_weighted_with_electrodes' : bc_weighted_with_electrodes_list,
					'bc_unweighted_no_electrodes' : bc_unweighted_no_electrodes_list,
					'bc_weighted_no_electrodes' : bc_weighted_no_electrodes_list }
			# converting datasets to np arrays for easier indexing
			for key in bcDict.keys():
				bcDict[key] = np.array(bcDict[key])

			#converting power list to np array
			powerList = np.array(powerList)
			# converting coordinate lists to np array for similar indexing
			nodesX = np.array(nodesX)
			nodesY = np.array(nodesY)
			edgesX = np.array(edgesX)
			edgesY = np.array(edgesY)

			# making bc plots
			for key, data in bcDict.items():
				plotData = np.zeros((numBins, numBins))
				prefix = folder + '_' + key
				for row in range(numBins):
					for col in range(numBins):
						# selecting the data corresponding to this point on the plot
						dataInds = np.logical_and( nodesX == col, nodesY == row)
						filteredData = data[dataInds]
						if len(filteredData) > 0:
							plotValue = np.average(filteredData)
						else:
							plotValue = 0
						plotData[row, col] = plotValue	
				
				# writing the gnuplot data file
				f = open(prefix + '_gnp_data.txt', 'w')
				for line in plotData:
					f.write(' '.join([str(val) for val in line]) + '\n')
				f.close()

				# writing the gnuplot script
				f = open(prefix + '_gnp_script.txt', 'w')
				f.write('set size square\n')
				f.write('set term png size 8192, 8192\n')
				f.write('set output \"'+ prefix + '_gnp_plot.png\"\n')
				f.write('set pm3d map\n')
				f.write('set title \"' + 'avg ' + prefix.replace('_',' ') + '\" font \",200\" offset 0,20\n')
				f.write('set lmargin at screen 0.05\n')
				f.write('set rmargin at screen 0.9\n')
				f.write('set xrange [0:' + str(numBins - 1) + ']\n')
				f.write('set yrange [0:' + str(numBins - 1) + ']\n')
				xInterval = str(round( width / numBins, 2))
				yInterval = str(round( height / numBins, 2))
				f.write('set xlabel \"Units of ' + xInterval + ' um\" font \",120\" offset 0,-30\n')
				f.write('set ylabel \"Units of ' + yInterval + ' um\" font \",120\" offset -30,0\n')
				f.write('set xtics font \",100\" offset 0,-10\n')
				f.write('set ytics font \",100\"\n')
				f.write('set cbtics font \",100\"\n')
				f.write('splot \"' + prefix + '_gnp_data.txt' + '\" matrix')
				f.close()

				system('gnuplot ' + prefix + '_gnp_script.txt')

			#making power plot
			print('making power plot')
			plotData = np.zeros((numBins, numBins))
			prefix = folder + '_power'
			for row in range(numBins):
				for col in range(numBins):
					# selecting the data corresponding to this point on the plot
					dataInds = np.logical_and( edgesX == col, edgesY == row)
					filteredData = powerList[dataInds]
					if len(filteredData) > 0:
						plotValue = np.average(filteredData)
					else:
						plotValue = 0
					plotData[row, col] = plotValue	
			
			# writing the gnuplot data file
			f = open(prefix + '_gnp_data.txt', 'w')
			for line in plotData:
				f.write(' '.join([str(val) for val in line]) + '\n')
			f.close()

			# writing the gnuplot script
			f = open(prefix + '_gnp_script.txt', 'w')
			f.write('set size square\n')
			f.write('set term png size 8192, 8192\n')
			f.write('set output \"'+ prefix + '_gnp_plot.png\"\n')
			f.write('set pm3d map\n')
			f.write('set title \"' + 'avg ' + prefix.replace('_',' ') + '\" font \",200\" offset 0,20\n')
			f.write('set lmargin at screen 0.05\n')
			f.write('set rmargin at screen 0.9\n')
			f.write('set xrange [0:' + str(numBins - 1) + ']\n')
			f.write('set yrange [0:' + str(numBins - 1) + ']\n')
			xInterval = str(round( width / numBins, 2))
			yInterval = str(round( height / numBins, 2))
			f.write('set xlabel \"Units of ' + xInterval + ' um\" font \",120\" offset 0,-30\n')
			f.write('set ylabel \"Units of ' + yInterval + ' um\" font \",120\" offset -30,0\n')
			f.write('set xtics font \",100\" offset 0,-10\n')
			f.write('set ytics font \",100\"\n')
			f.write('set cbtics font \",100\"\n')
			f.write('splot \"' + prefix + '_gnp_data.txt' + '\" matrix')
			f.close()

			system('gnuplot ' + prefix + '_gnp_script.txt')
