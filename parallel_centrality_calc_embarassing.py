import NanowireMesh as nwm
import networkx as nx
import os
from mpi4py import MPI
import time


basePath = '/home/amwt/TPV/2DNanowires/Data/centrality_networks'
folderList = os.listdir(basePath)
folderList.sort(reverse = True)
comm = MPI.COMM_WORLD
procNum = comm.Get_rank() 
xycePath = '/home/amwt/TPV/2DNanowires/installations/build_Xyce-6.11_parallel_openMPI/bin/Xyce'
print('Processor', procNum, ':', 'up and running')
time.sleep(5)
for folder in folderList:
	print('Processor', procNum,':', 'starting work on folder', folder)
	bottomNumber = procNum * 10
	topNumber = (procNum + 1) * 10
	for fileNumber in range(bottomNumber,topNumber):
		# Loading network from pickled file
		fileName = ''.join(['net', str(fileNumber).zfill(3),'.p'])
		print('Processor', procNum, ':', 'working on', fileName)
		time.sleep(5)
		g = nwm.NanowireMesh(inPickle = '/'.join([basePath, folder, fileName]))
		
		time.sleep(9)
			
		# Calculating Electrical Properties
		netlistName = fileName.replace('.p','')
		g.to_netlist(netlistName = '/'.join([basePath, folder, netlistName]))
		os.system(' '.join([xycePath, '/'.join([basePath, folder, netlistName])]))
		nwm.update_with_xyce_output(g, '/'.join([basePath, folder, netlistName + '.csv']))
		
		# calculating centrality properties
		print('Processor' , procNum, ':', 'calculating centralities')
		betweennessCentrality = nx.betweenness_centrality(g, k = 500, normalized = True, weight = 'resistance')
		percolationCentralityY = nx.percolation_centrality(g, attribute = 'y', weight = 'resistance')
		percolationCentralityVoltage = nx.percolation_centrality(g, attribute = 'voltage', weight = 'resistance')
		
		# writing to network and saving
		print('Processor' , procNum, ':', 'updating network with xyce results and saving')
		attributes = {node : {'betweennessCentrality' : betweennessCentrality[node], 'percolationCentralityY' : percolationCentralityY[node], 'percolationCentralityVoltage' : percolationCentralityVoltage[node]} for node in g.nodes}
		nx.set_node_attributes(g, attributes)
		nwm.to_pickle(g, outFileName = '/'.join([basePath, folder, fileName]))
		h = nwm.NanowireMesh(inPickle = '/'.join([basePath, folder, fileName]))
		print('Processor' , procNum, ':', 'completed work on file', fileName)
	print('Processor', procNum, ':', 'completed work in folder', folder)

print('Processor', procNum, ':', 'all work completed')
