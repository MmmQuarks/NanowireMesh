import NanowireMesh as nwm
import networkx as nx
import sys
import numpy as np
import pandas as pd
import time
import random
from copy import deepcopy

#results folder
resultsFolder = '/home/amwt/TPV/2DNanowires/Data/junctionvsinternal001/'

# making blank dataframe
df = pd.DataFrame(columns = ['taskID', 'percMultiple', 'nwLength', 'Rw', 'Rj', 'res', 'noJuncRes', 'noIntRes'])

# setting run duration
startTime = time.perf_counter()
now = time.perf_counter()
fortySixHours = 46 * 60 * 60
almostZero = 10**(-5)

# setting parameters
params = dict(width = 150,
		height = 150,
		nwLength = 10, # will be randomized
		nwLengthSD = 0,
		nwDiam = .15,
		nwDiamSD = 0,
		percMultiple = 1.5, # will be randomized
		removeNonPercolating = False,
		addInternalResistance = True,
		useFastJunctionFinder = True,
		rcMean = 10, # will be randomized
		rcSD = 0)

#nwLengthOptions = list(np.linspace(5,50,10))
#percMultipleOptions = [1.2, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20]
rcMeanOptions = [10**(2*n) for n in range(-2,3)]
iterations = 0

while now - startTime < fortySixHours:
	params['nwLength'] = np.random.uniform(5,50) #random.sample(nwLengthOptions,1)[0]
	params['percMultiple'] = np.random.uniform(1.2, 20) #random.sample(percMultipleOptions,1)[0]
	params['rcMean'] = random.sample(rcMeanOptions, 1)[0]

	params['percMultiple'] = 2
	# make network
	g = nwm.NanowireMesh(**params)

	# get data from internal resistor and scale up to find Rw
	for edge in g.edges:
		if g.edges[edge]['resistanceType'] == 'int':
			e = g.edges[edge]
			Rw = e['resistance'] * params['nwLength'] / e['length']
			break

	# removeing dangling ends
	nonConducting = [node for node in g if g.degree(node) in [0,1]]
	while nonConducting:
		g.remove_nodes_from(nonConducting)
		nonConducting = [node for node in g if g.degree(node) in [0,1]]
 
	try:
		g.solve_circuit_using_xyce(netlistName = 'netlist' + sys.argv[1],
			xycePath = '/home/amwt/TPV/2DNanowires/installations/build_Xyce-6.11_parallel_openMPI/bin/Xyce')
	except Exception:
		print('Xyce solve failure')
		continue

	res = deepcopy(g.sheetResistance)
	print('found orig res')
	originalResistances = {edge : g.edges[edge]['resistance'] for edge in g.edges}
		

	# calculating noJuncRes
	for edge in g.edges:
		if g.edges[edge]['resistanceType'] == 'cont':
			g.edges[edge]['resistance'] = almostZero
	try:
		g.solve_circuit_using_xyce(netlistName = 'netlist' + sys.argv[1],
			xycePath = '/home/amwt/TPV/2DNanowires/installations/build_Xyce-6.11_parallel_openMPI/bin/Xyce')
	except Exception:
		print('Xyce solve failure')
		continue

	noJuncRes = deepcopy(g.sheetResistance)

	# resetting all edges
	nx.set_edge_attributes(g, originalResistances, 'resistance')	

	# calculating noIntRes
	for edge in g.edges:
		if g.edges[edge]['resistanceType'] == 'int':
			g.edges[edge]['resistance'] = almostZero
	try:
		g.solve_circuit_using_xyce(netlistName = 'netlist' + sys.argv[1],
			xycePath = '/home/amwt/TPV/2DNanowires/installations/build_Xyce-6.11_parallel_openMPI/bin/Xyce')
	except Exception:
		print('Xyce solve failure')
		continue

	noIntRes = deepcopy(g.sheetResistance)

	# making dict to add to dataframe
	datum = dict(taskID = sys.argv[1],
			percMultiple = params['percMultiple'],
			nwLength = params['nwLength'],
			Rw = Rw,
			Rj = params['rcMean'],
			res= res,
			noJuncRes = noJuncRes,
			noIntRes = noIntRes)

	df = df.append(datum, ignore_index = True)
	if iterations % 10 == 0:
		df.to_csv(path_or_buf = resultsFolder + 'task' + sys.argv[1] + '.csv')

	iterations +=1
	now = time.perf_counter()
