import NanowireMesh as nwm
import pickle
import re
import pdb
import pandas as pd
import networkx as nx
import random
import numpy as np
from scipy import stats, interpolate, optimize
from copy import deepcopy
import ParallelCentrality as PC
import csv
import itertools
import os
import sys
import argparse
from pathlib import Path
from NanowireNetworkEvolutionFunctions import *

parser = argparse.ArgumentParser(
	description = 'Evolve Networks by Node Removal and compare to MCTS evolution'
	)

parser.add_argument('--results-folder', type = str, help = 'Path to folder where output will be stored')
parser.add_argument('--network-folder', type = str, help = 'Path to folder where starting networks are stored')
parser.add_argument('--network-file-pattern', type = str, help = 'The file pattern that identifies something as a series of states')
parser.add_argument('--centralities-to-test', type = str, nargs = '+', help = 'Which centralities to use in the simulation')
args = parser.parse_args(sys.argv[1:])


# allowing for customization without changing this script or the submit script 
resultsFolder = Path(args.results_folder)
network_folder = Path(args.network_folder)

xycePath = '/home/amwt/TPV/2DNanowires/installations/build_Xyce-6.11_parallel_openMPI/bin/Xyce'

centralities = dict(
		randCleaned = dict(name = 'random_cleaned',
			color = 'black',
			func = lambda g : {node : np.random.rand() for node in g},
			clean = True,
			pointtype = 6
			),
		lowBC = dict(name = 'betweenness_centrality',
				color = 'dark-violet',
				func = lambda g : nx.betweenness_centrality(g, weight = 'resistance'),
				clean = True,
				pointtype = 8
				)
		)

evolveData = []
possible_network_files = list(network_folder.glob('*'))
network_files = [f for f in possible_network_files if re.search(args.network_file_pattern, str(f))]
network_files = sorted(network_files)
# hard coding a few things
nwDiam = 0.15
for origNetwork in network_files: 
	# making sure that 
	# opening network for evolution
	with open(origNetwork, 'rb') as file:
		states = pickle.load(file)

	# getting the trial number
	trial = re.search(
		'(\d+)(\.p)', # the search pattern
		str(origNetwork) # the string to search
		).group(1) # showing that we want to return the item from the first group in parens
		#0 index returns the whole matched string
	
	# iterating through all states and getting evolution data
	for state in states:
		if state.isTerminal():
			break # we don't want to record the data for the terminal states since their resistance is infinite

		g = state.g
		# a kludge because the networks attributes are only properly stored in the graph of the first state
		if state == states[0]:
			nwDiam = g.nwDiam
			nwLength = g.nwLength
			width = g.width
			height = g.height
		g.solve_circuit_using_xyce(
			xycePath = xycePath,
			netlistName = str(
				Path(
					resultsFolder,
					'netlist_node_removal'
					)
				)
			)
		evolveData.append(
			dict(
				trial = trial,
				centrality = 'MCTS',
				percMultiple = get_perc_multiple(g),
				resistance = g.sheetResistance,
				nwLength = nwLength,
				nwDiam = nwDiam,
				width = width,
				height = height,
				lengthUnits = 'um'
				)
			)
			

	# getting the graph from the original network (which is taken from the first state)
	h = states[0].g
	# iterating through the betweenness based centrality methods
	for key in centralities.keys():
		g = deepcopy(h)
		counter = 0
		isPercolating = nx.has_path(g, g.topElectrode, g.bottomElectrode)
		evolveData.append(
			dict(
				trial = trial,
				centrality = centralities[key]['name'],
				percMultiple = get_perc_multiple(g),
				resistance = g.sheetResistance,
				nwLength = g.nwLength,
				nwDiam = g.nwDiam,
				width = g.width,
				height = g.height,
				lengthUnits = 'um'
			)
		)
		# calculate the given centrality if it's not already calculated
		for node in g:
			centrality_name = centralities[key]['name']
			try:
				g.nodes[node][centrality_name]
			except KeyError:
				values = centralities[key]['func'](g)
				nx.set_node_attributes(g, values = values, name = centrality_name)
				break

		while isPercolating:
			# return sorted dict of centrality : node pairs
			# sorted in increasing order of centrality
			sd = sort_nodes_by_attribute(g,
							attribute = centralities[key]['name'],
							includeElectrodes = False)
			c, node = sd.peekitem(index = 0)
			nodesToRemove = get_wire_segments(g, node) 

			# removing nodes from graph
			poppedNodes, poppedEdges = g.pop_nodes_with_edges(nodesToRemove)
			isPercolating = nx.has_path(g, g.topElectrode, g.bottomElectrode)
			
			# saving data periodically or if this is the last iteration
			if not isPercolating or (counter > 0 and counter % 5 == 0):
				# adding back in the last removed nodes so the circuit can be solved
				# if we have just destroyed conductivity
				if not isPercolating:
					g.add_nodes_from(poppedNodes)
					g.add_edges_from(poppedEdges)
				if centralities[key]['clean']:
					clean_network(g)
				# solving the circuit and adding the data to our list
				try:
					g.solve_circuit_using_xyce(
						xycePath = xycePath,
						netlistName = str(
							Path(
								resultsFolder,
								'netlist_node_removal'
								)
							)
						)
				except csv.Error as e:
					print(e)

				evolveData.append(
					dict(
						trial = trial,
						centrality = centralities[key]['name'],
						percMultiple = get_perc_multiple(g),
						resistance = g.sheetResistance,
						nwLength = g.nwLength,
						nwDiam = g.nwDiam,
						width = g.width,
						height = g.height,
						lengthUnits = 'um'

					)
				)

				# save to pickle
				picklePath = Path(
					resultsFolder,
					'_'.join([
						'trial', str(trial).zfill(2),
						'cent' , centralities[key]['name'],
						'iter' , str(counter).zfill(5) + '.p'
						])
					)

				g.to_pickle(outFileName = str(picklePath))

				# saving data from evolutions to CSV
				# this should be done frequently so if it crashes we don't lose everything
				pd.DataFrame(evolveData).to_csv(path_or_buf = Path(resultsFolder, 'evolveData.csv'), index = False)
				# breaking the while loop if the sample is not percolating
				if not isPercolating:
					break
		
			counter +=1
	

evolveData = pd.DataFrame(evolveData)
# calculating extinction coefficients for all present nanowire geometries
coefficients = evolveData[['nwDiam']].drop_duplicates()
coefficients['C_ext_per_length'] = coefficients.apply(
	lambda x : nwm.extinction_coefficient(
		radius = x['nwDiam'] / 2 * 1E-6,
		mr = 0.055 + 3.32j,
		wavelength = 550E-9
	),
	axis = 1
)
evolveData = evolveData.merge(
	coefficients,
	on = 'nwDiam',
	how = 'left'
	)
evolveData['C_ext_per_wire'] = evolveData['C_ext_per_length'] * evolveData['nwLength'] * 1E-6
# getting number density
evolveData['n_s'] = evolveData.apply(
	lambda x : effective_perc_multiple_to_density(
		x['percMultiple'],
		x['nwLength'] * 1E-6,
		x['height'] * 1E-6,
		x['width'] * 1E-6
		),
	axis = 1
	)
evolveData['transparency'] = evolveData.apply(
	lambda x : np.exp(-x['n_s'] * x['C_ext_per_wire']),
	axis = 1
	)
evolveData['fom'] = evolveData.apply(
	lambda x : fom(x['resistance'], x['transparency']),
	axis = 1
	)
evolveData.to_csv(
	Path(resultsFolder, 'evolveData.csv'),
	index = False
	)
