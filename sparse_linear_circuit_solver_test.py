import networkx as nx
import pandas as pd
import pdb
import numpy as np
import NanowireMesh as nwm
import sparse_linear_circuit_solver as circ
from scipy import sparse, linalg, stats, optimize
from copy import deepcopy
from time import time
from numba import jit


def calc_err_by_name(g, attrName):
	errs = np.array([g.nodes[node]['voltage'] - g.nodes[node][attrName] for node in g])
	errs = np.absolute(errs)
	badElements = errs > 0.1
	badCount = badElements.sum()
	rnorm = linalg.norm(errs)
	return rnorm, badCount 


def get_x0(h, voltage):
	return {node : h.nodes[node]['y'] / h.height * voltage for node in h}

def get_node_attr_or_na(g, node, attr):
	try:
		ans = g.nodes[node][attr]
	except KeyError:
		ans = np.nan
	return ans
	


def main():
	# testing various ways of improving accuracy
	numNetworks = 10
	networkSize = 50
	voltage = 1

	methods = ['xyce', 'solve', 'spsolve']

	
	df = pd.DataFrame(columns = ['trial', 'percMultiple', 'num nodes', 'method', 'rms err', 'num errs', 'time'])
	
	# running the numba code just to compile it the first time
#	print('Initializing Numba Function')
#	g = nwm.NanowireMesh()
#	circ.solve_with_admittances(g, solver = 'lstsq_numba')
	
	percMultipleList = list(np.linspace(1.2,5,10))
	for percMultiple in percMultipleList:
		for trial in range(numNetworks):
			g = nwm.NanowireMesh(width = networkSize,
					height = networkSize,
					percMultiple = percMultiple)
	
			for m in methods:
				print('Trial', trial, 'solving with', m)
				start = time()
				if m == 'xyce':
					g.solve_circuit_using_xyce(voltage = 1)
					elapsed = time() - start
					# giving the voltages a new name
					vdict = nx.get_node_attributes(g, 'voltage')
					nx.set_node_attributes(g, vdict, 'xyce')

				# initializing the numba functions if we are on the first iteration at the first percMultiple
				if 'numba' in m and percMultiple == percMultipleList[0] and trial == 0:
					print('Compiling numba solver function:', m)
					circ.solve_with_admittances(g,
						solver = m,
						currentInjections = {g.topElectrode : 1 / g.sheetResistance, g.bottomElectrode : -1 / g.sheetResistance},
						voltageName = m)
					start = time()

				if m != 'xyce':
					circ.solve_with_admittances(g,
						solver = m,
						currentInjections = {g.topElectrode : 1 / g.sheetResistance, g.bottomElectrode : -1 / g.sheetResistance},
						voltageName = m)
					elapsed = time() - start

				rmsErr, numErrs = calc_err_by_name(g, m)
				df = df.append({'trial' : trial,
						'percMultiple' : percMultiple,
						'num nodes' : len(g),
						'method' : m,
						'rms err' : rmsErr,
						'num errs' : numErrs,
						'time' : elapsed},
						ignore_index = True)
				df.to_csv('sparse_linear_circuit_solver_benchmarking.csv')
	
	pdb.set_trace()
if __name__ == '__main__':
	main()

