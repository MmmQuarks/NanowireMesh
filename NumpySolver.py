import NanowireMesh as nwm
import numpy as np
from pathlib import Path
import sys
import argparse
import networkx as nx
from copy import deepcopy

def solve(
	g, 
	current = None, 
	voltage = None,
	voltage_label = 'voltage', current_label = 'current', power_label = 'power',
	inplace = False
	):
	'''
	Solve circuit using numpy.linalg.solve(). Dangling ends and nodes that are not in the percolating cluster are treated as having undefined voltage.
	
	args:
		g = NanowireMesh object
		specify only one of current or voltage
			current = current injected to top electrode and extracted from bottom electrode in Amperes
			voltage = voltage difference between top and bottom electrodes. bottom is 0 by convention
		voltage_label = name of attribute in node dict where solve voltage is stored
		current_label = name of attribute in edge dict where solved current is stored
		power_label = name of attribute in edge dict where solved power is stored
		inplace = whether to modify the input graph in place or return a copy with new attribute values
	returns:
		if inplace = True:
			returns None, modified g in place to have new attribute values
		if inplace = False:
			returns deep copy of g with new attribute values
	'''
	both_none = current is None and voltage is None
	both_specified = current is not None and voltage is not None
	if both_none or both_specified:
		raise ValueError(
			'''Exactly one of current of voltage must be specified: inputs are\n 
			current = {}\n'
			voltage = {}'''.format(
				current,
				voltage
			)
		)

	if not inplace:
		g = deepcopy(g)
		
	if current is None:
		current = 1
	adm, I, nodelist = nwm.linear_system(g, current = current)

	# solving for voltages
	v = np.linalg.solve(adm.todense(), I)
	v = v - v.min()
	
	# setting sheetResistance attribute (we can do this before we scale voltages and currents)
	g.sheetResistance = (v.max() - v.min()) / current
	
	# scaling voltages if voltage was specified at the beginning
	if voltage is not None:
		v *= voltage/(v.max() - v.min())

	# assigning attributes to graph elements that are in nodelist
	for (node, voltage_value) in zip(nodelist, v):
		g.nodes[node][voltage_label] = voltage_value
	for edge in g.edges:
		n1, n2 = edge
		if voltage_label in g.nodes[n1].keys() and voltage_label in g.nodes[n2].keys():
			edge_voltage = abs(g.nodes[n1][voltage_label] - g.nodes[n2][voltage_label])
			edge_resistance = g.edges[edge]['resistance']
			edge_current = edge_voltage / edge_resistance
			g.edges[edge][current_label] = edge_current
			g.edges[edge][power_label] = edge_current**2 * edge_resistance
	


	if not inplace:
		return g






if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		prog = 'NumpySolver',
		description = 'Solve nanowire mesh with numpy'
		)
	parser.add_argument(
		'path',
		help = 'Path to pickle of nanowiremesh object'
		)
	parser.add_argument(
		'current',
		help = 'Current [A] to inject into topElectrode and extract from bottomElectrode'
		)
	args = parser.parse_args(sys.argv[1:])
	g = nwm.NanowireMesh(
		makeEmpty = False,
		inPickle = args.path
		)
	solve(g, current = args.current)
	
