import NanowireMesh as nwm
from copy import deepcopy
import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
import subprocess
pi = np.pi

# options
width = 10
height = 10
addInternalResistance = True
nwDiam = 0.15
rcMean = 10
rcSD = 0




g = nwm.NanowireMesh(width = width,
			height = height,
			addInternalResistance = addInternalResistance,
			nwDiam = nwDiam,
			rcMean = 10,
			rcSD = 0,
			disableTQDM = True)

nodesToRemove = set(g.nodes) - {g.topElectrode, g.bottomElectrode}
g.remove_nodes_from(nodesToRemove)

# the nodes will be added according to the pattern layed out here: https://www.desmos.com/calculator/eugmudqvqb 
# node 10 is top electrode and node 0 is bottom.
# wires are numbered 1-7 in the order of equations shown on desmos


# there are 9 wires to be added so top and bottom electrodes should be numbered 0 and 8 respectively
mapping = {g.topElectrode : 10, g.bottomElectrode : 0}
# note that copy = False is vital otherwise this command appears to have a bug in it
nx.relabel_nodes(g, mapping, copy = False)
g.topElectrode = 10
g.bottomElectrode = 0

# adding node endpoints
endpoints = {1 : [(0,4), (4.2, 10.3)],
		2: [(4.4, 3), (7.25, 10.125)],
		3: [(4, 8), (0, 5.2)],
		4: [(8, 4), (5, 8.5)],
		5: [(2.625, 7.5), (2.95, 1)],
		6: [(5, 3), (10, 5.5)],
		7: [(8.3, 0), (9.4, 5.5)],
		8: [(1, 2.2), (7, -.2)],
		9: [(0.5, 5), (9.5, 3.2)]
	}


# adding voltages as calculated from this circuit simulator
# https://www.falstad.com/circuit/circuitjs.html?cct=$+1+0.000005+10.20027730826997+50+5+43%0Av+640+480+-32+480+0+0+40+1+0+0+0.5%0Aw+-32+480+-64+256+0%0Aw+640+304+656+256+0%0Aw+-64+256+-32+192+2%0Ar+-32+192+64+144+0+10%0Aw+64+144+96+128+2%0Ar+96+128+208+16+0+10%0Aw+208+16+240+16+2%0Ar+240+16+352+96+0+10%0Aw+352+96+384+96+2%0Ar+384+96+560+128+0+10%0Aw+560+128+592+144+2%0Ar+592+144+656+256+0+10%0Ax+-96+216+-70+219+4+24+10%0Ax+676+301+689+304+4+24+0%0Ax+52+121+65+124+4+24+1%0Ax+216+-23+229+-20+4+24+3%0Ax+364+59+377+62+4+24+5%0Ax+588+110+601+113+4+24+8%0Aw+224+208+272+208+2%0Ax+239+168+252+171+4+24+9%0Aw+640+304+640+480+0%0Ar+96+128+224+208+0+10%0Ar+272+208+384+96+0+10%0Aw+384+320+432+320+2%0Ax+406+271+419+274+4+24+7%0Ar+384+320+272+208+0+10%0Ar+640+304+432+320+0+10%0Aw+240+400+288+400+2%0Ar+288+400+384+320+0+10%0Ax+257+373+270+376+4+24+6%0Ar+240+400+224+208+0+10%0Aw+64+416+112+416+2%0Ar+112+416+240+400+0+10%0Ax+84+385+97+388+4+24+4%0Aw+0+304+48+304+2%0Ar+64+416+48+304+0+10%0Ax+17+275+30+278+4+24+2%0Ar+48+304+224+208+0+10%0Ar+0+304+-64+256+0+10%0A
voltages = {0 : 0,
		1 : 0.708,
		2 : 0.718,
		3 : 0.579,
		4 : 0.610,
		5 : 0.450,
		6 : 0.501,
		7 : 0.349,
		8 : 0.225,
		9 : 0.545,
		10 : 1 }

# adding nodes tograph
g.add_nodes_from(list(endpoints.keys()))
# adding endpoints
nx.set_node_attributes(g, endpoints, 'endpoints')
# calculating and adding properties to nodes
attrs = {}
for node in g.nodes:
	if node not in [g.topElectrode, g.bottomElectrode]:
		p0, p1 = g.nodes[node]['endpoints']
		x0, y0 = p0
		x1, y1 = p1
		attrs[node] = {'x' :  (x1 + x0) / 2,
				'y' : (y1 + y0) / 2,
				'diam' : 0.15,
				'length' : np.sqrt( (y1 - y0)**2 + (x1 - x0)**2),
				'voltage' : voltages[node],
				'temp' : g.initialTemp,
				'angle' : np.arctan( (y1 - y0) / (x1 - x0) ) }
		# adding wire mass which can be done most easily once the other attrs are there
		diam = attrs[node]['diam']
		length = attrs[node]['length']
		mass = g.calc_wire_mass(diam, length)
		attrs[node]['mass'] = mass
#adding voltages for electrodes manually
attrs[10] = {'voltage' : 1}
attrs[0] = {'voltage' : 0}

# writing properties to graph
nx.set_node_attributes(g, attrs)


#making the new graph that will run the algorithms and check against the answers known to be correct
h = deepcopy(g)

# calculating the properties of h using the initialization functions

# adding junctions
h.fast_add_junctions(disableTQDM = True)
h.isPercolating = nx.has_path(h, h.topElectrode, h.bottomElectrode)
if not h.isPercolating:
	print('Error: NanowireMesh.py thinks that this sample network is not percolating')
h.find_percolating_cluster()




# note there should be 16 contact edges
# manually adding the points of intersection and edges
edgeAttrs = { (1, 10) : {'x' : 4,
			'y' : 10,
			'resistance' : rcMean,
			'resistanceType' : 'cont',
			'mass' : 0,
			'temp' : g.initialTemp,
			'power' : 0,
			'length' : 'N/A',
			'diam' : 'N/A' },
		(1, 3) : {'x' : 1.5,
			'y' : 6.25,
			'resistance' : rcMean,
			'resistanceType' : 'cont',
			'mass' : 0,
			'temp' : g.initialTemp,
			'power' : 0,
			'length' : 'N/A',
			'diam' : 'N/A' },
		(1, 9) : {'x' : 0.647,
			'y' : 4.971,
			'resistance' : rcMean,
			'resistanceType' : 'cont',
			'mass' : 0,
			'temp' : g.initialTemp,
			'power' : 0,
			'length' : 'N/A',
			'diam' : 'N/A' },
		(2, 10) : {'x' : 7.2,
			'y' : 10,
			'resistance' : rcMean,
			'resistanceType' : 'cont',
			'mass' : 0,
			'temp' : g.initialTemp,
			'power' : 0,
			'length' : 'N/A',
			'diam' : 'N/A' },
		(2, 4) : {'x' : 6,
			'y' : 7,
			'resistance' : rcMean,
			'resistanceType' : 'cont',
			'mass' : 0,
			'temp' : g.initialTemp,
			'power' : 0,
			'length' : 'N/A',
			'diam' : 'N/A' },
		(2, 9) : {'x' : 4.852,
			'y' : 4.13,
			'resistance' : rcMean,
			'resistanceType' : 'cont',
			'mass' : 0,
			'temp' : g.initialTemp,
			'power' : 0,
			'length' : 'N/A',
			'diam' : 'N/A' },
		(3, 5) : {'x' : 2.647,
			'y' : 7.053,
			'resistance' : rcMean,
			'resistanceType' : 'cont',
			'mass' : 0,
			'temp' : g.initialTemp,
			'power' : 0,
			'length' : 'N/A',
			'diam' : 'N/A' },
		(4, 6) : {'x' : 7.75,
			'y' : 4.375,
			'resistance' : rcMean,
			'resistanceType' : 'cont',
			'mass' : 0,
			'temp' : g.initialTemp,
			'power' : 0,
			'length' : 'N/A',
			'diam' : 'N/A' },
		(5, 9) : {'x' : 2.773,
			'y' : 4.545,
			'resistance' : rcMean,
			'resistanceType' : 'cont',
			'mass' : 0,
			'temp' : g.initialTemp,
			'power' : 0,
			'length' : 'N/A',
			'diam' : 'N/A' },
		(5, 8) : {'x' : 2.929,
			'y' : 1.429,
			'resistance' : rcMean,
			'resistanceType' : 'cont',
			'mass' : 0,
			'temp' : g.initialTemp,
			'power' : 0,
			'length' : 'N/A',
			'diam' : 'N/A' },
		(6, 7) : {'x' : 9.333,
			'y' : 5.167,
			'resistance' : rcMean,
			'resistanceType' : 'cont',
			'mass' : 0,
			'temp' : g.initialTemp,
			'power' : 0,
			'length' : 'N/A',
			'diam' : 'N/A' },
		(6, 9) : {'x' : 6.571,
			'y' : 3.786,
			'resistance' : rcMean,
			'resistanceType' : 'cont',
			'mass' : 0,
			'temp' : g.initialTemp,
			'power' : 0,
			'length' : 'N/A',
			'diam' : 'N/A' },
		(7, 9) : {'x' : 8.962,
			'y' : 3.308,
			'resistance' : rcMean,
			'resistanceType' : 'cont',
			'mass' : 0,
			'temp' : g.initialTemp,
			'power' : 0,
			'length' : 'N/A',
			'diam' : 'N/A' },
		(7, 0) : {'x' : 8.3,
			'y' : 0,
			'resistance' : rcMean,
			'resistanceType' : 'cont',
			'mass' : 0,
			'temp' : g.initialTemp,
			'power' : 0,
			'length' : 'N/A',
			'diam' : 'N/A' },
		(8, 0) : {'x' : 6.5,
			'y' : 0,
			'resistance' : rcMean,
			'resistanceType' : 'cont',
			'mass' : 0,
			'temp' : g.initialTemp,
			'power' : 0,
			'length' : 'N/A',
			'diam' : 'N/A' },
	}

# adding edge names
for edge in edgeAttrs.keys():
	name = '_'.join([ 'rcont', str(int(edge[0])), str(int(edge[1]))])
	edgeAttrs[edge]['name'] = name

g.add_edges_from(list(edgeAttrs.keys()))
nx.set_edge_attributes(g, edgeAttrs)
g.find_percolating_cluster()

print('Comparing Intersections in manually generated network and automatically generated one.')
noErrors = True
for edge in h.edges:
	hData = h.edges[edge]
	try:
		gData = g.edges[edge]
		xMatch = np.abs(gData['x'] - hData['x']) <= 0.01
		yMatch = np.abs(gData['y'] - hData['y']) <= 0.01
		noErrors = noErrors and xMatch and yMatch
		if xMatch and yMatch:
			pass
		else:
			print('ERROR: Coordinates of junction between', edge, 'off by more than 0.01.') 
	except KeyError:
		print('ERROR: The junction between nodes', edge, 'exists in the automatically created graph but not the manual one.')	
		noErrors = False

for edge in g.edges:
	gData = g.edges[edge]
	try:
		hData = h.edges[edge]
	except KeyError:
		print('ERROR: The junction between nodes', edge, 'exists in the manually created graph but not the automatic one.') 
		noErrors = False

if noErrors:
	print('No Errors found in comparison of intersections.')


print('Visualizing Network so Circuit can be constructed')
ax = plt.axes()
nx.draw_networkx(g, with_labels = True, ax = ax)
plt.savefig('validation_equiv_network.png', dpi = 300)
subprocess.check_call(['open', 'validation_equiv_network.png'])

print('Comparing node voltages before adding internal resistance')


# comparing voltages at nodes
h.solve_circuit_using_xyce()

voltagesMatch = True
for node in g.nodes:
	match = np.abs(g.nodes[node]['voltage'] - h.nodes[node]['voltage']) <= 0.01
	voltagesMatch = voltagesMatch and match
	if not match:
		print('Voltage mismatch found on node', node)
		noErrors = False

if voltagesMatch:
	print('All node voltages Match to within 0.01 before adding internal resistance')
else:
	print('ERROR: Not all node voltages match')


print('Internal Resistance Testing not yet Implemented')
# add internal resistance
#if addInternalResistance:
#	h.add_internal_resistance()
#
#h.find_percolating_cluster()


